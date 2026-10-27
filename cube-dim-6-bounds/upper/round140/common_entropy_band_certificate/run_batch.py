#!/usr/bin/env python3
"""Dual-gated finite common-entropy certificates; reuse C134 corner values."""
from pathlib import Path
from fractions import Fraction as F
import hashlib,itertools,json,math,os,sys,time,traceback
sys.dont_write_bytecode=True
P=Path(__file__).resolve().parent
DEN=2**24
def read(p):return json.loads(Path(p).read_text())
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def require(ok,msg):
    if not ok:raise RuntimeError(msg)
def canon(x):return (json.dumps(x,sort_keys=True,separators=(',',':'))+'\n').encode()
def key(k,p):return hashlib.sha256(canon({'semantic_key':k,'point':p})).hexdigest()
def points(box):return [['1']+[box[i][b]for i,b in enumerate(bits)]for bits in itertools.product((0,1),repeat=5)]
def write_new(name,x):
    with (P/name).open('x')as f:json.dump(x,f,indent=2);f.write('\n');f.flush();os.fsync(f.fileno())
def release_requirements():
    return {'round':'140','authorized':True,'scope':'ONLY_C134_ELEVEN_DEPTH6_UNRESOLVED_LEAVES_COMMON_TWO_BAND_ENTROPY','target':'1/63500','precision_bits':512,'workers':1,'nice_minimum':10,'time_limit':None,'selected_leaf_count':11,'pair_count':2816,'weight_denominator':DEN,'max_float_objective_calls_per_pair':70,'max_Arb_logs':14176,'native_calls_permitted':0,'heat_kernel_calls_permitted':0,'source_point_exp_calls_permitted':0,'final_exp_calls_permitted':0,'old_paid_leaf_retries_permitted':0,'run_batch_sha256':sha(P/'run_batch.py'),'proof_sha256':sha(P/'PROOF.md'),'selected_metadata_sha256':sha(P/'selected_metadata.json'),'input_hashes_sha256':sha(P/'input_hashes.json'),'prepared_packet_sha256':sha(P/'prepared_packet.json')}
def release_gate():
    r=read(P/'parent_release.json')
    for k,v in release_requirements().items():require(r.get(k)==v,'RELEASE_MISMATCH:'+k)
    require(r.get('main_approved')is True and r.get('independent_review_approved')is True,'DUAL_GATE_REQUIRED')
    require(r.get('C134_point_cache_independently_accepted')is True,'PARENT_SOURCE_CACHE_ACCEPTANCE_REQUIRED')
def preflight():
    pins=read(P/'input_hashes.json');selection=read(P/'selected_metadata.json')
    for s in pins['frozen_sources']:require(sha(s['path'])==s['sha256'],'FROZEN_HASH:'+s['path'])
    source_dir=Path(pins['source_directory_allowlist'])
    cache_path=source_dir/'point_cache.json'
    cache=read(cache_path)['entries'];used=set()
    require(selection['target']=='1/63500'and len(selection['selected_leaves'])==11,'FIXED_SELECTION')
    cfg=read(source_dir/'fixed_roots_and_atlas.json')
    roots={r['original_node_id']:r for r in cfg['roots']};expected=set();paid=0
    for rid,root in roots.items():
        tree=read(source_dir/'roots'/(rid+'.json'))
        for nid in tree['closed_partition_leaf_ids']:
            node=tree['nodes'][nid]
            if node['status']=='LOCAL_LEAF_PASS':paid+=1
            else:
                require(node['status']=='DEPTH6_UNRESOLVED'and node['depth']==6,'ONLY_DEPTH6_UNRESOLVED')
                expected.add((rid,nid))
    require(paid==62 and expected=={(r['root_id'],r['leaf_id'])for r in selection['selected_leaves']},'EXACT_UNRESOLVED_SELECTION')
    for leaf in selection['selected_leaves']:
        root=roots[leaf['root_id']];node=read(source_dir/'roots'/(leaf['root_id']+'.json'))['nodes'][leaf['leaf_id']]
        require(leaf['box']==node['box']and leaf['source_keys']==[b['candidate_keys']for b in root['bands'][:2]]and leaf['windows']==[b['window']for b in root['bands'][:2]],'ORIGINAL_BOX_FIXED_ATLAS_BINDING')
        require([len(k)for k in leaf['source_keys']]==[16,16],'FIXED_ATLAS_SIZE')
        for keys in leaf['source_keys']:
            for k in keys:
                for point in points(leaf['box']):
                    ck=key(k,point);r=cache[ck]
                    require(r['point']==point and r['semantic_key']==k and F(0)<F(r['lower'])<=F(r['upper']),'SAVED_POSITIVE_BOUND_IDENTITY')
                    used.add(ck)
    require(len(used)==selection['unique_required_source_point_entries']==8544,'REQUIRED_SOURCE_POINTS')
    return selection,cache
def choose_weight(a,b):
    calls=0;candidates=[]
    def score(w):
        nonlocal calls
        calls+=1;require(calls<=70,'FLOAT_CALL_CAP')
        h=0.0 if w in (0.0,1.0)else-w*math.log(w)-(1-w)*math.log1p(-w)
        y=h+min(w*x+(1-w)*z for x,z in zip(a,b));require(math.isfinite(y),'FINITE_FLOAT_SELECTOR')
        candidates.append((y,w));return y
    score(0.0);score(1.0);score(0.5)
    lo,hi=0.0,1.0;ratio=(math.sqrt(5)-1)/2
    x=hi-ratio*(hi-lo);y=lo+ratio*(hi-lo);sx=score(x);sy=score(y)
    for _ in range(64):
        if sx<sy:lo=x;x=y;sx=sy;y=lo+ratio*(hi-lo);sy=score(y)
        else:hi=y;y=x;sy=sx;x=hi-ratio*(hi-lo);sx=score(x)
    score((lo+hi)/2);require(calls==70,'LOCKED_FLOAT_CALL_COUNT')
    best=max(candidates,key=lambda t:(t[0],-t[1]))
    # Exactly two endpoint objective calls avoid entropy logs; every other
    # call uses one log and one log1p. These are selector-only operations.
    return min(DEN,max(0,round(best[1]*DEN))),calls,best[0],2*(calls-2)
def run():
    release_gate();require(not(P/'started.json').exists(),'NO_AUTOMATIC_RESTART')
    os.nice(max(0,10-os.getpriority(os.PRIO_PROCESS,0)))
    require(all(os.environ.get(k)=='1'for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS']),'THREAD_ENV_ONE')
    selection,cache=preflight()
    import flint
    from flint import arb,ctx
    require(flint.__version__=='0.9.0'and flint.__FLINT_VERSION__=='3.6.0','PINNED_FLINT_VERSION')
    ctx.prec=512;target=F(1,63500);logs={};entropy={};counts={'Arb_source_logs':0,'Arb_entropy_logs':0,'Arb_corner_score_combinations':0,'float_source_logs':0,'float_entropy_logs':0,'float_objective_calls':0,'pair_certificates':0,'native_calls':0,'heat_kernel_calls':0,'source_point_exp_calls':0,'final_exp_calls':0,'old_paid_leaf_retries':0}
    def ar(q):q=F(q);return arb(q.numerator)/q.denominator
    def endpoints(v):
        require(v.is_finite(),'FINITE_ARB_SCORE');return F(str(v.lower().fmpq())),F(str(v.upper().fmpq()))
    def getlog(k,point):
        ck=key(k,point)
        if ck not in logs:
            ratio=F(cache[ck]['lower'])/target;fv=float(ratio);require(math.isfinite(fv)and fv>0,'FINITE_POSITIVE_FLOAT_RATIO')
            counts['Arb_source_logs']+=1;counts['float_source_logs']+=1
            require(counts['Arb_source_logs']<=8544,'SOURCE_LOG_CAP')
            v=ar(ratio).log();require(v.is_finite(),'FINITE_SOURCE_LOG');logs[ck]=(v,math.log(fv))
        return logs[ck]
    def H(n):
        if n not in entropy:
            if n in (0,DEN):entropy[n]=arb(0)
            else:
                w=ar(F(n,DEN));z=ar(F(DEN-n,DEN));counts['Arb_entropy_logs']+=2
                require(counts['Arb_entropy_logs']<=5632,'ENTROPY_LOG_CAP');entropy[n]=-w*w.log()-z*z.log()
        return entropy[n]
    write_new('started.json',{'round':'140','pid':os.getpid(),'nice':os.getpriority(os.PRIO_PROCESS,0),'release_sha256':sha(P/'parent_release.json'),'precision_bits':512,'time_limit':None,'source_cache_read_only':True})
    start=time.monotonic();results=[]
    try:
        for leaf in selection['selected_leaves']:
            vertices=points(leaf['box']);arrays=[[[getlog(k,p)for p in vertices]for k in keys]for keys in leaf['source_keys']]
            pairs=[];best=None
            for i,a in enumerate(arrays[0]):
                for j,b in enumerate(arrays[1]):
                    n,calls,float_score,float_entropy_logs=choose_weight([x[1]for x in a],[x[1]for x in b]);counts['float_objective_calls']+=calls;counts['float_entropy_logs']+=float_entropy_logs
                    w=ar(F(n,DEN));z=ar(F(DEN-n,DEN));ent=H(n)
                    scores=[endpoints(ent+w*x[0]+z*y[0])for x,y in zip(a,b)]
                    counts['Arb_corner_score_combinations']+=32;counts['pair_certificates']+=1
                    lower=min(v[0]for v in scores);upper=min(v[1]for v in scores)
                    row={'band0_source_index':i,'band1_source_index':j,'weight_numerator':n,'weight_denominator':DEN,'score_min_lower':str(lower),'score_min_upper':str(upper),'all_corner_score_lowers_nonnegative':lower>=0,'float_selector_score_diagnostic_only':float_score}
                    pairs.append(row)
                    if best is None or lower>F(best['score_min_lower']):best={**row,'band0_semantic_key':leaf['source_keys'][0][i],'band1_semantic_key':leaf['source_keys'][1][j],'corner_score_bounds':[[str(l),str(u)]for l,u in scores]}
            results.append({'root_id':leaf['root_id'],'leaf_id':leaf['leaf_id'],'box':leaf['box'],'whole_leaf_pass':F(best['score_min_lower'])>=0,'selected_best_certificate':best,'fixed_pair_results':pairs})
        require(counts['pair_certificates']==2816 and counts['float_objective_calls']==197120 and counts['float_entropy_logs']==382976,'EXACT_FINITE_BATCH_COUNTS')
        require(counts['Arb_corner_score_combinations']==90112 and counts['Arb_source_logs']+counts['Arb_entropy_logs']<=14176,'ARITHMETIC_BUDGET')
        root_completion=[{'root_id':ref['root_id'],'all_original_partition_leaves_paid':all(r['whole_leaf_pass']for r in results if r['root_id']==ref['root_id']),'old_paid_leaves_referenced_only':True}for ref in selection['original_closed_partition_references']]
        result={'status':'FINITE_COMMON_ENTROPY_BATCH_COMPLETE','target':'1/63500','selected_leaf_count':11,'passed_leaf_count':sum(r['whole_leaf_pass']for r in results),'leaf_results':results,'original_seven_local_box_completion':root_completion,'actual_counts':counts,'wall_seconds':time.monotonic()-start,'source_input_hashes_sha256':sha(P/'input_hashes.json'),'release_sha256':sha(P/'parent_release.json'),'global_g_certificate':False,'new_C6_bound':False,'scope':'only C134 eleven unresolved leaves; 62 paid leaves reference only; no family optimality or global claim'}
        write_new('result.json',result)
        write_new('exit_receipt.json',{'status':'FINITE_COMMON_ENTROPY_BATCH_COMPLETE','result_sha256':sha(P/'result.json'),'actual_counts':counts,'global_g_certificate':False,'new_C6_bound':False})
        print(json.dumps({'status':result['status'],'passed_leaf_count':result['passed_leaf_count'],'actual_counts':counts,'wall_seconds':result['wall_seconds']}))
    except BaseException:
        write_new('exit_receipt.json',{'status':'FAILED_STOP_NO_RETRY','failure':traceback.format_exc(),'actual_counts':counts,'global_g_certificate':False,'new_C6_bound':False});raise
if __name__=='__main__':
    if sys.argv[1:]==['--saved-only-preflight']:
        selection,cache=preflight();print(json.dumps({'status':'SAVED_ONLY_PREFLIGHT_PASS','leaves':len(selection['selected_leaves']),'source_cache_entries':len(cache),'new_log_calls':0,'new_exp_calls':0}))
    elif not sys.argv[1:]:run()
    else:raise SystemExit('usage: run_batch.py [--saved-only-preflight]')
