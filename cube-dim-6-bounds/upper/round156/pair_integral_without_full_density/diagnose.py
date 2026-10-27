#!/usr/bin/env python3
"""One gated whole-joint-law bound family; no density compiler or6D grid."""
from pathlib import Path
from fractions import Fraction as Q
from itertools import combinations,product
from math import prod
import gzip,hashlib,json,os,sys,time,traceback
sys.set_int_max_str_digits(0)
P=Path(__file__).resolve().parent
def read(p):return json.loads(Path(p).read_text())
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def req(ok,msg):
    if not ok:raise ArithmeticError(msg)
def save(n,x):
    with(P/n).open('x')as f:json.dump(x,f,indent=2);f.write('\n');f.flush();os.fsync(f.fileno())
def multiply(intervals):
    rows=list(intervals)
    return prod(z[0]for z in rows),prod(z[1]for z in rows)
def min_interval(rows):return min(z[0]for z in rows),min(z[1]for z in rows)
def main():
    inp=read(P/'inputs.json');reg=read(P/'registration.json');rel=read(P/'parent_release.json')
    req(rel.get('round')=='156'and rel.get('authorized')is True and rel.get('main_approved')is True and rel.get('independent_review_approved')is True,'DUAL_GATE')
    for field,name in [('diagnose_sha256','diagnose.py'),('proof_sha256','PROOF.md'),('inputs_sha256','inputs.json'),('source_manifest_sha256','source_manifest.json')]:req(reg[field]==sha(P/name)==rel.get(field),'FROZEN_'+name)
    req(rel.get('registration_sha256')==sha(P/'registration.json'),'REGISTERED_HASH')
    source_rows=read(P/'source_manifest.json')['sources']
    for source in source_rows:req(sha(source['path'])==source['sha256'],'IMMUTABLE_SOURCE')
    source_paths={Path(row['path']).name:Path(row['path'])for row in source_rows}
    req(inp['cones']==read(source_paths['cone_table.json'])['cones'],'TRUE_FROZEN_CONE_TABLE')
    req(inp['p_interval']==read(source_paths['upstream.json'])['p_intervals'][5]and inp['kappa_interval']==read(source_paths['results.json'])['pi_factor_interval'],'EXACT_SAVED_SCALAR_INTERVALS')
    req(not(P/'started.json').exists(),'NO_RESTART');os.nice(max(0,10-os.getpriority(os.PRIO_PROCESS,0)))
    req(all(os.environ.get(k)=='1'for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS']),'ONE_THREAD')
    save('started.json',{'round':'156','pid':os.getpid(),'nice':os.getpriority(os.PRIO_PROCESS,0),'workers':1,'time_or_CPU_limit':None,'release_sha256':sha(P/'parent_release.json')});start=time.monotonic();cpu=time.process_time();sqrt_cache={};counts={'raw_pairs':0,'caps_evaluated':0,'exact_column_determinants':0,'four_seed_pairing_products':0,'cauchy_sqrt_cache_misses':0,'cauchy_sqrt_pruned_by_exact_range':0,'new_I0_closed':0,'new_positive_k_anchor':0}
    try:
        import flint
        from flint import arb,ctx
        req(flint.__version__=='0.9.0'and flint.__FLINT_VERSION__=='3.6.0','PINNED_FLINT');ctx.prec=512
        def ar(q):q=Q(q);return arb(q.numerator)/q.denominator
        def ends(v):req(v.is_finite(),'FINITE_SCALAR');return Q(str(v.lower().fmpq())),Q(str(v.upper().fmpq()))
        h=ar(inp['h']);pL,pU=map(Q,inp['p_interval']);req(0<pL<=pU<1,'EXACT_SAVED_P');p=ar(pL).union(ar(pU));s=h.sqrt();pi=arb.pi()
        I0=h*(1+p)-p/2+s*(-h).exp()/pi.sqrt();iL,iU=ends(I0);counts['new_I0_closed']=1;req(0<iL<=iU,'POSITIVE_I0')
        e2=(2*s).exp();plus=(ar(Q(1,2))+s).erfc();minus=(ar(Q(1,2))-s).erfc();qstar=e2*plus;kstar=(minus-qstar)/2;req(kstar>0,'POSITIVE_TRUE_ANCHOR_K');tau1=qstar/kstar;C1=kstar*(tau1*s).exp();counts['new_positive_k_anchor']=1
        taux=[ends((1-p)/p),ends(tau1)];constants=[(pL,pU),ends(C1)];req(all(z[0]>0 for z in taux+constants),'POSITIVE_TANGENT_PARAMETERS')
        pp=(pL**2,pU**2);ii=(iL**2,iU**2)
        def tangent(b,c,r,t,normalized=False):
            lo=constants[r][0]*constants[t][0]/(1+taux[r][1]*b+taux[t][1]*c)
            hi=constants[r][1]*constants[t][1]/(1+taux[r][0]*b+taux[t][0]*c)
            return (lo/pp[1],hi/pp[0])if normalized else(lo,hi)
        def seed2(D):return min(pp[0],ii[0]/D),min(pp[1],ii[1]/D)
        def seed4(P4):
            if iL**4>=16*P4*pU**4:counts['cauchy_sqrt_pruned_by_exact_range']+=1;return pp
            if P4 not in sqrt_cache:
                counts['cauchy_sqrt_cache_misses']+=1;req(counts['cauchy_sqrt_cache_misses']<=156240,'FINITE_CAUCHY_ROOT_BUDGET');sqrt_cache[P4]=ends(ar(P4).sqrt())
            lo,hi=sqrt_cache[P4];req(lo>0 and lo**2<=P4<=hi**2,'EXACT_POSITIVE_CAUCHY_ROOT_CERTIFICATE');return min(pp[0],ii[0]/(4*hi)),min(pp[1],ii[1]/(4*lo))
        sums={'flat':(Q(0),Q(0)),'finite_best':(Q(0),Q(0))};headers=[];alphasum=Q(0);(P/'pairs').mkdir(exist_ok=False)
        for cone in inp['cones']:
            alpha=Q(cone['alpha']);req(alpha>0,'POSITIVE_ALPHA');alphasum+=alpha;nodes=[[Q(v)for v in row]for row in cone['nodes_by_body']];req(len(nodes)==32 and all(len(v)==6 and min(v)>=0 for v in nodes),'TRUE_NONNEGATIVE_2X6_COLUMNS');records=[];cone_sums={'flat':(Q(0),Q(0)),'finite_best':(Q(0),Q(0))}
            for i,j in combinations(range(32),2):
                b,c=nodes[i],nodes[j];D={}
                for a,d in combinations(range(6),2):D[(a,d)]=abs(b[a]*c[d]-b[d]*c[a]);counts['exact_column_determinants']+=1
                maxD=max(D.values());flat=seed2(maxD)if maxD>0 else pp
                singles=[]
                for a in range(6):
                    row=(Q(1),Q(1))if b[a]==c[a]==0 else min_interval([(Q(1),Q(1))]+[tangent(b[a],c[a],r,t,True)for r,t in product(range(2),repeat=2)])
                    singles.append(row)
                allS=multiply(singles);caps=[(pp[0]*allS[0],pp[1]*allS[1])];names=['all_single_columns']
                for key,det in D.items():
                    if det==0:continue
                    seed=seed2(det);rest=multiply(singles[a]for a in range(6)if a not in key);caps.append((seed[0]*rest[0],seed[1]*rest[1]));names.append('seed2:'+','.join(map(str,key)))
                for four in combinations(range(6),4):
                    a,d,e,f=four;pairings=[((a,d),(e,f)),((a,e),(d,f)),((a,f),(d,e))];products=[D[x]*D[y]for x,y in pairings];counts['four_seed_pairing_products']+=3;bestP=max(products)
                    if bestP==0:continue
                    seed=seed4(bestP);rest=multiply(singles[a]for a in range(6)if a not in four);caps.append((seed[0]*rest[0],seed[1]*rest[1]));names.append('seed4:'+','.join(map(str,four)))
                for r,t in product(range(2),repeat=2):
                    lo=constants[r][0]*constants[t][0];hi=constants[r][1]*constants[t][1]
                    for a in range(6):lo/=1+taux[r][1]*b[a]+taux[t][1]*c[a];hi/=1+taux[r][0]*b[a]+taux[t][0]*c[a]
                    caps.append((lo,hi));names.append('global_tangent:'+str(r)+str(t))
                req(len(caps)<=35 and all(0<z[0]<=z[1]for z in caps),'FIXED_POSITIVE_CAP_FAMILY');best=min_interval(caps);req(best[1]<=flat[1],'SEED_REMAINDERS_DO_NOT_WORSEN_FLAT');counts['raw_pairs']+=1;counts['caps_evaluated']+=len(caps)
                for name,val in [('flat',flat),('finite_best',best)]:old=cone_sums[name];cone_sums[name]=(old[0]+val[0],old[1]+val[1])
                records.append({'body_pair_zero_based':[i,j],'best_seed_determinant':str(maxD),'flat_budget_interval':list(map(str,flat)),'finite_best_budget_interval':list(map(str,best)),'caps_evaluated':len(caps),'selected_lower_cap':names[min(range(len(caps)),key=lambda n:caps[n][0])],'selected_upper_cap':names[min(range(len(caps)),key=lambda n:caps[n][1])]})
            for name,val in cone_sums.items():old=sums[name];sums[name]=(old[0]+alpha*val[0],old[1]+alpha*val[1])
            filename='pairs/cone_'+str(cone['id']).zfill(2)+'.json.gz';path=P/filename
            with path.open('xb')as raw:
                with gzip.GzipFile(filename='',mode='wb',fileobj=raw,mtime=0)as z:z.write(json.dumps(records,separators=(',',':')).encode())
            headers.append({'cone_id':cone['id'],'alpha':str(alpha),'raw_pairs':len(records),'file':filename,'sha256':sha(path),'unweighted_budget_sums':{k:list(map(str,v))for k,v in cone_sums.items()}})
        req(counts['raw_pairs']==10416 and 46080*alphasum==Q(inp['R0']),'COMPLETE_GROUP_REDUCED_FAN_AND_PAIR_ACCOUNT')
        kL,kU=map(Q,inp['kappa_interval']);goals={'goal7':Q(inp['goal7']),'goal8':Q(inp['goal8']),'saved146_legal_cap':Q(inp['saved146_legal_beta2_upper'])};budgets={}
        for name,(lo,hi)in sums.items():
            BL,BU=kL*lo,kU*hi;req(BU-BL<=Q(inp['global_budget_width_gate']),'OUTWARD_GLOBAL_WIDTH_GATE');budgets[name]={'budget_expression_interval':[str(BL),str(BU)],'valid_true_beta2_upper':str(BU),'budget_lower_is_actual_mu2_lower':False,'upper_comparisons':{k:BU<v for k,v in goals.items()},'this_fixed_budget_family_cannot_reach':{k:BL>=v for k,v in goals.items()}}
        cache_headers=[];cache_items=list(sqrt_cache.items())
        for offset in range(0,len(cache_items),8192):
            name='pairs/cauchy_roots_'+str(offset//8192).zfill(2)+'.json.gz';path=P/name;data={str(q):list(map(str,v))for q,v in cache_items[offset:offset+8192]}
            with path.open('xb')as raw:
                with gzip.GzipFile(filename='',mode='wb',fileobj=raw,mtime=0)as z:z.write(json.dumps(data,separators=(',',':')).encode())
            cache_headers.append({'file':name,'sha256':sha(path),'entries':len(data)})
        result={'status':'COMPLETE_ONE_TRUE_JOINT_LAW_BOUND_FAMILY_DIAGNOSTIC','h':'1/128','budgets':budgets,'p_interval':inp['p_interval'],'I0_interval':list(map(str,(iL,iU))),'tangent_tau_intervals':[list(map(str,v))for v in taux],'tangent_C_intervals':[list(map(str,v))for v in constants],'anchor_k_interval':list(map(str,ends(kstar))),'counts':counts,'scalar_calls':{'new_I0_closed':1,'new_positive_k_anchor':1,'exp':3,'erfc':2,'erf':0,'sqrt_h':1,'pi_constant':1,'sqrt_pi':1},'pair_files':headers,'source_alpha_sum':str(alphasum),'kappa_interval':inp['kappa_interval'],'all_support_integrated_no_tail_or_probability_renormalization':True,'density_generations':0,'old343456_functional_replays':0,'old6D_grid_calls':0,'MC_native_optimizer_calls':0,'actual_mu2_lower_bound':False,'new_M_upper':False,'global_g_proved':False,'new_C6_bound':False,'diagnose_sha256':sha(Path(__file__)),'release_sha256':sha(P/'parent_release.json'),'wall_before_receipt_seconds':time.monotonic()-start,'CPU_before_receipt_seconds':time.process_time()-cpu}
        result['cauchy_root_certificate_files']=cache_headers
        save('result.json',result);save('exit_receipt.json',{'status':result['status'],'result_sha256':sha(P/'result.json'),'counts':counts,'density_generation_old6D_MC_native_calls':0,'new_C6_bound':False});print(json.dumps({'status':result['status'],'budget_comparisons':{k:v['upper_comparisons']for k,v in budgets.items()},'raw_pairs':10416,'new_C6_bound':False}))
    except BaseException:
        save('exit_receipt.json',{'status':'FAILED_STOP_NO_RETRY','failure':traceback.format_exc(),'counts':counts,'new_C6_bound':False});raise
if __name__=='__main__':main()
