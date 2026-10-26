#!/usr/bin/env python3
"""Prepared156 only:292 new factors,97907 means, fixed-tree min recursion.

Never run without primary authorization of the frozen plan. Fifteen workers
stream disjoint source shards. Only integer products are used for node means.
No H/R/log/gradient/exp/power-root function is imported or evaluated.
"""
from pathlib import Path
from fractions import Fraction as F
from datetime import datetime,timezone
import argparse,json,os,resource,subprocess,sys,time
from common import sha,read,put,lines,endpoint_arrays,lockcheck,ceildiv

HERE=Path(__file__).resolve().parent
N=1<<160;D=1<<192

def identity(plan_sha,plan):
    assert sha(HERE/'consumer.py')==plan['consumer_sha256'] and sha(HERE/'common.py')==plan['common_sha256']
    assert sha(HERE/'sources.json')==plan['sources_sha256']
    binding=plan['semantic_binding'];assert sha(binding['path'])==binding['sha256']
    b=read(binding['path']);assert b['all_shards_and_topology_semantically_bound'] and b['preflight_results_sha256']==plan['preflight_results']['sha256']
    assert b['sources_sha256']==plan['sources_sha256'] and b['new_factor_product_mean_gain_UH_evaluations']==0
    return {'plan_sha256':plan_sha,'consumer_sha256':plan['consumer_sha256'],'common_sha256':plan['common_sha256'],
        'sources_sha256':plan['sources_sha256'],'preflight_results_sha256':plan['preflight_results']['sha256'],'semantic_binding_sha256':binding['sha256']}

def new_factors(out,plan,pre,src):
    geos=read(src['112_geometries']['path']);old=read(src['112_factors']['path']);cat=read(src['112_ratios']['path'])
    times=pre['times'];factors={'g002':read(src['150_factors']['path'])};made=[]
    assert len(factors['g002'])==73
    for g in ['g000','g001','g003','g004']:
        prof=geos[g]['profile'];h,mu=F(prof['h']),F(prof['mu_infinity']);factors[g]=[]
        def endpoint(t):
            if t is None:return F(1),F(1)
            a=F(t)
            return F(cat[str(h/a)]['upper']),1+mu/a+F(25,24)*F(old[g][t]['J_upper'])/a**2
        for j,astr in enumerate(times):
            bstr=times[j+1] if j<72 else None;za=1/F(astr);zb=F(0) if bstr is None else 1/F(bstr)
            assert za>zb>=0
            w=(F(4,5)*(za**5-zb**5)/(za**4-zb**4)-zb)/(za-zb);assert 0<w<1
            ra,ta=endpoint(astr);rb,tb=endpoint(bstr);rb,tb=min(ra,rb),min(ta,tb)
            R=w*ra+(1-w)*rb;T=w*ta+(1-w)*tb;oldC=F(old[g][astr]['factor_upper'])
            assert (oldC*D).denominator==1
            exact=min(oldC,R,T);C=F(ceildiv((exact*D).numerator,(exact*D).denominator),D)
            assert 1<=C<=oldC
            item={'a':astr,'b':bstr,'theta':str(w),'range_time_upper':str(R),'taylor_time_upper':str(T),'old_upper':str(oldC),'new_upper':str(C)}
            factors[g].append(item);made.append({'geometry':g,**item})
            put(out/'factor_state.json',{'stage':'factors','completed_new_factors':len(made),'new_factors':made},True)
    assert len(made)==292 and sum(map(len,factors.values()))==365
    put(out/'factors.json',{'reused_g002_source_sha256':src['150_factors']['sha256'],'new_factor_count':292,'geometry_factors':factors})

def mean_worker(out,i,reviewed_sha):
    plan=read(HERE/'plan.json');assert sha(HERE/'plan.json')==reviewed_sha
    ident=identity(reviewed_sha,plan);reg=read(out/'registration.json');assert reg['identity']==ident
    pre=read(plan['preflight_results']['path']);assert sha(plan['preflight_results']['path'])==plan['preflight_results']['sha256']
    qreg=read(pre['external_files']['registration.json']['path']);idx=qreg['shards'][i];assert sha(idx['path'])==idx['sha256']
    fr=read(out/'factor_ready.json');assert sha(out/'factors.json')==fr['factors_sha256']
    catalog=read(out/'factors.json')['geometry_factors'];coeff={}
    for g,rows in catalog.items():
        coeff[g]=[];assert len(rows)==73
        for r in rows:
            x=F(r['new_upper'])*D;assert x.denominator==1;coeff[g].append(int(x))
    dest=out/f'worker{i:02d}';dest.mkdir(exist_ok=False);start=time.monotonic();completed=0;chunks=[];buffer=[];improved=0
    put(dest/'registration.json',{'identity':ident,'worker':i,'pid':os.getpid(),'source_index':idx,'factors_sha256':fr['factors_sha256'],'chunk_size':256,'no_automatic_retry':True})
    def checkpoint(partial=False):
        nonlocal buffer
        if not buffer:return
        p=dest/f"chunk{len(chunks):04d}.json";put(p,buffer)
        chunks.append({'path':str(p),'sha256':sha(p),'count':len(buffer),'partial_on_fault':partial});buffer=[]
        put(dest/'state.json',{'identity':ident,'completed':completed,'chunks':chunks,'mean_improved':improved,'status':'prefix'},True)
    try:
        for row in lines(idx['path']):
            p=Path(row['packet_path']);assert sha(p)==row['packet_sha256']
            lo,up,E=endpoint_arrays(read(p),pre['times'],row['point'],row['key'])
            cf=coeff[row['geometry']]
            total=sum((up[j]-(lo[j+1] if j<72 else 0))*cf[j] for j in range(73))
            newmi=min(row['old_mean_i'],ceildiv(total,1<<(E+32)))
            assert row['old_L_i']<=newmi<=row['old_mean_i']<=row['old_U_i']
            gain=(row['weight']*(row['old_mean_i']-newmi))>>(4*row['level']+4)
            newself=row['old_selfU']-gain;assert row['old_selfL']<=newself<=row['old_selfU']
            buffer.append({'key':row['key'],'old_mean_i':row['old_mean_i'],'new_mean_i':newmi,'old_selfU':row['old_selfU'],'new_selfU':newself,
                'gain_units':gain,'packet_sha256':row['packet_sha256']})
            completed+=1;improved+=newmi<row['old_mean_i']
            if len(buffer)==256:checkpoint()
        checkpoint();assert completed==idx['count']
        for row in lines(idx['path']):assert sha(row['packet_path'])==row['packet_sha256']
        assert sha(idx['path'])==idx['sha256'] and sha(out/'factors.json')==fr['factors_sha256']
        result={'identity':ident,'worker':i,'pid':os.getpid(),'status':'complete','completed':completed,'mean_improved':improved,'chunks':chunks,
            'cached_band_products':73*completed,'source_packets_before_after_unchanged':True,'seconds':time.monotonic()-start,
            'maxrss_native_units':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
        put(dest/'result.json',result);put(dest/'state.json',result,True)
    except BaseException as ex:
        partial=None
        if buffer:
            p=dest/'fault_prefix.json';put(p,buffer);partial={'path':str(p),'sha256':sha(p),'count':len(buffer)}
        put(dest/'fault.json',{'pid':os.getpid(),'exception':repr(ex),'completed':completed,'chunks':chunks,'additional_fault_prefix':partial,
            'action':'stop;prefix preserved;no automatic retry'})
        raise

def aggregate(out,plan,pre):
    qreg=read(pre['external_files']['registration.json']['path']);tp=qreg['topology'];assert sha(tp['path'])==tp['sha256'];tree=read(tp['path'])
    nodes={n['key']:n for n in tree['nodes']};updates={};worker_results=[];ident=read(out/'registration.json')['identity']
    for i in range(15):
        wr=read(out/f'worker{i:02d}/result.json');assert wr['status']=='complete' and wr['identity']==ident and wr['worker']==i
        worker_results.append(wr);count=0;expected=iter(lines(qreg['shards'][i]['path']))
        for chunk in wr['chunks']:
            assert sha(chunk['path'])==chunk['sha256'];data=read(chunk['path']);assert len(data)==chunk['count'];count+=len(data)
            for r in data:
                source=next(expected);assert r['key']==source['key'] and r['packet_sha256']==source['packet_sha256'] and r['old_mean_i']==source['old_mean_i']
                key=r['key'];assert key not in updates and r['old_selfU']==nodes[key]['selfU'];updates[key]=r['new_selfU']
        assert count==wr['completed']==qreg['shards'][i]['count']
    assert len(updates)==97907 and sum(x['cached_band_products'] for x in worker_results)==7147211
    reused=read(qreg['reuse150']['path']);assert sha(qreg['reuse150']['path'])==qreg['reuse150']['sha256']
    for r in reused:
        key=r['key'];assert key not in updates and nodes[key]['parent'] is None and not nodes[key]['children']
        assert nodes[key]['selfU']==r['old_selfU'];updates[key]=r['new_selfU']
    assert len(reused)==64 and len(updates)==97971
    eff={};changes=[]
    for key in sorted(nodes,key=lambda k:(-nodes[k]['level'],k)):
        n=nodes[key];selfU=updates.get(key,n['selfU']);children=n['children']
        if children:
            assert all(nodes[ch]['parent']==key and nodes[ch]['level']==n['level']+1 for ch in children)
            value=min(selfU,sum(eff[ch] for ch in children))
        else:value=selfU
        assert n['selfL']<=selfU<=n['selfU'] and value<=n['effU']
        eff[key]=value
        if value<n['effU']:changes.append({'key':key,'old_effU':n['effU'],'new_effU':value})
    units=sum(eff[k] for k in tree['roots']);old=F(tree['old_complete_UH'])*N;prev=F(read(read(HERE/'sources.json')['source_locks']['150_results']['path'])['new_UH'])*N
    assert old.denominator==prev.denominator==1 and units<=prev<=old
    put(out/'effective_changes.json',changes)
    result={'task':'upper_integral156','status':'finite_fixedtree_candidate_pending_independent_primary','pid':os.getpid(),'workers':15,'worker_pids':[w['pid'] for w in worker_results],
        'new_time_factors':292,'reused_time_factors':73,'new_mean_nodes':97907,'reused150_roots':64,'unchanged_self_nodes':1,'new_band_products':7147211,
        'mean_improved_nodes':sum(w['mean_improved'] for w in worker_results),'effective_improved_nodes':len(changes),'original_tree_nodes':97972,'roots':52360,
        'old112_UH':str(F(int(old),N)),'accepted150_UH':str(F(int(prev),N)),'new_UH':str(F(units,N)),
        'exact_gain_vs112':str(F(int(old)-units,N)),'exact_incremental_gain_vs150':str(F(int(prev)-units,N)),
        'new_H_R_gradient_log_exp_powerroot_calls':0,'old_baseline_recomputations':0,'global_Cn_claim':False}
    put(out/'results.json',result)

def run(out,reviewed_sha):
    assert sha(HERE/'plan.json')==reviewed_sha;plan=read(HERE/'plan.json');ident=identity(reviewed_sha,plan)
    pre=read(plan['preflight_results']['path']);assert sha(plan['preflight_results']['path'])==plan['preflight_results']['sha256']
    src=read(HERE/'sources.json')['source_locks'];lockcheck(src)
    for x in pre['external_files'].values():assert sha(x['path'])==x['sha256']
    assert pre['counts']=={'new_mean_nodes':97907,'reused150':64,'unchanged':1} and pre['new_factors_products_means_gains_UH']==0
    out.mkdir(exist_ok=False);start=time.monotonic()
    put(out/'registration.json',{'task':'upper_integral156_unique_fixedtree_batch','identity':ident,'pid':os.getpid(),'utc':datetime.now(timezone.utc).isoformat(),
        'workers':15,'new_time_factors':292,'reused150_factors':73,'new_mean_nodes':97907,'cached_band_products':7147211,'reused150_roots':64,'wall_cap':None,
        'chunk_size':256,'new_H_R_gradient_log_exp_powerroot_calls':0,'old_baseline_recomputations':0,'global_Cn_claim':False})
    procs=[]
    try:
        new_factors(out,plan,pre,src);put(out/'factor_ready.json',{'factors_sha256':sha(out/'factors.json')})
        for i in range(15):
            log=(out/f'worker{i:02d}.log').open('xb');p=subprocess.Popen([sys.executable,'-B',__file__,'worker','--output',str(out),'--reviewed-plan-sha',reviewed_sha,'--worker',str(i)],stdout=log,stderr=subprocess.STDOUT);log.close();procs.append(p)
        put(out/'worker_pids.json',{'pids':[p.pid for p in procs]})
        while any(p.poll() is None for p in procs):
            if any(p.poll() not in (None,0) for p in procs):raise RuntimeError('Worker failure:stop fixed batch')
            time.sleep(.2)
        assert [p.returncode for p in procs]==[0]*15
        aggregate(out,plan,pre);lockcheck(src);assert identity(reviewed_sha,plan)==ident
        assert sha(HERE/'plan.json')==reviewed_sha and sha(plan['preflight_results']['path'])==plan['preflight_results']['sha256']
        for x in pre['external_files'].values():assert sha(x['path'])==x['sha256']
        result=read(out/'results.json');result.update({'source_locks_before_after_unchanged':True,'seconds':time.monotonic()-start,
            'registration_sha256':sha(out/'registration.json'),'factors_sha256':sha(out/'factors.json')})
        # Completion is an additional file; never overwrite candidate results.
        put(out/'completion.json',result);print(json.dumps(result,sort_keys=True))
    except BaseException as ex:
        for p in procs:
            if p.poll() is None:p.terminate()
        for p in procs:p.wait()
        put(out/'fault.json',{'pid':os.getpid(),'exception':repr(ex),'completed_prefix_chunks_preserved':True,'action':'stop;no automatic retry,refinement,new H or scope expansion'})
        raise

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('mode',choices=['run','worker']);ap.add_argument('--output',type=Path,required=True)
    ap.add_argument('--reviewed-plan-sha',required=True);ap.add_argument('--worker',type=int);a=ap.parse_args()
    if a.mode=='run':run(a.output,a.reviewed_plan_sha)
    else:mean_worker(a.output,a.worker,a.reviewed_plan_sha)
