"""Prepared binary-root/secant/long-division/full-tree audit;no old consumer."""
from pathlib import Path
from datetime import datetime,timezone
import argparse,gzip,json,os,subprocess,sys,time
from roots175 import witness
from division175 import units160,ratio160
from io175 import sha,read,put,lines,lock
HERE=Path(__file__).resolve().parent

def setup(psha,csha):
    assert sha(HERE/'plan.json')==psha and sha(HERE/'candidate_binding.json')==csha
    p,c=read(HERE/'plan.json'),read(HERE/'candidate_binding.json')
    assert p['maximum_actual_audits']==1 and p['workers']==15
    assert p['candidate_binding_sha256']==csha and c['execution_ready'] is True
    assert sha(HERE/'source_binding.json')==p['source_binding_sha256']==c['source_binding_sha256']
    assert c['producer_plan_sha256']==p['producer_plan_sha256']
    alias=p['typed_producer_identity_alias']
    assert c['typed_input_preflight']['sha256']==alias['actual_input_preflight_sha256']
    assert c['typed_small_selfcheck']['sha256']==alias['small_selfcheck_sha256']
    assert c['producer_identity']['preflight_sha256']==alias['historical_field_sha256']==alias['small_selfcheck_sha256']
    for f,h in p['code_locks'].items():assert sha(HERE/f)==h,f
    return read(HERE/'source_binding.json'),c,{'plan_sha256':psha,'candidate_binding_sha256':csha,'code_locks':p['code_locks']}

def worker(out,i,psha,csha):
    b,c,ident=setup(psha,csha);idx=b['shards'][i];assert sha(idx['path'])==idx['sha256']
    dest=out/f'worker{i:02d}';dest.mkdir(exist_ok=False);start=time.monotonic();count=improved=selected=0;updates=[]
    pr=read(Path(c['candidate_dir'])/f'worker{i:02d}/result.json')
    assert pr['status']=='complete' and pr['worker']==i and pr['identity']==c['producer_identity']
    put(dest/'registration.json',{'pid':os.getpid(),'worker':i,'identity':ident,'source_shard':idx})
    put(dest/'ready.json',{'pid':os.getpid(),'before_any_new_root':True})
    while not (out/'START').exists():
        if (out/'STOP').exists():raise RuntimeError('peer failure before start')
        time.sleep(.02)
    def published():
        for chunk in pr['chunks']:
            assert sha(chunk['path'])==chunk['sha256'];rows=read(chunk['path']);assert len(rows)==chunk['count'];yield from rows
    actual=iter(published());f=gzip.open(dest/'root_slope_gain_checks.jsonl.gz','xt')
    try:
        for s in lines(idx['path']):
            if (out/'STOP').exists():raise RuntimeError('peer failure;close successful gzip prefix')
            assert s['strict112_active'] is True and s['self112']<s['certificate112_old_selfU']
            result,proof=witness(s);assert next(actual)==result
            f.write(json.dumps({**result,'independent_proof':proof},sort_keys=True,separators=(',',':'))+'\n')
            updates.append([s['key'],result['new_selfU']]);count+=1
            improved+=result['strict_self_gain_units']>0;selected+=result['selected_new_secant']
            if count%256==0:f.flush();put(dest/f'prefix{count:05d}.json',{'completed':count,'new_fourth_roots':2*count,'old_math_replays':0})
        assert next(actual,None) is None and count==idx['count']==pr['completed']
        assert pr['new_fourth_roots']==2*count and pr['new_gains']==count
        assert improved==pr['self_improved_vs166'] and selected==pr['selected_new_secant']
        f.close();assert sha(idx['path'])==idx['sha256'];put(dest/'updates.json',updates)
        put(dest/'result.json',{'status':'complete','pid':os.getpid(),'worker':i,'identity':ident,'completed':count,
            'new_fourth_roots':2*count,'new_slopes_gains':count,'improved':improved,'selected':selected,
            'seconds':time.monotonic()-start,'proof_sha256':sha(dest/'root_slope_gain_checks.jsonl.gz'),'updates_sha256':sha(dest/'updates.json')})
    except BaseException as e:
        f.close();put(dest/'fault.json',{'pid':os.getpid(),'completed':count,'exception':repr(e),'successful_prefix_preserved':True});raise

def aggregate(out,b,c,ident):
    assert sha(b['tree']['path'])==b['tree']['sha256'];tree=read(b['tree']['path']);nodes=tree['nodes'];updates={};workers=[]
    for i in range(15):
        wd=out/f'worker{i:02d}';r=read(wd/'result.json');assert r['status']=='complete' and r['identity']==ident
        assert sha(wd/'updates.json')==r['updates_sha256'] and sha(wd/'root_slope_gain_checks.jsonl.gz')==r['proof_sha256']
        values=read(wd/'updates.json');assert len(values)==r['completed']
        for k,v in values:assert k not in updates;updates[k]=v
        workers.append(r)
    assert len(updates)==97971 and sum(r['new_fourth_roots'] for r in workers)==195942
    assert set(nodes)-set(updates)=={'5:31,31,31,31'}
    roots=tree['roots'];assert len(roots)==len(set(roots))==52360
    assert set(roots)=={k for k,n in nodes.items() if n['parent'] is None}
    for k,n in nodes.items():
        assert len(n['children'])==len(set(n['children']))
        for ch in n['children']:assert nodes[ch]['parent']==k and nodes[ch]['level']==n['level']+1
        if n['parent'] is not None:assert k in nodes[n['parent']]['children']
    visiting=set();visited=set();eff={};selves={}
    def dfs(k):
        assert k not in visiting
        if k in visited:return eff[k]
        visiting.add(k);n=nodes[k];s=updates.get(k,n['self112']);assert n['selfL']<=s<=n['self166']
        value=min(s,sum(dfs(ch) for ch in n['children'])) if n['children'] else s
        assert 0<=value<=n['eff166'];eff[k]=value;selves[k]=s;visiting.remove(k);visited.add(k);return value
    units=sum(dfs(k) for k in roots);assert len(visited)==len(nodes)==97972
    raw=read(Path(c['candidate_dir'])/'effective_nodes.json');actual={r['key']:r for r in raw};assert len(actual)==len(raw)==97972
    expected={k:{'key':k,'new_selfU':selves[k],'new_effU':eff[k],'baseline166_selfU':n['self166'],'baseline166_effU':n['eff166'],
        'strict_self_gain_units':n['self166']-selves[k],'strict_eff_gain_units':n['eff166']-eff[k]} for k,n in nodes.items()}
    assert actual==expected
    r=read(Path(c['candidate_dir'])/'results.json');base=units160(tree['baseline166_UH'])
    assert units<=base and r['new_UH']==ratio160(units) and r['baseline166_UH']==tree['baseline166_UH']
    assert r['exact_gain_vs166']==ratio160(base-units)
    assert r['self_improved_vs166']==sum(w['improved'] for w in workers)
    assert r['new_secant_selected']==sum(w['selected'] for w in workers)
    assert r['effective_improved_vs166']==sum(eff[k]<n['eff166'] for k,n in nodes.items())
    assert r['roots_improved_vs166']==sum(eff[k]<nodes[k]['eff166'] for k in roots)
    assert r['new_fourth_roots']==195942 and r['new_gains']==97971 and r['active_nodes']==97971
    assert r['nodes']==97972 and r['roots']==52360 and r['unchanged_self_nodes']==1 and r['nonactive_key']=='5:31,31,31,31'
    assert r['gain_baseline']=='112 self minus one total floor relative112 mean;not166 self' and r['global_Cn_claim'] is False
    put(out/'tree_checks.json',{'nodes':97972,'roots':52360,'root_sum_units160':units,'effective_nodes':expected})
    return {'new_UH':ratio160(units),'exact_gain_vs166':ratio160(base-units),
        'self_improved_vs166':r['self_improved_vs166'],'effective_improved_vs166':r['effective_improved_vs166'],'roots_improved_vs166':r['roots_improved_vs166']}

def stop_failed_workers(out,procs):
    """Fault cleanup only. Normal mathematical execution has no deadline."""
    put(out/'STOP',{'reason':'peer/controller fault;preserve prefix'})
    def wait_phase(seconds):
        end=time.monotonic()+seconds
        while any(p.poll() is None for p in procs) and time.monotonic()<end:time.sleep(.01)
    errors=[];wait_phase(.25)
    for operation,seconds in [('terminate',2.0),('kill',2.0)]:
        for p in procs:
            if p.poll() is None:
                try:getattr(p,operation)()
                except OSError as e:errors.append({'pid':p.pid,'operation':operation,'exception':repr(e)})
        wait_phase(seconds)
    for p in procs:
        try:p.wait(timeout=.01)
        except subprocess.TimeoutExpired:pass
    return {'exitcodes':[p.returncode for p in procs],'still_alive_pids':[p.pid for p in procs if p.poll() is None],
        'signal_errors':errors,'fault_grace_terminate_kill_seconds':[.25,2.0,2.0],'final_join_timeout_each':.01}

def run(out,psha,csha):
    b,c,ident=setup(psha,csha);lock(b['source_locks']);lock(c['candidate_locks'])
    for row in b['shards']+[b['tree']]:assert sha(row['path'])==row['sha256']
    out.mkdir(exist_ok=False);start=time.monotonic();procs=[]
    put(out/'registration.json',{'task':'175_unique_independent_binaryroot_audit','pid':os.getpid(),'utc':datetime.now(timezone.utc).isoformat(),
        'identity':ident,'workers':15,'new_fourth_roots':195942,'new_slopes_gains':97971,'old_math_calls':0,'wall_cap':None})
    try:
        for i in range(15):
            log=(out/f'worker{i:02d}.log').open('xb');proc=subprocess.Popen([sys.executable,'-B',__file__,'worker','--output',str(out),
                '--plan-sha',psha,'--candidate-binding-sha',csha,'--worker',str(i)],stdout=log,stderr=subprocess.STDOUT);log.close();procs.append(proc)
        while not all((out/f'worker{i:02d}/ready.json').exists() for i in range(15)):
            assert all(p.poll() is None for p in procs),'worker failed before ready';time.sleep(.02)
        pids=[read(out/f'worker{i:02d}/ready.json')['pid'] for i in range(15)];assert pids==[p.pid for p in procs] and len(set(pids))==15
        put(out/'worker_pids.json',{'pids':pids,'all15_ready_before_first_root':True});put(out/'START',{'go':True})
        while any(p.poll() is None for p in procs):
            assert all(p.poll() in (None,0) for p in procs),'worker failure;preserve prefix';time.sleep(.05)
        assert [p.returncode for p in procs]==[0]*15
        answer=aggregate(out,b,c,ident);lock(b['source_locks']);lock(c['candidate_locks']);setup(psha,csha)
        for row in b['shards']+[b['tree']]:assert sha(row['path'])==row['sha256']
        put(out/'results.json',{'task':'175_independent','status':'all_new_roots_slopes_gains_and_complete_tree_qualified_pending_primary',
            'pid':os.getpid(),'seconds':time.monotonic()-start,**answer,'new_fourth_roots':195942,'new_slopes_gains':97971,'nodes':97972,'roots':52360,
            'worker_pids':pids,'worker_exitcodes':[p.returncode for p in procs],'all_source_hashes_before_after_equal':True,
            'old_H_log_mean_slope_witness_power_DFS_calls':0,'global_Cn_claim':False})
        print('PASS175 binary-root/secant/long-division/full-tree independent audit')
    except BaseException as e:
        cleanup=stop_failed_workers(out,procs)
        put(out/'fault.json',{'pid':os.getpid(),'exception':repr(e),'worker_exitcodes':cleanup['exitcodes'],'cleanup':cleanup,
            'successful_prefix_preserved':True,'no_automatic_retry':True});raise

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mode',choices=['run','worker']);p.add_argument('--output',type=Path,required=True)
    p.add_argument('--plan-sha',required=True);p.add_argument('--candidate-binding-sha',required=True);p.add_argument('--worker',type=int);a=p.parse_args()
    if a.mode=='run':run(a.output,a.plan_sha,a.candidate_binding_sha)
    else:worker(a.output,a.worker,a.plan_sha,a.candidate_binding_sha)
