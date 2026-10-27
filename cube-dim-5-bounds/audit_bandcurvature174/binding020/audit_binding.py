#!/usr/bin/env python3
"""455 identity/qualification checks only. No173/174 arithmetic imports.

Fraction is used only to compare existing endpoints/points,not to form a
band mass,weight,K,q,log-gradient,corner,power or global margin.
"""
from pathlib import Path
from fractions import Fraction as F
from collections import Counter
import argparse,gzip,hashlib,json,multiprocessing as mp,os,queue,resource,time,traceback
H=Path(__file__).resolve().parent;T=Path('/Users/zhengzhihao/Desktop/T');B=T/'Math-Research-cube-dim-5/cube-dim-5-bounds';U=T/'Math-Research-cube5-upper/cube-dim-5-bounds'
def sha(p):
    h=hashlib.sha256()
    with Path(p).open('rb') as f:
        for b in iter(lambda:f.read(1<<20),b''):h.update(b)
    return h.hexdigest()
def digest(x):return hashlib.sha256(json.dumps(x,sort_keys=True,separators=(',',':')).encode()).hexdigest()
def read(p):
    p=Path(p)
    with (gzip.open(p,'rt') if p.suffix=='.gz' else p.open()) as f:return json.load(f)
def lines(p):
    p=Path(p)
    with (gzip.open(p,'rt') if p.suffix=='.gz' else p.open()) as f:
        for x in f:yield json.loads(x)
def put(p,x):
    with Path(p).open('x') as f:json.dump(x,f,sort_keys=True,separators=(',',':'));f.write('\n');f.flush();os.fsync(f.fileno())
def point(a,b):return len(a)==len(b)==5 and all(F(x)==F(y) for x,y in zip(a,b))
def one(cid,s,jet,co,event,common):
    old=s['old155'];ref=s['original_g151_jet_ref'];bd=s['donor_binding'];md=s['source171_metadata']['new_mass_donor'];model=bd['model_binding'];stage=model['stage'];reg=common['regs'][stage]
    assert co['canonical_id']==cid==digest(co['core']) and jet['jet_id']==co['jet_id']==ref['jet_id']
    payload={k:v for k,v in jet.items() if k!='minimal_projection_sha256'};assert digest(payload)==jet['minimal_projection_sha256']
    assert jet['epsilon']==ref['epsilon']==co['core']['epsilon']==old['epsilon']=='1/64' and jet['model']==ref['model']==common['model151']
    assert point(jet['point'],old['point']) and point(co['core']['point'],old['point'])
    assert jet['raw_record_sha256']==co['raw_record_sha256']==ref['raw_record_sha256']==old['raw_record_sha256']
    assert len(jet['values'])==5 and all(len(v)==2 and F(v[0])<=F(v[1]) for v in jet['values'])
    assert old['loggradient_directed']==co['core']['loggradient_directed'] and old['Delta110_lower']==co['core']['Delta110_lower']
    body=event['body'];assert digest(body)==old['source155_body_sha256'] and event['hash']==old['source155_event_hash']
    for n in ('epsilon','point','Delta110_lower','Delta115_lower','loggradient_directed','raw_record_sha256'):assert old[n]==body[n]
    assert old['oldqualified155q']==body['q_new'] and old['source155_q_rule']==body['q_new_rule']
    assert reg['model']==model['model'] and reg['table_sha256']==model['bodytable_sha256']==common['table_sha'] and reg['time_endpoints']==common['times'] and reg['epsilon']=='1/512'
    assert sha(md['packet_path'])==md['packet_sha256'];raw=read(md['packet_path']);assert raw['key']==bd['packet_key']
    assert point(raw['point'],old['point']) and point(raw['point'],md['point'])
    ends=raw['endpoints'];assert len(ends)==73 and [e['t'] for e in ends]==common['times']
    extra={}
    if bd['kind']=='legacy092_compact':
        nr=common['range092'][md['packet_path']];assert nr['key']==raw['key'] and nr['packet_sha256']==md['packet_sha256']
        assert nr['original_raw_path']==bd['original_raw_path'] and nr['original_raw_sha256']==bd['original_raw_sha256']==raw['original_raw_sha256']
        assert sha(bd['original_raw_path'])==bd['original_raw_sha256'];original=read(bd['original_raw_path'])
        assert point(original.get('point',original.get('anchor')),raw['point']) and original['endpoints']==ends
        assert stage=='upper_global080';rank=original['rank'];cell=common['cells080'][str(rank)];root=common['roots080'][str(rank)]
        assert cell['evidence_sha256']==bd['original_raw_sha256'] and original['indices']==root['indices']==nr['indices'] and point(root['point'],raw['point'])
        no_runtime={k:v for k,v in original.items() if k!='runtime'};assert digest(no_runtime)==cell['mathematical_sha256']
        policy=common['policy080'];expected=policy['old_fingerprint'] if rank<policy['transition_rank'] else common['newfp080'];assert original['fingerprint']==expected
        extra={'original_raw_sha256':bd['original_raw_sha256'],'old080_rank':rank,'old080_mathematical_payload_sha256':cell['mathematical_sha256']}
    else:
        assert bd['kind']=='full_adaptive_evidence' and stage in ('upper_adaptive095','upper_adaptive103')
        expected=hashlib.sha256((reg['run_sha256']+reg['kernel_sha256']+sha(model['registration_path'])).encode()).hexdigest()
        assert expected==model['expected_packet_fingerprint']==raw['fingerprint']
        j=common['jobs'][stage][raw['key']];cert=common['certs'][stage][raw['key']]
        assert cert['evidence_sha256']==md['packet_sha256'] and cert['bounds']==raw['bounds'] and point(j['point'],raw['point'])
        for n in ('level','indices','parent','weight'):assert raw[n]==j[n]
        for e in ends:
            if e.get('origin')=='exact_registered_cache':assert j['reuse_H'][e['t']]==e['H']
        extra={'production_fingerprint':expected,'job_key':raw['key'],'original_certificate_sha256':digest(cert)}
    for e in ends:assert 0<F(e['H']['lower'])<=F(e['H']['upper'])
    gl,gu=map(F,jet['values'][0]);dl,du=(F(ends[24]['H'][n]) for n in ('lower','upper'));assert 0<gl<=gu and max(gl,dl)<=min(gu,du)
    assert 0<=F(ends[32]['H']['lower'])<gu
    # Ordering only;do not construct even a band mass or any weight/product.
    assert all(F(ends[j]['H']['upper'])>F(ends[j+1]['H']['lower']) for j in range(24,72))
    for p,h in [(md['packet_path'],md['packet_sha256']),*([(bd['original_raw_path'],bd['original_raw_sha256'])] if 'original_raw_path' in bd else [])]:assert sha(p)==h
    return {'canonical_id':cid,'status':'PASS_identity_and_existing_qualification_binding','stage':stage,'donor_packet_sha256':md['packet_sha256'],'endpoint_array_sha256':digest(ends),
      'original_jet_projection_sha256':jet['minimal_projection_sha256'],'original_jet_numeric_source_id':jet['numeric_source_id'],'original_raw_jet_sha256':jet['raw_record_sha256'],
      'all73_positive_intervals':True,'H64_overlap_nonempty':True,'H32_lower_vs_original_g64_upper_consistent':True,'positive_band_upper_ordering_only':True,'original155_Delta_gamma_q_copied':True,
      'physical_H_derivative_component_order':['H','d_tail1','d_tail2','d_tail3','d_tail4'],'gradient_numeric_recomputed':False,'H_numeric_truth_shared':True,**extra}
def worker(sid,tasks,common,out,start,stop,ready,result):
    dest=Path(out)/f'worker{sid:02d}';dest.mkdir();t=time.monotonic();done=[]
    try:
        ready.put({'sid':sid,'pid':os.getpid(),'jobs':len(tasks)});start.wait()
        with (dest/'source_checks.jsonl').open('xb') as f:
            for args in tasks:
                if stop.is_set():break
                r=one(*args,common);f.write((json.dumps(r,sort_keys=True,separators=(',',':'))+'\n').encode());f.flush();os.fsync(f.fileno());done.append(r)
        u=resource.getrusage(resource.RUSAGE_SELF);r={'sid':sid,'pid':os.getpid(),'status':'complete' if len(done)==len(tasks) else 'stopped_successful_prefix','count':len(done),'seconds':time.monotonic()-t,'cpu_user':u.ru_utime,'cpu_system':u.ru_stime,'peak_rss_bytes':u.ru_maxrss,'source_checks_sha256':sha(dest/'source_checks.jsonl')};put(dest/'result.json',r);result.put(r)
    except BaseException as ex:
        stop.set();r={'sid':sid,'pid':os.getpid(),'status':'fault','completed':len(done),'exception':repr(ex),'traceback':traceback.format_exc()};put(dest/'fault.json',r);result.put(r);raise
def baseline(plan,registry):
    q=read(B/'audit_reuse168/run001/result.json');assert q['status']=='all_actual158_corners_firstpass_best_and_complete_domain_independently_qualified' and q['complete_old111201_domain_preserved'] and q['pending115']==1120
    bind=read(B/'audit_reuse168/candidate_binding001.json');prod=read(plan['base158_candidate_path']);assert prod['full115_passed']==110081 and prod['target115_pending']==1120
    for p in (plan['base158_candidate_path'],plan['base158_updates_path']):
        h=bind['candidate_locks'][p];assert sha(p)==(h['sha256'] if isinstance(h,dict) else h)
    base={r['leaf_index']:r for r in read(plan['base158_updates_path'])};leafs={};needprior158={};prior155=[];count=Counter();seen=set();receivers=set();used=set()
    for sp in plan['phase2_corner_shards']:
        d=read(sp['input_path']);assert d['shard_id']==sp['shard_id']
        for k,r in d['receivers'].items():
            i=int(k);assert i not in receivers and i%15==sp['shard_id'];receivers.add(i);assert r['base_best158']==base[r['leaf_index']]['best'] and r['node']==base[r['leaf_index']]['node']
        for pair in d['pairs']:
            i=pair['pair_index'];cid=pair['canonical_id'];assert i not in seen and cid in registry;seen.add(i);used.add(cid)
            assert pair['node']==d['receivers'][str(pair['receiver_index'])]['node'] and d['sources'][cid]==registry[cid]['source171_metadata']
            pr=pair['prior_test_same_old155_q']
            if pr:
                assert pr['new_q_strict_decrease_required'];count[pr['stage']]+=1
                if pr['stage']=='158':needprior158.setdefault(pr['journal_path'],[]).append((pair,pr))
                else:assert pr['stage']=='155';prior155.append((pair,pr));leafs[pr['source155_leaf_index']]=None
    assert seen==set(range(35872)) and receivers==set(range(1120)) and used==set(registry) and dict(count)=={'155':13,'158':58}
    for row in lines(plan['source155_margin_part']):
        if row['leaf_index'] in leafs:leafs[row['leaf_index']]=row
    for pair,pr in prior155:
        r=leafs[pr['source155_leaf_index']];assert r['leaf_index']==pair['leaf_index'] and r['node']==pair['node'] and r['canonical_id']==pair['canonical_id']
        assert r['margin110_new']==pr['margin110'] and r['margin115']==pr['margin115']
    for path,items in needprior158.items():
        assert sha(path)==items[0][1]['journal_sha256'];events={r['seq']:r for r in lines(path)}
        for pair,pr in items:
            e=events[pr['event_seq']];assert e['hash']==pr['event_hash'] and e['canonical_id']==pair['canonical_id'] and e['kind']=='evaluated'
            c=e['certificate'];assert c['node']==pair['node'] and c['margin110']==pr['margin110'] and c['margin115']==pr['margin115']
    return {'status':'PASS_existing168_baseline_metadata_linkage','receivers':1120,'fixed_pairs':35872,'guarded_old_pairs':71,'guard_stages':dict(count),'new_corner_or_old_margin_recalculations':0}
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--plan-sha',required=True);ap.add_argument('--output',type=Path,required=True);a=ap.parse_args();assert sha(H/'plan.json')==a.plan_sha;p=read(H/'plan.json')
    for n,h in p['code_locks'].items():assert sha(H/n)==h
    for path,h in p['source_locks'].items():assert sha(path)==h,path
    a.output.mkdir(exist_ok=False);began=time.monotonic();registry=read(p['registry_path']);finite=read(p['finite_plan_path']);assert len(registry)==455
    ids=set(registry);coeff={r['canonical_id']:r for r in lines(U/'upper_global151/data/coefficients-0001.jsonl.gz') if r['canonical_id'] in ids};wanted={s['original_g151_jet_ref']['jet_id']:cid for cid,s in registry.items()}
    jets={wanted[r['jet_id']]:r for r in lines(U/'upper_global151/data/jets-0001.jsonl.gz') if r['jet_id'] in wanted};events={}
    for path in sorted({s['old155']['source155_original_path'] for s in registry.values()}):
        for r in lines(path):
            cid=r['body']['canonical_id']
            if cid in ids:rr={k:v for k,v in r.items() if k!='hash'};assert digest(rr)==r['hash'];events[cid]=r
    assert set(coeff)==set(jets)==set(events)==ids
    regs={st:read(B/st/'registration.json') for st in ('upper_global080','upper_adaptive095','upper_adaptive103')};policy=read(B/'upper_global080/execution_policy_v2.json')
    assert sha(B/'upper_global080/engine.py')==regs['upper_global080']['engine_sha256']==policy['old_engine_sha256']
    assert sha(B/'upper_global080/engine_v2.py')==policy['new_engine_sha256'] and sha(B/'upper_global080/roots.json.gz')==regs['upper_global080']['roots_sha256']
    for st in ('upper_adaptive095','upper_adaptive103'):
        assert sha(B/st/'kernel.py')==regs[st]['kernel_sha256'] and sha(B/st/'run.py')==regs[st]['run_sha256']
    cells=read(B/'upper_global080/initial_cells.json.gz');ranks={int(Path(s['donor_binding']['original_raw_path']).stem.split('.')[0]) for s in registry.values() if s['donor_binding']['kind']=='legacy092_compact'}
    needed_packets={s['source171_metadata']['new_mass_donor']['packet_path'] for s in registry.values()};meta151=read(U/'upper_global151/metadata.json')
    assert meta151['bodytable_document_id']==sha(B/'upper_global080/geometry017.json.gz')=='8dd2e9895d7465ac599d43feb1b5099506f4e43d1248215eba0a77e11f4d13f9'
    common={'times':finite['time_endpoints'],'model151':meta151['model'],'table_sha':sha(B/'upper_global080/geometry017.json.gz'),'regs':regs,'policy080':policy,
      'newfp080':hashlib.sha256((policy['new_engine_sha256']+sha(B/'upper_global080/registration.json')+sha(B/'upper_global080/execution_policy_v2.json')).encode()).hexdigest(),
      'cells080':{str(row[0]):dict(zip(cells['columns'],row)) for row in cells['rows'] if row[0] in ranks},'roots080':{str(r['rank']):r for r in read(B/'upper_global080/roots.json.gz') if r['rank'] in ranks},
      'range092':{r['packet_path']:r for r in read(B/'upper_range092/nodes.json.gz') if r['packet_path'] in needed_packets},'jobs':{},'certs':{}}
    for st in ('upper_adaptive095','upper_adaptive103'):
        keys={s['donor_binding']['packet_key'] for s in registry.values() if s['donor_binding']['model_binding']['stage']==st};common['jobs'][st]={r['key']:r for r in read(B/st/'jobs.json.gz') if r['key'] in keys};common['certs'][st]={r['key']:r for r in read(B/st/'children_certificates.json.gz') if r['key'] in keys}
    tasks=[[] for _ in range(15)]
    for i,cid in enumerate(sorted(registry)):tasks[i%15].append((cid,registry[cid],jets[cid],coeff[cid],events[cid]))
    ctx=mp.get_context('spawn');ready=ctx.Queue();result=ctx.Queue();start=ctx.Event();stop=ctx.Event();procs=[];info={};outcomes={};fault=None
    try:
        for sid in range(15):
            q=ctx.Process(target=worker,args=(sid,tasks[sid],common,str(a.output),start,stop,ready,result));q.start();procs.append(q)
        while len(info)<15:
            try:r=ready.get(timeout=.1);info[r['sid']]=r
            except queue.Empty:pass
            if any(q.exitcode is not None for q in procs):raise RuntimeError('audit child exited before gate')
        put(a.output/'pool_registration.json',{'master_pid':os.getpid(),'workers':[info[i] for i in range(15)],'plan_sha256':a.plan_sha,'before_endpoint_checks':True});print(json.dumps({'master_pid':os.getpid(),'worker_pids':[q.pid for q in procs]}),flush=True);start.set()
        while len(outcomes)<15:
            try:r=result.get(timeout=.1);outcomes[r['sid']]=r
            except queue.Empty:pass
            if stop.is_set() or any(q.exitcode not in (None,0) for q in procs):raise RuntimeError('source binding fault;stop,preserve successful prefix')
        for q in procs:q.join();assert q.exitcode==0
        assert all(r['status']=='complete' for r in outcomes.values()) and sum(r['count'] for r in outcomes.values())==455
        link=baseline(finite,registry);put(a.output/'baseline168_link.json',link)
        for path,h in p['source_locks'].items():assert sha(path)==h,path
        for n,h in p['code_locks'].items():assert sha(H/n)==h
    except BaseException as ex:
        fault={'exception':repr(ex),'traceback':traceback.format_exc()};stop.set();start.set()
        for q in procs:
            if q.is_alive():q.terminate()
        for q in procs:
            q.join(timeout=5)
            if q.is_alive():q.kill();q.join(timeout=5)
        put(a.output/'fault.json',fault)
    finally:
        for q in procs:
            if q.is_alive():q.terminate();q.join(timeout=5)
        ready.close();result.close()
    u=resource.getrusage(resource.RUSAGE_SELF);ch=resource.getrusage(resource.RUSAGE_CHILDREN)
    receipt={'status':'PASS_all455_source_and_existing_qualification_binding' if fault is None else 'FAULT_prefix_preserved','master_pid':os.getpid(),'workers':[{'pid':q.pid,'exitcode':q.exitcode} for q in procs],
      'outcomes':outcomes,'seconds':time.monotonic()-began,'master_cpu_user':u.ru_utime,'master_cpu_system':u.ru_stime,'master_peak_rss_bytes':u.ru_maxrss,'children_cpu_user':ch.ru_utime,'children_cpu_system':ch.ru_stime,
      'plan_sha256':a.plan_sha,'fault':fault,'new_q_bandproducts_corners_H_log_gradient_calls':0,'H_enclosure_truth_shared_existing_accepted_dependency':True}
    put(a.output/'result.json',receipt);print(json.dumps({k:v for k,v in receipt.items() if k!='outcomes'}),flush=True)
    if fault:raise SystemExit(1)
if __name__=='__main__':main()
