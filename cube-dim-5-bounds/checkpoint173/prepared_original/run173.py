"""Future one455source + fixed35872pair program; no automatic H or children."""
from pathlib import Path
from fractions import Fraction as F
from collections import defaultdict
import argparse,json,gzip,multiprocessing as mp,os,queue,time,traceback,resource
from util173 import sha,digest,read,atomic,create
from worker173 import execute
H=Path(__file__).resolve().parent;W=Path('/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/upper_bandcurvature173')
def run_stage(phase,specs,O,reviewed,runtime,resume,attempt):
 ctx=mp.get_context('spawn');ready=ctx.Queue();result=ctx.Queue();start=ctx.Event();stop=ctx.Event();processes=[];infos={};summaries={};fault=None
 try:
  for spec in specs:
   p=ctx.Process(target=execute,args=(phase,spec,str(O),reviewed,runtime,ready,result,start,stop,resume),name=f'173-{phase}{spec["shard_id"]:02d}');p.start();processes.append(p)
  while len(infos)<15 and not stop.is_set():
   try:r=ready.get(timeout=.1);assert r['shard_id'] not in infos;infos[r['shard_id']]=r
   except queue.Empty:pass
   if any(p.exitcode is not None for p in processes):raise RuntimeError('child exited before gate')
  assert not stop.is_set() and len(infos)==15 and len({r['pid'] for r in infos.values()})==15;atomic(O/f'actual_attempt{attempt:04d}_{phase}_pool_registration.json',{'master_pid':os.getpid(),'workers':[infos[i] for i in range(15)],'phase':phase,'before_new_band_weights_or_corners':True});print(json.dumps({'event':'registered15_before_gate','phase':phase,'master_pid':os.getpid(),'worker_pids':[infos[i]['pid'] for i in range(15)]}),flush=True);start.set()
  while len(summaries)<15:
   try:r=result.get(timeout=.1);assert r['shard_id'] not in summaries;summaries[r['shard_id']]=r
   except queue.Empty:pass
   if any(r['status'] not in ('complete','stopped_successful_prefix') for r in summaries.values()):stop.set()
   for i,p in enumerate(processes):
    if p.exitcode not in (None,0) and i not in summaries:stop.set();summaries[i]={'status':'process_fault','shard_id':i,'pid':p.pid}
   if all(p.exitcode is not None for p in processes) and len(summaries)<15:
    while True:
     try:r=result.get_nowait();summaries.setdefault(r['shard_id'],r)
     except queue.Empty:break
    for i,p in enumerate(processes):summaries.setdefault(i,{'status':'missing_result_fault','shard_id':i,'pid':p.pid})
 except BaseException as e:fault={'exception':repr(e),'traceback':traceback.format_exc()};stop.set();start.set()
 finally:
  if stop.is_set():start.set()
  while any(p.is_alive() for p in processes):
   try:r=result.get(timeout=.1);summaries.setdefault(r['shard_id'],r)
   except queue.Empty:pass
   for p in processes:p.join(timeout=0)
  for p in processes:p.join()
  while True:
   try:r=result.get_nowait();summaries.setdefault(r['shard_id'],r)
   except queue.Empty:break
  ready.close();result.close();ready.join_thread();result.join_thread()
 evidence=[{'pid':p.pid,'exitcode':p.exitcode,'alive':p.is_alive()} for p in processes];complete=not fault and len(summaries)==15 and all(r['status']=='complete' for r in summaries.values()) and all(p.exitcode==0 for p in processes);receipt={'phase':phase,'master_pid':os.getpid(),'complete':complete,'workers':evidence,'summaries':summaries,'fault':fault};atomic(O/f'actual_attempt{attempt:04d}_{phase}_exit_receipt.json',receipt);return receipt

def original_jets_once(plan,sources,O):
 byjet={s['original_g151_jet_ref']['jet_id']:cid for cid,s in sources.items()};assert len(byjet)==455;found={}
 for part in plan['original_g151_jet_parts']:
  assert sha(part['path'])==part['sha256']
  with gzip.open(part['path'],'rt') as f:
   for r in map(json.loads,f):
    if r['jet_id'] in byjet:
     h=r['minimal_projection_sha256'];payload={k:v for k,v in r.items() if k!='minimal_projection_sha256'};assert digest(payload)==h;cid=byjet[r['jet_id']];assert cid not in found and r['model']==plan['g151_model'] and r['epsilon']=='1/64';assert r['raw_record_sha256']==sources[cid]['old155']['raw_record_sha256'] and r['point']==sources[cid]['old155']['point'];found[cid]=r
 assert set(found)==set(sources);runtime=[]
 for spec in plan['phase1_mass_shards']:
  ids=read(spec['path'])['source_ids'];path=O/f'original_g_input_shard{spec["shard_id"]:02d}.json';data={cid:found[cid] for cid in ids}
  if path.exists():assert read(path)==data
  else:create(path,data)
  runtime.append({'g_input_path':str(path),'g_input_sha256':sha(path)})
 return runtime

def assemble_full_domain(plan,O,qregistry,sources,C):
 baseupdates={r['leaf_index']:r for r in read(plan['base158_updates_path'])};basecandidate=read(plan['base158_candidate_path']);assert basecandidate['full115_passed']==110081 and basecandidate['target115_pending']==1120;newbest={};receivers={};evaluated=firstskip=oldskip=0
 for spec in plan['phase2_corner_shards']:
  data=read(spec['input_path']);receivers.update({int(k):r for k,r in data['receivers'].items()});out=O/'corners'/f'shard{spec["shard_id"]:02d}';best=read(out/'receiver_best.json');ss=read(out/'summary.json');assert digest(best)==ss['best_sha256'];evaluated+=ss['counts']['evaluated'];firstskip+=ss['counts']['skip_after_firstPASS'];oldskip+=ss['counts']['skip_oldpair_q_unchanged']
  for k,b in best.items():assert int(k) not in newbest;newbest[int(k)]=b
 assert set(newbest)==set(receivers)==set(range(1120)) and evaluated+firstskip+oldskip==35872;by_leaf={receivers[i]['leaf_index']:(receivers[i],b) for i,b in newbest.items()};assert len(by_leaf)==1120;seen=set();roots=defaultdict(list);minimum=None;base_min=None;pending=[];updates=[];newpass=improved=0
 assert sha(plan['source155_margin_part'])==plan['source155_margin_sha256']
 with gzip.open(plan['source155_margin_part'],'rt') as f:
  for row in map(json.loads,f):
   i=row['leaf_index'];assert i not in seen;seen.add(i);root,path=row['node'];roots[root].append(path);v=F(row['margin110_new']);source={'kind':'original155','canonical_id':row['canonical_id']}
   if i in baseupdates:
    b=baseupdates[i];assert b['node']==row['node'] and F(b['best']['best_margin110'])>=v;v=F(b['best']['best_margin110']);source={'kind':b['best']['best_source_kind'],'canonical_id':b['best']['best_canonical_id']}
   assert v>0;base_min=v if base_min is None else min(base_min,v)
   if i in by_leaf:
    r,b=by_leaf[i];assert r['node']==row['node'] and r['base_best158']==baseupdates[i]['best'] and F(r['base_best158']['best_margin110'])==v;assert F(b['best_margin110'])>=v and F(b['best_margin115'])==F(b['best_margin110'])-C;improved+=F(b['best_margin110'])>v;v=F(b['best_margin110']);source={'kind':b['best_source_kind'],'canonical_id':b['best_canonical_id']};newpass+=b['first_pass115'] is not None;updates.append({'leaf_index':i,'node':row['node'],'best':b})
    if b['first_pass115'] is None:pending.append({'leaf_index':i,'node':row['node'],'best_canonical_id':b['best_canonical_id'],'best_margin110':b['best_margin110'],'best_margin115':b['best_margin115']})
   if i not in by_leaf:assert v-C>0
   key=(v,i);minimum=(key,row['node'],source) if minimum is None or key<minimum[0] else minimum
 assert seen==set(range(111201)) and len(roots)==52360 and newpass+len(pending)==1120 and F(basecandidate['minimum_margin110_updated'])==base_min and minimum[0][0]>=base_min
 for root in range(52360):
  paths=sorted(roots[root]);assert all(not b.startswith(a) for a,b in zip(paths,paths[1:]));d=max(map(len,paths));assert sum(1<<(d-len(p)) for p in paths)==1<<d
 atomic(O/'receiver_best_updates.json',sorted(updates,key=lambda r:r['leaf_index']));atomic(O/'target115_pending.json',sorted(pending,key=lambda r:r['leaf_index']));candidate={'full111201_closed_domain_retained':True,'unchanged158_leaves':110081,'updated_receivers':1120,'best_rule':'max(base158,all actually tested new single-source wholebox margins including failed)','base158_minimum_margin110':str(base_min),'minimum_margin110_updated':str(minimum[0][0]),'minimum_ge_base158':True,'witness_leaf_index':minimum[0][1],'witness_node':minimum[1],'witness_source':minimum[2],'minimum_margin115':str(minimum[0][0]-C),'candidate_linear_floor':str(F(11,100000)*(1+minimum[0][0])),'new115_passes':newpass,'target115_pending':len(pending),'full115_passed':110081+newpass,'improved_receivers':improved,'source_q_strict_decreases':sum(r['certificate']['q_strictly_decreased'] for r in qregistry.values()),'evaluated_pairs':evaluated,'firstPASS_skips':firstskip,'oldpair_q_unchanged_skips':oldskip,'main_P8_down192_separate':True,'formal_g_C5_update':False};atomic(O/'complete_domain_candidate.json',candidate);return candidate

def main():
 began=time.monotonic();ap=argparse.ArgumentParser();ap.add_argument('--dispatch-reviewed',action='store_true');ap.add_argument('--reviewed-manifest-sha',required=True);ap.add_argument('--authorization',required=True);ap.add_argument('--authorization-sha',required=True);ap.add_argument('--resume-reviewed',action='store_true');a=ap.parse_args();assert a.dispatch_reviewed;mpath=H/'prepared_manifest.json';assert sha(mpath)==a.reviewed_manifest_sha;m=read(mpath)
 for name,h in m['repository_files'].items():assert sha(H/name)==h
 for p,h in {**m['source_locks'],**m['prepared_external_files']}.items():assert sha(p)==h
 assert sha(a.authorization)==a.authorization_sha;auth=read(a.authorization);assert auth['reviewed_manifest_sha256']==a.reviewed_manifest_sha and auth['sources']==455 and auth['fixed_pairs']==35872 and auth['workers']==15
 for key in ('qualified168','qualified172_174','qualified_donor_endpoint_provenance'):
  gate=auth[key];assert gate['accepted'] is True and sha(gate['path'])==gate['sha256']
 plan=read(W/'finite_plan.json');assert plan['weighted_band_terms']==22295 and plan['sources']==455 and plan['fixed_pairs']==35872 and plan['max_concurrent_workers']==15;O=W/'run001'
 if O.exists():assert a.resume_reviewed
 else:assert not a.resume_reviewed;O.mkdir()
 assert not (O/'complete_domain_candidate.json').exists();attempt=1
 while (O/f'actual_attempt{attempt:04d}_registration.json').exists():attempt+=1
 atomic(O/f'actual_attempt{attempt:04d}_registration.json',{'master_pid':os.getpid(),'reviewed_manifest_sha256':a.reviewed_manifest_sha,'authorization_sha256':a.authorization_sha,'sources':455,'weighted_band_terms':22295,'fixed_pairs':35872,'receivers':1120,'max_concurrent_workers':15,'phase1_unique_sources_then_phase2_receiver_owned':True,'no_new_H_log_gradient_mean_split':True,'wall_deadline':None});fault=None;candidate=None;exits=[]
 try:
  assert sha(plan['source_registry_path'])==plan['source_registry_sha256'];sources=read(plan['source_registry_path']);ginputs=original_jets_once(plan,sources,O);runtime={k:plan[k] for k in ('source_registry_path','source_registry_sha256','old_common_log_path','old_common_log_sha256','time_endpoints')}
  # Per-shard original g input metadata is dispatched only to its own source worker.
  massspec=plan['phase1_mass_shards']
  # Runtime wrapper takes shard-specific g input; no mass or corner function has run yet.
  mass_runtime={**runtime,'g_input_by_shard':ginputs};r=run_stage('mass',massspec,O,a.reviewed_manifest_sha,mass_runtime,a.resume_reviewed,attempt);exits.append(r);assert r['complete'],'mass fault;preserve all source prefixes and stop'
  qregistry={}
  for spec in massspec:
   d=read(O/'mass'/f'shard{spec["shard_id"]:02d}/source_certificates.json');assert set(d)==set(read(spec['path'])['source_ids'])
   for cid,q in d.items():assert cid not in qregistry and q['certificate']['canonical_id']==cid and len(q['certificate']['bands'])==49 and 0<F(q['certificate']['q_new'])<=F(q['certificate']['oldqualified155q']);qregistry[cid]=q
  assert set(qregistry)==set(sources) and len(qregistry)==455
  for source in sources.values():
   md=source['source171_metadata']['new_mass_donor'];binding=source['donor_binding'];assert sha(md['packet_path'])==md['packet_sha256']
   if binding['kind']=='legacy092_compact':assert sha(binding['original_raw_path'])==binding['original_raw_sha256']
  qp=O/'source_q_registry.json'
  if qp.exists():assert read(qp)==qregistry
  else:create(qp,qregistry)
  C=F(read(plan['old_common_log_path'])['upper']);corner_runtime={**runtime,'q_registry_path':str(qp),'q_registry_sha256':sha(qp)};r=run_stage('corners',plan['phase2_corner_shards'],O,a.reviewed_manifest_sha,corner_runtime,a.resume_reviewed,attempt);exits.append(r);assert r['complete'],'corner fault;preserve prefixes and stop';candidate=assemble_full_domain(plan,O,qregistry,sources,C)
 except BaseException as e:fault={'exception':repr(e),'traceback':traceback.format_exc(),'successful_prefixes_preserved':True,'automatic_restart':False};atomic(O/f'fault_attempt{attempt:04d}.json',fault)
 u=resource.getrusage(resource.RUSAGE_SELF);c=resource.getrusage(resource.RUSAGE_CHILDREN);summary={'status':'complete_fixed455source_35872pair_candidate' if candidate is not None and fault is None else 'fault_prefix_preserved','master_pid':os.getpid(),'phase_exit_receipts':[{'phase':e['phase'],'complete':e['complete'],'workers':e['workers']} for e in exits],'candidate':candidate,'fault':fault,'new_H_log_gradient_mean_split_calls':0,'weighted_band_terms_if_complete':22295 if candidate else None,'master_cpu_user_seconds':u.ru_utime,'master_cpu_system_seconds':u.ru_stime,'children_cpu_user_seconds':c.ru_utime,'children_cpu_system_seconds':c.ru_stime,'master_peak_rss_bytes':u.ru_maxrss,'wall_seconds':time.monotonic()-began,'formal_acceptance':False};atomic(O/f'actual_attempt{attempt:04d}_summary.json',summary);print(json.dumps(summary),flush=True)
if __name__=='__main__':main()
