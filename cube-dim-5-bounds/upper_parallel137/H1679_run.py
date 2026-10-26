#!/usr/bin/env python3
"""Reviewed fixed1679 residual only, fifteen explicit nonrespawning processes."""
from pathlib import Path
from fractions import Fraction as F
from datetime import datetime,timezone
import argparse,importlib.util,json,multiprocessing as mp,os,queue,sys,time,traceback
from common import sha,digest,write
import residual_certificate
import residual_run as parent_residual137
H=Path(__file__).resolve().parent;S=Path('/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch');W=S/'upper_parallel137/next1679_H_preflight';LEGACY=H.parent/'upper_continue123';TARGET=F(11,100000);WORKERS=15
def legacy_load(package):
 for p,h in package['legacy_code_locks'].items():assert sha(p)==h
 sys.path.insert(0,str(LEGACY));spec=importlib.util.spec_from_file_location('locked123_dispatch_helpers',LEGACY/'run.py');m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);assert m.TARGET==TARGET;return m
def worker_loop(package,tasks,returned,ready):
 try:
  old=legacy_load(package);old.worker.initialize(package['table_path'],package['table_sha256'],str(TARGET),ready)
  while True:
   task=tasks.get()
   if task is None:return
   got=old.worker.evaluate_job(task);returned.put({'ordinal':task['batch_ordinal'],'returned':got})
 except BaseException as exc:returned.put({'fatal_worker_pid':os.getpid(),'exception':repr(exc),'traceback':traceback.format_exc()})
def source_seed(package):
 seed=parent_residual137.source_seed(package);last=json.load(open(S/'upper_parallel137/residual_preflight/future_batch0002/state.json'));seed['passed']+=last['new_passed'];seed['splits']+=last['new_splits']
 for p,r in last['new_cache_by_epsilon']['1/128'].items():assert p not in seed['cache_by_epsilon']['1/128'];seed['cache_by_epsilon']['1/128'][p]=r
 state=json.load(open(S/'upper_parallel137/next2350_preflight/run001/master_attempt0001_state.json'))
 for k,c in state['qualified_children'].items():seed['passed'].append({'node':c['child_node'],'certificate':c,'proof_source':{'path':str(S/'upper_parallel137/next2350_preflight/run001/master_attempt0001_state.json'),'qualified_child_index':k}})
 seed['queue']=state['pending_children'];assert len(seed['passed'])==109405 and len(seed['splits'])==58724 and len(seed['queue'])==1679 and len(seed['cache_by_epsilon']['1/128'])==14801 and len(seed['cache_by_epsilon']['1/64'])==184796;return seed
def main():
 began=time.monotonic();ap=argparse.ArgumentParser();ap.add_argument('--dispatch-reviewed',action='store_true');ap.add_argument('--batch',type=int,choices=(1,),required=True);ap.add_argument('--target',required=True);ap.add_argument('--reviewed-manifest-sha',required=True);a=ap.parse_args();assert a.dispatch_reviewed and F(a.target)==TARGET
 manifest=H/'next1679_H_prepared_manifest.json';assert sha(manifest)==a.reviewed_manifest_sha;package=json.load(open(manifest))
 for n,h in package['repository_files'].items():assert sha(H/n)==h
 for p,h in {**package['source_locks'],**package['external_files'],**package['legacy_code_locks']}.items():assert sha(p)==h
 old=legacy_load(package);seed=source_seed(package);plan=json.load(open(W/'summary.json'));selection=plan['finite_batches'][a.batch-1];tasks=json.load(open(selection['path']));assert sha(selection['path'])==selection['sha256'] and len(tasks)==selection['count'] and len(tasks)<=4096 and all(not t['cached'] for t in tasks)
 current=W/'CURRENT.json';parent=None
 if a.batch==1:
  assert not current.exists();ordered=json.load(open(W/'tasks.json'));assert {tuple(t['node']) for t in ordered}=={tuple(n) for n in seed['queue']} and len(ordered)==1679
  ov={'seed_source_manifest':str(H/'next2350_results_manifest.json'),'seed_source_manifest_sha256':package['result_manifest_sha256'],'target':str(TARGET),'pending_evaluation_epsilon':'1/128','lower_bound_object':'G','new_cache_by_epsilon':{'1/128':{}},'new_passed':[],'new_splits':[],'queue':[t['node'] for t in ordered],'inflight':{},'obstructions_H128':[],'faults':[],'counts':{'new_anchor_starts':0,'new_complete':0,'kernel_segment_calls':0,'cache_hits':0},'status':'initialized'}
 else:
  parent=json.load(open(current));assert parent['batch']==1 and sha(parent['state'])==parent['state_sha256'] and sha(parent['journal'])==parent['journal_sha256'];ov=json.load(open(parent['state']));assert ov['status']=='finite_original_prefix_complete' and not ov['faults'] and not ov['obstructions_H128'] and not ov['inflight']
 assert ov['queue'][:len(tasks)]==[t['node'] for t in tasks] and all(t['batch']==a.batch and t['batch_ordinal']==i and t['target']==str(TARGET) and t['epsilon']=='1/128' for i,t in enumerate(tasks))
 assert all(t['pair_cache_key'][1] not in seed['cache_by_epsilon']['1/128'] and t['pair_cache_key'][1] not in ov['new_cache_by_epsilon']['1/128'] for t in tasks)
 bdir=W/f'future_batch{a.batch:04d}';assert not bdir.exists() and not (W/'STOP').exists() and all(not Path(t['point_checkpoint']).exists() for t in tasks);bdir.mkdir();oldcache=digest(seed['cache_by_epsilon']);oldpassed=digest(seed['passed']);before=dict(ov['counts']);beforepass=len(ov['new_passed']);beforesplit=len(ov['new_splits']);initialtree=old.tree(seed,ov)
 ov.update({'batch':a.batch,'next_merge':0,'journal_seq':0,'journal_hash':'0'*64,'status':'registered_finite_residual'});reg={'master_pid':os.getpid(),'registered_utc':datetime.now(timezone.utc).isoformat(),'batch':a.batch,'reviewed_manifest_sha256':a.reviewed_manifest_sha,'target':str(TARGET),'epsilon':'1/128','selected':len(tasks),'maximum_new_anchor_starts':len(tasks),'tasks_sha256':selection['sha256'],'source_cache_canonical_sha256':oldcache,'source_passed_canonical_sha256':oldpassed,'parent_checkpoint':parent,'initial_tree':initialtree,'workers':15,'wall_deadline':None,'peak_memory_estimate_GiB':[4,8],'memory_estimate_not_cap':True,'no_newborn_child_dispatch':True};write(bdir/'registration.json',reg);write(bdir/'initial_overlay.json',ov);journal=open(bdir/'journal.jsonl','ab')
 def emit(kind,body):
  event={'seq':ov['journal_seq']+1,'previous':ov['journal_hash'],'kind':kind,'body':body};event['hash']=digest(event);journal.write((json.dumps(event,sort_keys=True,separators=(',',':'))+'\n').encode());journal.flush();os.fsync(journal.fileno());old.apply(ov,event,seed['cache_by_epsilon'])
  if ov['journal_seq']%128==0 or kind in ('coordinator_fault','finish'):write(bdir/'state.json',ov)
 def receive(item):
  if 'fatal_worker_pid' in item:raise RuntimeError(item)
  i=item['ordinal'];assert i in pending;got=item['returned'];assert got['epsilon']=='1/128' and got['target']==str(TARGET) and got['pair_cache_key']==tasks[i]['pair_cache_key'];emit('received',item);pending.remove(i)
  if got['kind']!='complete':write(W/'STOP',{'ordinal':i,'kind':got['kind'],'target':str(TARGET),'epsilon':'1/128','scope':'fixedH128 only;or diagnostic;not G/C5 counterexample'});return True
  return False
 def merge_ready():
  while str(ov['next_merge']) in ov['inflight'] and 'returned' in ov['inflight'][str(ov['next_merge'])]:
   i=ov['next_merge'];ent=ov['inflight'][str(i)];body={'ordinal':i}
   if ent['returned']['kind']=='complete':body['certificate']=residual_certificate.make(old.physical(ent['task']['node'],seed['roots']),ent['returned']['record'],TARGET,ent['task'],old)
   emit('merge',body)
 ctx=mp.get_context('spawn');tq=ctx.Queue();rq=ctx.Queue();ready=ctx.Queue();processes=[];pending=set();nextdispatch=0;stop=False;sentinels=False
 try:
  for i in range(WORKERS):p=ctx.Process(target=worker_loop,args=(package,tq,rq,ready),name=f'H128-residual137-{i:02d}');p.start();processes.append(p)
  infos=[]
  while len(infos)<WORKERS:
   try:infos.append(ready.get(timeout=.1))
   except queue.Empty:pass
   if any(p.exitcode is not None for p in processes):raise RuntimeError('worker exited before all PID registration')
  assert len({r['pid'] for r in infos})==15 and all(r['target']==str(TARGET) and r['epsilon']=='1/128' and r['kernel_sha256']==package['kernel_sha256'] for r in infos);write(bdir/'pool_registration.json',{'master_pid':os.getpid(),'workers':infos,'before_any_H_dispatch':True,'registration_sha256':sha(bdir/'registration.json')});print(json.dumps({'event':'15_actual_workers_registered_before_H','batch':a.batch,'master_pid':os.getpid(),'worker_pids':[r['pid'] for r in infos],'selected':len(tasks)}),flush=True)
  while pending or nextdispatch<len(tasks) and not stop:
   stop|=(W/'STOP').exists()
   while not stop and nextdispatch<len(tasks) and len(pending)<15:
    task=tasks[nextdispatch];assert not Path(task['point_checkpoint']).exists();emit('dispatch',{'task':task});pending.add(nextdispatch);tq.put(task);nextdispatch+=1
   try:stop|=receive(rq.get(timeout=.05))
   except queue.Empty:pass
   merge_ready()
   if any(p.exitcode is not None for p in processes):raise RuntimeError('worker death/no replacement;preserve successfulpoint/coneprefix')
  for _ in processes:tq.put(None)
  sentinels=True
 except BaseException as exc:
  stop=True;write(W/'STOP',{'kind':'coordinator_fault','exception':repr(exc),'target':str(TARGET),'epsilon':'1/128'});emit('coordinator_fault',{'exception':repr(exc),'traceback':traceback.format_exc(),'action':'preserve point/cone/checkpoint prefixes;no blind restart'})
 finally:
  if not sentinels:
   for _ in processes:tq.put(None)
  # Drain already-dispatched finite jobs while stopping; avoid Queue feeder/join deadlock.
  while any(p.is_alive() for p in processes):
   try:item=rq.get(timeout=.05)
   except queue.Empty:continue
   if 'fatal_worker_pid' in item:emit('coordinator_fault',item)
   elif item['ordinal'] in pending:
    try:receive(item);merge_ready()
    except BaseException as exc:emit('coordinator_fault',{'exception':repr(exc),'scope':'shutdown receipt;raw point checkpoints preserved'})
  for p in processes:p.join()
  tq.close();rq.close();ready.close();tq.join_thread();rq.join_thread();ready.join_thread()
 evidence=[{'pid':p.pid,'exitcode':p.exitcode,'alive':p.is_alive()} for p in processes];status='H128_g110_target_obstruction' if ov['obstructions_H128'] else 'fault_or_partial' if ov['faults'] else 'finite_original_prefix_complete' if ov['next_merge']==len(tasks) and not ov['inflight'] else 'explicit_stop_partial';emit('finish',{'status':status,'pool_evidence':evidence,'wall_seconds':time.monotonic()-began});journal.close();t=old.tree(seed,ov);assert digest(seed['cache_by_epsilon'])==oldcache and digest(seed['passed'])==oldpassed
 for p,h in package['source_locks'].items():assert sha(p)==h
 points=[]
 for task in tasks[:nextdispatch]:
  p=Path(task['point_checkpoint'])
  if p.exists():cp=json.load(open(p));assert cp['point']==task['point'] and cp['epsilon']=='1/128' and cp['target']==str(TARGET) and cp['kernel_sha256']==package['kernel_sha256'];points.append({'ordinal':task['batch_ordinal'],'global_ordinal':task['ordinal'],'node':task['node'],'epsilon':'1/128','target':str(TARGET),'point':task['point'],'point_checkpoint':str(p),'point_checkpoint_sha256':sha(p),'stage':cp['stage'],'worker_pid':cp['pid'],'partial':cp['partial'],'record':cp.get('record')})
 write(bdir/'point_manifest.json',{'records':points,'epsilon':'1/128','target':str(TARGET),'kernel_sha256':package['kernel_sha256']});summary={'status':status,'batch':a.batch,'master_pid':os.getpid(),'selected':len(tasks),'new_counts':{k:ov['counts'][k]-before[k] for k in before},'new_passed':len(ov['new_passed'])-beforepass,'new_splits':len(ov['new_splits'])-beforesplit,'tree':t,'faults':len(ov['faults']),'obstructions':len(ov['obstructions_H128']),'inflight':len(ov['inflight']),'queue':len(ov['queue']),'pool_evidence':evidence,'state_sha256':sha(bdir/'state.json'),'journal_sha256':sha(bdir/'journal.jsonl'),'wall_seconds':time.monotonic()-began,'globalg110_qualified':False,'S_credit':0,'no_automatic_next_batch':True};write(bdir/'summary.json',summary);write(bdir/'manifest.json',{'reviewed_manifest_sha256':a.reviewed_manifest_sha,'files':{p.name:sha(p) for p in bdir.iterdir() if p.is_file() and p.name!='manifest.json'}})
 if status=='finite_original_prefix_complete' and not ov['faults'] and not ov['obstructions_H128'] and not ov['inflight']:write(current,{'batch':a.batch,'state':str(bdir/'state.json'),'state_sha256':sha(bdir/'state.json'),'journal':str(bdir/'journal.jsonl'),'journal_sha256':sha(bdir/'journal.jsonl')})
 print(json.dumps(summary),flush=True)
if __name__=='__main__':main()
