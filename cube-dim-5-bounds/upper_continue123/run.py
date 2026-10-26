#!/usr/bin/env python3
"""Reviewed seed123 -> exactly two finite ORIGINAL-frontier input batches.

Worker/certificate/cache core is byte-identical to107v2. Geometry children
remain pending and are never selected in these frozen4096+2141 batches.
"""
from pathlib import Path
from fractions import Fraction as F
from collections import defaultdict
from datetime import datetime,timezone
import argparse,hashlib,json,multiprocessing as mp,os,sys,time
from binding import KERNEL_PATH,KERNEL_SHA,MODEL
from cache_model import validate_record,insert_new
import worker,certificate,reuse
H=Path(__file__).resolve().parent;W=Path('/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/upper_continue123')
SEED=W/'seed.json';TABLE=Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube5-astra-upper/cube-dim-5-bounds/astra_vertex017/boundary_triangulation.json.gz');TSHA='8dd2e9895d7465ac599d43feb1b5099506f4e43d1248215eba0a77e11f4d13f9'
TARGET=F(11,100000);WORKERS=8
ENC=json.JSONEncoder(sort_keys=True,separators=(',',':'))
def sha(p):
 h=hashlib.sha256()
 with open(p,'rb') as f:
  while b:=f.read(8*1024*1024):h.update(b)
 return h.hexdigest()
def digest(x):
 h=hashlib.sha256()
 for b in ENC.iterencode(x):h.update(b.encode())
 return h.hexdigest()
def write(p,x):
 p=Path(p);t=p.with_suffix(p.suffix+'.tmp')
 with open(t,'w') as f:
  for b in ENC.iterencode(x):f.write(b)
  f.write('\n');f.flush();os.fsync(f.fileno())
 os.replace(t,p)
def physical(node,roots):
 lo=[F(x,32) for x in roots[node[0]]];hi=[x+F(1,32) for x in lo]
 for bit in node[1]:
  j=max(range(4),key=lambda j:hi[j]-lo[j]);mid=(lo[j]+hi[j])/2
  if bit=='0':hi[j]=mid
  else:assert bit=='1';lo[j]=mid
 return lo,hi
def tree(seed,ov):
 by=defaultdict(list);seen=set();inner=set()
 for n in [r['node'] for r in seed['passed']+ov['new_passed']]+ov['queue']+[e['task']['node'] for e in ov['inflight'].values()]:
  r,p=n;assert (r,p) not in seen;seen.add((r,p));by[r].append(p)
 assert len(by)==52360
 for r in range(52360):
  ps=sorted(by[r]);assert ps and all(not b.startswith(a) for a,b in zip(ps,ps[1:]));d=max(map(len,ps));assert sum(1<<(d-len(p)) for p in ps)==1<<d
 for row in seed['splits']+ov['new_splits']:
  r,p=row['node'];assert (r,p) not in inner;inner.add((r,p));assert row['axis']==len(p)%4 and row['children']==[[r,p+'0'],[r,p+'1']]
 assert not seen&inner
 for r,p in seen|inner:
  if p:assert (r,p[:-1]) in inner
 assert len(seen)==52360+len(inner)
 return {'exact_closed_Kraft_roots':52360,'passed':len(seed['passed'])+len(ov['new_passed']),'queue_and_inflight':len(ov['queue'])+len(ov['inflight']),'splits':len(inner),'active_leaves':len(seen),'antichain':True}
def cached_result(task,seed):
 """Normalize only the per-job view, never rewrite the old cache object."""
 raw=seed['cache_by_epsilon']['1/128'][task['pair_cache_key'][1]];assert digest(raw)==task['cache_source_record_sha256']
 view={**raw,'model':MODEL,'epsilon':'1/128','qualification_target':task['target'],'pair_cache_key':task['pair_cache_key']};validate_record(view,task['target'])
 kind='H128_g110_target_obstruction' if F(raw['values'][0][1])<F(task['target']) else 'complete'
 cp={'stage':'cached_fulljet_consumed','epsilon':'1/128','target':task['target'],'pair_cache_key':task['pair_cache_key'],'point':task['point'],'kernel_sha256':KERNEL_SHA,'record':view,'source_raw_record_sha256':digest(raw),'partial':None,'calls':0,'attempts':[],'source_attempts':raw.get('attempts',[]),'pid':os.getpid()}
 assert not Path(task['point_checkpoint']).exists();write(Path(task['point_checkpoint']),cp)
 return {'kind':kind,'epsilon':'1/128','target':task['target'],'pair_cache_key':task['pair_cache_key'],'record':view,'source_raw_record_sha256':digest(raw),'pid':os.getpid(),'calls':0,'attempts':[],'cached_fulljet':True,'scope':'existing complete EPS128 jet,zero new H;per-job metadata view only'}
def apply(ov,event,seedcache):
 b=event['body'];kind=event['kind'];c=ov['counts']
 if kind=='dispatch':
  task=b['task'];assert task['target']==ov['target'] and task['epsilon']=='1/128' and ov['queue'].pop(0)==task['node'];i=str(task['batch_ordinal']);assert i not in ov['inflight'];ov['inflight'][i]={'task':task};c['new_anchor_starts']+=not task['cached']
 elif kind=='received':
  ent=ov['inflight'][str(b['ordinal'])];got=b['returned'];assert got['target']==ov['target'] and got['epsilon']=='1/128' and got['pair_cache_key']==ent['task']['pair_cache_key']
  if 'record' in got:
   validate_record(got['record'],ov['target']);assert got['record']['point']==ent['task']['point'] and got['record']['pair_cache_key']==ent['task']['pair_cache_key']
  assert 'returned' not in ent;ent['returned']=got;c['kernel_segment_calls']+=got.get('calls',0)
 elif kind=='merge':
  i=b['ordinal'];assert i==ov['next_merge'];ov['next_merge']+=1;ent=ov['inflight'][str(i)];got=ent['returned'];task=ent['task'];node=task['node']
  if got['kind']=='complete':
   rec=got['record'];assert rec['point']==task['point'] and rec['pair_cache_key']==task['pair_cache_key']
   if task['cached']:
    raw=seedcache['1/128'][task['pair_cache_key'][1]];assert digest(raw)==task['cache_source_record_sha256']==got['source_raw_record_sha256'] and raw['values']==rec['values'] and raw['point']==rec['point'];c['cache_hits']+=1
   else:insert_new(ov['new_cache_by_epsilon']['1/128'],seedcache,rec,ov['target']);c['new_complete']+=1
   cert=b['certificate'];assert cert['target']==ov['target'] and cert['leaf_epsilon']=='1/128' and cert['point_key']==rec['pair_cache_key'][1]
   if cert['passed']:ov['new_passed'].append({'node':node,'leaf_epsilon':'1/128','lower_bound_object':'G','target':ov['target'],'certificate':cert})
   else:
    children=[[node[0],node[1]+'0'],[node[0],node[1]+'1']];ov['new_splits'].append({'node':node,'axis':len(node[1])%4,'children':children,'leaf_epsilon':'1/128','target':ov['target'],'failed_certificate':cert});ov['queue'].extend(children)
   del ov['inflight'][str(i)]
  elif got['kind']=='H128_g110_target_obstruction':ov['obstructions_H128'].append({'node':node,'epsilon':'1/128','target':ov['target'],'returned':got,'scope':'fixedH128 only,not G/C5 refutation'})
  else:ov['faults'].append({'node':node,'returned':got})
 elif kind=='coordinator_fault':ov['faults'].append(b)
 elif kind=='finish':ov.update(b)
 else:raise ValueError(kind)
 ov['journal_seq']=event['seq'];ov['journal_hash']=event['hash']
def main():
 began=time.monotonic();ap=argparse.ArgumentParser();ap.add_argument('--batch',type=int,required=True,choices=(1,2));ap.add_argument('--target',required=True);ap.add_argument('--reviewed-manifest-sha',required=True);a=ap.parse_args();target=F(a.target);assert target==TARGET
 manifest=H/'production_prepared_manifest.json';assert sha(manifest)==a.reviewed_manifest_sha;package=json.load(open(manifest))
 for name,h in package['repository_files'].items():assert sha(H/name)==h
 for path,h in package['external_files'].items():assert sha(path)==h
 for path,h in package['source_locks'].items():assert sha(path)==h
 plan=json.load(open(W/'finite_plan.json'));assert sha(SEED)==plan['seed_sha256'] and not SEED.stat().st_mode&0o222;seed=json.load(open(SEED));assert len(seed['passed'])==98007 and len(seed['queue'])==6237 and len(seed['splits'])==51884
 assert len(seed['cache_by_epsilon']['1/64'])==184796 and len(seed['cache_by_epsilon']['1/128'])==4290 and seed['target']==str(target) and not seed['inflight'] and not seed['resume_jobs']
 selection=plan['batches'][a.batch-1];assert selection['batch']==a.batch and selection['worker_count']==WORKERS and selection['count'] in (4096,2141);tasks=json.load(open(selection['tasks_path']));assert sha(selection['tasks_path'])==selection['tasks_sha256'] and len(tasks)==selection['count'] and digest([r['node'] for r in tasks])==selection['selected_nodes_sha256']
 parent=None;current=W/'CURRENT.json'
 if a.batch==1:
  assert not current.exists();ordered=json.load(open(W/'planned_tasks.json'));assert {tuple(r['node']) for r in ordered}=={tuple(n) for n in seed['queue']} and len(ordered)==6237
  ov={'seed':str(SEED),'seed_sha256':plan['seed_sha256'],'target':str(target),'pending_evaluation_epsilon':'1/128','lower_bound_object':'G','new_cache_by_epsilon':{'1/128':{}},'new_passed':[],'new_splits':[],'queue':[r['node'] for r in ordered],'inflight':{},'obstructions_H128':[],'faults':[],'counts':{'new_anchor_starts':0,'new_complete':0,'kernel_segment_calls':0,'cache_hits':0},'status':'initialized'}
 else:
  parent=json.load(open(current));assert parent['batch']==1 and sha(parent['state'])==parent['state_sha256'];ov=json.load(open(parent['state']));assert ov['status']=='finite_original_prefix_complete' and not ov['inflight'] and not ov['faults'] and not ov['obstructions_H128']
  assert ov['seed_sha256']==plan['seed_sha256'] and sha(parent['journal'])==parent['journal_sha256']
 assert ov['queue'][:len(tasks)]==[r['node'] for r in tasks] and all(r['batch']==a.batch and r['batch_ordinal']==i for i,r in enumerate(tasks))
 assert all(r['cached'] or r['pair_cache_key'][1] not in seed['cache_by_epsilon']['1/128'] and r['pair_cache_key'][1] not in ov['new_cache_by_epsilon']['1/128'] for r in tasks)
 assert not (W/'STOP').exists();bdir=W/f'batch{a.batch:04d}';assert not bdir.exists() and all(not Path(r['point_checkpoint']).exists() for r in tasks);bdir.mkdir(exist_ok=False)
 oldcachehash=digest(seed['cache_by_epsilon']);oldpassedhash=digest(seed['passed']);oldcounts=dict(ov['counts']);oldlens={k:len(ov[k]) for k in ('new_passed','new_splits')};initialtree=tree(seed,ov)
 ov.update({'next_merge':0,'journal_seq':0,'journal_hash':'0'*64,'batch':a.batch,'status':'registered_finite_original_prefix'})
 reg={'task':'upper_continue123_actual','batch':a.batch,'master_pid':os.getpid(),'registered_utc':datetime.now(timezone.utc).isoformat(),'reviewed_manifest_sha256':a.reviewed_manifest_sha,'target':str(target),'epsilon':'1/128','worker_count':WORKERS,'selected_original_pending_count':len(tasks),'maximum_new_anchor_starts':selection['max_new_anchor_starts'],'selected_nodes_sha256':selection['selected_nodes_sha256'],'tasks_sha256':selection['tasks_sha256'],'seed_sha256':plan['seed_sha256'],'source_cache_canonical_sha256':oldcachehash,'source_passed_canonical_sha256':oldpassedhash,'parent_checkpoint':parent,'parent_counts':oldcounts,'parent_new_lengths':oldlens,'initial_tree':initialtree,'wall_deadline':None,'no_geometry_children_in_frozen_two_batches':True,'initial_source_schema':'one readonly complete seed with passed/splits/roots/cache_by_epsilon and exactqueue;overlay appends new objects only'}
 write(bdir/'registration.json',reg);write(bdir/'initial_overlay.json',ov);write(bdir/'state.json',ov);journal=open(bdir/'journal.jsonl','ab');saved=json.load(open(W/'saved_fulljet110_proofs.json'));assert sha(W/'saved_fulljet110_proofs.json')==plan['saved_fulljet110_proofs_sha256']
 def emit(kind,body):
  event={'seq':ov['journal_seq']+1,'previous':ov['journal_hash'],'kind':kind,'body':body};event['hash']=digest(event);journal.write((json.dumps(event,sort_keys=True,separators=(',',':'))+'\n').encode());journal.flush();os.fsync(journal.fileno());apply(ov,event,seed['cache_by_epsilon'])
  if ov['journal_seq']%128==0 or kind in ('finish','coordinator_fault'):write(bdir/'state.json',ov)
 context=mp.get_context('spawn');ready=context.Queue();pool=context.Pool(WORKERS,initializer=worker.initialize,initargs=(str(TABLE),TSHA,str(target),ready));info=[];pending={};nextdispatch=0;stop=False
 try:
  for _ in range(WORKERS):info.append(ready.get(timeout=120))
  assert len({r['pid'] for r in info})==WORKERS and all(r['target']==str(target) and r['epsilon']=='1/128' and r['kernel_sha256']==KERNEL_SHA for r in info)
  write(bdir/'pool_registration.json',{'master_pid':os.getpid(),'workers':info,'before_any_H_dispatch':True,'registration_sha256':sha(bdir/'registration.json')});print(json.dumps({'event':'pool_registered','batch':a.batch,'master_pid':os.getpid(),'worker_pids':[r['pid'] for r in info],'selected':len(tasks),'max_new_starts':selection['max_new_anchor_starts'],'registration_sha256':sha(bdir/'registration.json')}),flush=True)
  while pending or nextdispatch<len(tasks) and not stop:
   if (W/'STOP').exists():stop=True
   while not stop and nextdispatch<len(tasks) and len(pending)<WORKERS:
    task=tasks[nextdispatch];assert not Path(task['point_checkpoint']).exists();emit('dispatch',{'task':task})
    if task['cached']:
     got=cached_result(task,seed);emit('received',{'ordinal':nextdispatch,'returned':got})
     if got['kind']!='complete':stop=True;write(W/'STOP',{'ordinal':nextdispatch,'kind':got['kind'],'epsilon':'1/128','target':str(target)})
    else:pending[nextdispatch]=pool.apply_async(worker.evaluate_job,(task,))
    nextdispatch+=1
   progress=False
   for i,job in list(pending.items()):
    if job.ready():
     got=job.get();assert got['target']==str(target) and got['epsilon']=='1/128' and got['pair_cache_key']==tasks[i]['pair_cache_key'];emit('received',{'ordinal':i,'returned':got});del pending[i];progress=True
     if got['kind']!='complete':stop=True;write(W/'STOP',{'ordinal':i,'kind':got['kind'],'epsilon':'1/128','target':str(target)})
   while str(ov['next_merge']) in ov['inflight'] and 'returned' in ov['inflight'][str(ov['next_merge'])]:
    i=ov['next_merge'];ent=ov['inflight'][str(i)];got=ent['returned'];body={'ordinal':i}
    if got['kind']=='complete':body['certificate']=reuse.make_certificate(physical(ent['task']['node'],seed['roots']),got['record'],target,saved[ent['task']['pair_cache_key'][1]]) if ent['task']['cached'] else certificate.make_certificate(physical(ent['task']['node'],seed['roots']),got['record'],target)
    emit('merge',body)
    if ov['next_merge']%128==0:print(json.dumps({'event':'progress','batch':a.batch,'merged':ov['next_merge'],'selected':len(tasks),'cumulative_new_passed':len(ov['new_passed']),'cumulative_new_splits':len(ov['new_splits']),'queue':len(ov['queue']),'seconds':time.monotonic()-began}),flush=True)
   for process in pool._pool:
    if process.pid not in {r['pid'] for r in info} or process.exitcode not in (None,0):raise RuntimeError('worker death/replacement;preserve point/cone fingerprints')
   if not progress:time.sleep(.05)
  pool.close();pool.join()
 except BaseException as exc:
  emit('coordinator_fault',{'exception':repr(exc),'pid':os.getpid(),'epsilon':'1/128','target':str(target),'action':'save prefix,no blind restart'});pool.terminate();pool.join()
 evidence=[{'pid':p.pid,'exitcode':p.exitcode,'alive':p.is_alive()} for p in pool._pool]
 status='H128_g110_target_obstruction' if ov['obstructions_H128'] else 'fault_or_partial' if ov['faults'] else 'finite_original_prefix_complete' if ov['next_merge']==len(tasks) else 'explicit_stop_partial'
 emit('finish',{'status':status,'pool_evidence':evidence,'wall_seconds':time.monotonic()-began});journal.close();t=tree(seed,ov)
 assert digest(seed['cache_by_epsilon'])==oldcachehash and digest(seed['passed'])==oldpassedhash and sha(SEED)==plan['seed_sha256']
 for path,h in package['source_locks'].items():assert sha(path)==h
 points=[]
 for task in tasks[:nextdispatch]:
  path=Path(task['point_checkpoint'])
  if path.exists():
   cp=json.load(open(path));assert cp['point']==task['point'] and cp['epsilon']=='1/128' and cp['target']==str(target) and cp['kernel_sha256']==KERNEL_SHA
   points.append({'ordinal':task['batch_ordinal'],'global_ordinal':task['ordinal'],'node':task['node'],'epsilon':'1/128','target':str(target),'cached':task['cached'],'point':task['point'],'pair_cache_key':task['pair_cache_key'],'point_checkpoint':str(path),'point_checkpoint_sha256':sha(path),'stage':cp['stage'],'worker_pid':cp['pid'],'record':cp.get('record'),'source_raw_record_sha256':cp.get('source_raw_record_sha256')})
 write(bdir/'point_manifest.json',{'model':MODEL,'epsilon':'1/128','target':str(target),'kernel_sha256':KERNEL_SHA,'records':points,'scope':'new production values and explicit oldfulljet references;independent audit pending'})
 summary={'status':status,'batch':a.batch,'master_pid':os.getpid(),'target':str(target),'epsilon':'1/128','selected_original_pending':len(tasks),'new_counts':{k:ov['counts'][k]-oldcounts[k] for k in oldcounts},'cumulative_counts':ov['counts'],'new_passed':len(ov['new_passed'])-oldlens['new_passed'],'new_splits':len(ov['new_splits'])-oldlens['new_splits'],'total_passed':98007+len(ov['new_passed']),'queue':len(ov['queue']),'inflight':len(ov['inflight']),'faults':len(ov['faults']),'obstructions_H128':len(ov['obstructions_H128']),'tree':t,'pool_evidence':evidence,'state_sha256':sha(bdir/'state.json'),'journal_sha256':sha(bdir/'journal.jsonl'),'point_manifest_sha256':sha(bdir/'point_manifest.json'),'wall_seconds':time.monotonic()-began,'no_children_dispatch_in_this_plan':True,'globalg110_qualified':False,'S_credit':0}
 write(bdir/'summary.json',summary);write(bdir/'manifest.json',{'task':'upper_continue123_batch','batch':a.batch,'reviewed_manifest_sha256':a.reviewed_manifest_sha,'files':{p.name:sha(p) for p in sorted(bdir.iterdir()) if p.is_file() and p.name!='manifest.json'}})
 if status=='finite_original_prefix_complete' and not ov['inflight'] and not ov['faults'] and not ov['obstructions_H128']:
  write(current,{'batch':a.batch,'state':str(bdir/'state.json'),'state_sha256':summary['state_sha256'],'journal':str(bdir/'journal.jsonl'),'journal_sha256':summary['journal_sha256'],'selected_original_count_so_far':sum(b['count'] for b in plan['batches'][:a.batch])})
 print(json.dumps(summary),flush=True)
if __name__=='__main__':main()
