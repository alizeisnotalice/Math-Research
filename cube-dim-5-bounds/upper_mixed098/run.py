#!/usr/bin/env python3
"""One fresh finite194 H128 batch. No continuation or blind point restart.

Frozen preparation/kernel/worker files are imported read-only. The journal
contains every transition; replay needs no numerical kernel calls.
"""
from pathlib import Path
from fractions import Fraction as F
from collections import defaultdict
from datetime import datetime, timezone
import hashlib,json,multiprocessing as mp,os,sys,time
import jet_job
import certificate
from cache_model import pair_key,insert

HERE=Path(__file__).resolve().parent
WORK=Path('/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/upper_mixed098')
SOURCE=WORK/'initial_mixed_state.json'
SOURCE_SHA='19294a9bd493524fdc23ec2abe2bb3e70577f04b762b08d4fc57e78af5377480'
TABLE=Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube5-astra-upper/cube-dim-5-bounds/astra_vertex017/boundary_triangulation.json.gz')
TABLE_SHA='8dd2e9895d7465ac599d43feb1b5099506f4e43d1248215eba0a77e11f4d13f9'
KERNEL_SHA='128d98c94e55632336edafe024aa3a27c1486898670a52be89cb39cc22ed6cc0'
PREP_SHA='fe03bfa691b6947a24f13a2f94885962332a77edfce087e2f9e01bc49a3ace15'
TASK_SHA='16cd0a1c303ec8f5f0c0247f7017721bd09e93032905da11f279d3f0a0914f6e'
THEORY=Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube5-astra-upper/cube-dim-5-bounds/audit_mixed099/manifest.json')
THEORY_SHA='5df7eaf9511a3b91a5a878edd7d6fb5f45ef1c69bdda4ae6886f21145ba49b3c'
PRIMARY=Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube-dim-5/cube-dim-5-bounds/round013/mixed098_primary_review.json')
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
  j=max(range(4),key=lambda j:hi[j]-lo[j]);m=(lo[j]+hi[j])/2
  assert bit in '01'
  if bit=='0':hi[j]=m
  else:lo[j]=m
 return lo,hi
def partition(source,ov):
 byroot=defaultdict(list);leaves=set();internals=set()
 nodes=[x['node'] for x in source['passed']+ov['new_passed']]+ov['queue']+[e['task']['node'] for e in ov['inflight'].values()]
 for r,p in nodes:
  assert (r,p) not in leaves;leaves.add((r,p));byroot[r].append(p)
 assert len(byroot)==len(source['roots'])==52360
 for r in range(52360):
  ps=sorted(byroot[r]);assert ps and all(not b.startswith(a) for a,b in zip(ps,ps[1:]))
  d=max(map(len,ps));assert sum(1<<(d-len(p)) for p in ps)==1<<d
 for row in source['splits']+ov['new_splits']:
  r,p=row['node'];assert (r,p) not in internals;internals.add((r,p))
  assert row['axis']==len(p)%4 and row['children']==[[r,p+'0'],[r,p+'1']]
 assert not leaves&internals
 for r,p in leaves|internals:
  if p:assert (r,p[:-1]) in internals
 assert len(leaves)==52360+len(internals)
 return {'roots_exact_Kraft_one':52360,'leaves':len(leaves),'splits':len(internals),'antichain':True,'closed_bisection':True}
def apply(ov,event):
 b=event['body'];kind=event['kind'];c=ov['counts']
 if kind=='dispatch':
  task=b['task'];k=str(task['ordinal']);assert ov['queue'].pop(0)==task['node'] and k not in ov['inflight']
  ov['inflight'][k]={'task':task};c['new_anchor_starts']+=1
 elif kind=='received':
  ent=ov['inflight'][str(b['ordinal'])];assert 'returned' not in ent;ent['returned']=b['returned']
  c['kernel_segment_calls']+=b['returned'].get('calls',0)
 elif kind=='merge':
  ordinal=b['ordinal'];assert ordinal==ov['next_merge'];ov['next_merge']+=1
  ent=ov['inflight'][str(ordinal)];got=ent['returned'];task=ent['task'];node=task['node']
  if got['kind']=='complete':
   rec=got['record'];assert rec['epsilon']=='1/128' and rec['point']==task['point']
   p=pair_key(F(1,128),rec['point'])[1];assert p not in ov['new_cache_H128'];ov['new_cache_H128'][p]=rec;c['new_complete']+=1
   cert=b['certificate'];assert cert['leaf_epsilon']=='1/128' and cert['penalty_coefficient']=='1024'
   if cert['passed']:
    ov['new_passed'].append({'node':node,'certificate':cert,'leaf_epsilon':'1/128','lower_bound_object':'G','component_lower_function':'H_epsilon','target':'21/200000'})
   else:
    children=[[node[0],node[1]+'0'],[node[0],node[1]+'1']]
    ov['new_splits'].append({'node':node,'axis':len(node[1])%4,'children':children,'failed_certificate':cert,'leaf_epsilon':'1/128'})
    ov['queue'].extend(children)
   del ov['inflight'][str(ordinal)]
  elif got['kind']=='H128_target_obstruction':
   ov['obstructions_H128'].append({'node':node,'returned':got,'scope':'fixed H128 target only,not full G/C5 refutation'})
  else:ov['faults'].append({'node':node,'returned':got})
 elif kind=='finish':ov.update(b)
 else:raise ValueError(kind)
 ov['journal_seq']=event['seq'];ov['journal_hash']=event['hash']
def main():
 began=time.monotonic()
 assert '--execute-authorized-first194' in sys.argv,'Explicit main authorization flag required'
 assert sha(SOURCE)==SOURCE_SHA and sha(HERE/'manifest.json')==PREP_SHA
 assert sha(HERE/'epsilon_kernel.py')==KERNEL_SHA and sha(TABLE)==TABLE_SHA
 assert sha(THEORY)==THEORY_SHA and PRIMARY.is_file()
 prep=json.load(open(HERE/'manifest.json'))
 # The preparation manifest binds its exact source files; never update it.
 for name,expected in prep.get('files',{}).items():
  if isinstance(expected,str) and (HERE/name).is_file():assert sha(HERE/name)==expected
 assert sha(WORK/'execution_tasks_194.json')==TASK_SHA
 overlap=json.load(open(WORK/'overlap095_verification.json'));assert overlap['intersection_count']==0 and overlap['status']=='exact_registry_intersection_empty'
 assert sha(overlap['095_registration_path'])==overlap['095_registration_sha256']
 tasks=json.load(open(WORK/'execution_tasks_194.json'));assert len(tasks)==194
 assert len({tuple(t['pair_cache_key']) for t in tasks})==194
 assert not (WORK/'STOP').exists(),'Read STOP cause; do not blindly restart'
 source=json.load(open(SOURCE));assert len(source['passed'])==100804 and len(source['queue'])==194 and len(source['splits'])==48638
 assert len(source['cache_by_epsilon']['1/64'])==184796 and not source['cache_by_epsilon']['1/128'] and not source['inflight'] and not source['resume_jobs']
 assert [t['node'] for t in tasks]==source['queue'] and all(t['epsilon']=='1/128' for t in tasks)
 oldcache=digest(source['cache_by_epsilon']['1/64']);oldpassed=digest(source['passed']);oldother={k:digest(v) for k,v in source.items() if k not in ('cache_by_epsilon','passed')}
 bdir=WORK/'batch0001'
 # Fresh-only ownership gate: never overwrite a prior cone/point checkpoint.
 assert not bdir.exists(),'Existing batch: preserve it and inspect epsilon/kernel/point/precision fingerprints before any separately authorized recovery'
 assert all(not Path(t['point_checkpoint']).exists() for t in tasks)
 bdir.mkdir(exist_ok=False)
 ov={'source':str(SOURCE),'source_sha256':SOURCE_SHA,'new_cache_H128':{},'new_passed':[],'new_splits':[],'queue':list(source['queue']),
  'inflight':{},'obstructions_H128':[],'faults':[],'counts':{'new_anchor_starts':0,'new_complete':0,'kernel_segment_calls':0,'H128_cache_hits':0},
  'status':'registered_finite194','next_merge':0,'journal_seq':0,'journal_hash':'0'*64,'epsilon':'1/128','target':'21/200000','lower_bound_object':'G'}
 initialtree=partition(source,ov)
 code={p.name:sha(p) for p in HERE.glob('*.py')}
 reg={'task':'upper_mixed098_actual_first194','master_pid':os.getpid(),'registered_utc':datetime.now(timezone.utc).isoformat(),
  'source':str(SOURCE),'source_sha256':SOURCE_SHA,'preparation_manifest_sha256':PREP_SHA,'code_sha256':code,
  'epsilon':'1/128','target':'21/200000','lower_bound_object':'G','Taylor_penalty':'1024','closed_corners':16,
  'source_H64_cache_canonical_sha256':oldcache,'source_H64_passed_canonical_sha256':oldpassed,'source_other_canonical_sha256':oldother,
  'workers':12,'finite_selected_nodes':source['queue'],'finite_selected_sha256':digest(source['queue']),'maximum_new_anchor_starts':194,
  'new_H128_cache_hits':0,'H64_values_not_reused_as_H128':True,'execution_tasks_sha256':TASK_SHA,'wall_deadline':None,'alarm':False,
  'kernel_sha256':KERNEL_SHA,'table':str(TABLE),'table_sha256':TABLE_SHA,'theory099_manifest':str(THEORY),'theory099_manifest_sha256':THEORY_SHA,
  'primary_review':str(PRIMARY),'primary_review_sha256':sha(PRIMARY),'overlap095_verification_sha256':sha(WORK/'overlap095_verification.json'),
  'stops':['finite194 exhausted','strict H128 upper<target:retain leaf','actual numerical/process/state fault:save prefix','explicit STOP'],
  'restart_policy':'fresh-only; existing point checkpoint rejected; inspect epsilon/kernel/point/precision before separately authorized recovery',
  'point_checkpoint_schema':'frozen jet_job explicit epsilon,pair_cache_key,point,kernel_sha256,precision_index,partial cone prefix,attempts',
  'overlay_schema':'source immutable+new_cache_H128,new_passed,new_splits,queue,inflight,obstructions_H128,faults,counts',
  'journal_schema':'canonical JSON hash chain dispatch/received/merge/finish; initial_overlay replay without H',
  'initial_tree':initialtree,'production_numeric_scope':'standalone DD/Arb production,not independent numeric audit','python':sys.executable}
 write(bdir/'registration.json',reg);write(bdir/'initial_overlay.json',ov);write(bdir/'state.json',ov)
 journal=open(bdir/'journal.jsonl','ab')
 def emit(kind,body):
  e={'seq':ov['journal_seq']+1,'previous':ov['journal_hash'],'kind':kind,'body':body};e['hash']=digest(e)
  journal.write((json.dumps(e,sort_keys=True,separators=(',',':'))+'\n').encode());journal.flush();os.fsync(journal.fileno())
  apply(ov,e);write(bdir/'state.json',ov)
 context=mp.get_context('spawn');ready=context.Queue();pool=context.Pool(12,initializer=jet_job.initialize,initargs=(str(TABLE),TABLE_SHA,KERNEL_SHA,ready))
 workerinfo=[];pending={};nextdispatch=0;stop=False
 try:
  for _ in range(12):workerinfo.append(ready.get(timeout=120))
  assert len({x['pid'] for x in workerinfo})==12 and all(x['epsilon']=='1/128' and x['kernel_sha256']==KERNEL_SHA for x in workerinfo)
  write(bdir/'pool_registration.json',{'master_pid':os.getpid(),'workers':workerinfo,'registered_before_any_H128_dispatch':True,'run_registration_sha256':sha(bdir/'registration.json')})
  print(json.dumps({'event':'pool_registered','master_pid':os.getpid(),'worker_pids':[x['pid'] for x in workerinfo],'source_sha256':SOURCE_SHA,'registration_sha256':sha(bdir/'registration.json')}),flush=True)
  while pending or (nextdispatch<194 and not stop):
   if (WORK/'STOP').exists():stop=True
   while not stop and nextdispatch<194 and len(pending)<12:
    task=tasks[nextdispatch];assert not Path(task['point_checkpoint']).exists()
    emit('dispatch',{'task':task});pending[nextdispatch]=pool.apply_async(jet_job.evaluate_job,(task,));nextdispatch+=1
   progress=False
   for ordinal,job in list(pending.items()):
    if job.ready():
     got=job.get();task=tasks[ordinal];assert got['epsilon']=='1/128' and got['pair_cache_key']==task['pair_cache_key']
     emit('received',{'ordinal':ordinal,'returned':got});del pending[ordinal];progress=True
     if got['kind']!='complete':stop=True;write(WORK/'STOP',{'cause':'noncomplete H128 result','ordinal':ordinal,'kind':got['kind'],'source_sha256':SOURCE_SHA})
   while str(ov['next_merge']) in ov['inflight'] and 'returned' in ov['inflight'][str(ov['next_merge'])]:
    ordinal=ov['next_merge'];ent=ov['inflight'][str(ordinal)];got=ent['returned'];body={'ordinal':ordinal}
    if got['kind']=='complete':body['certificate']=certificate.make_certificate(physical(ent['task']['node'],source['roots']),got['record'])
    emit('merge',body)
    if ov['next_merge']%12==0:print(json.dumps({'event':'progress','merged':ov['next_merge'],'new_passed':len(ov['new_passed']),'new_splits':len(ov['new_splits']),'queued':len(ov['queue']),'seconds':time.monotonic()-began}),flush=True)
   for process in pool._pool:
    if process.pid not in {x['pid'] for x in workerinfo} or process.exitcode not in (None,0):raise RuntimeError('unexpected worker death/replacement; preserve checkpoints')
   if not progress:time.sleep(.1)
  pool.close();pool.join()
 except BaseException as exc:
  write(bdir/'coordinator_fault.json',{'exception':repr(exc),'pid':os.getpid(),'epsilon':'1/128','kernel_sha256':KERNEL_SHA,'action':'preserve all checkpoints; no blind restart'})
  pool.terminate();pool.join();ov['faults'].append({'coordinator_exception':repr(exc)})
  write(bdir/'state.json',ov)
 finally:
  evidence=[{'pid':p.pid,'exitcode':p.exitcode,'alive':p.is_alive()} for p in pool._pool]
 status='H128_target_obstruction' if ov['obstructions_H128'] else 'fault_or_partial' if ov['faults'] else 'finite194_complete' if ov['next_merge']==194 else 'explicit_stop_partial'
 emit('finish',{'status':status,'pool_evidence':evidence,'wall_seconds':time.monotonic()-began});journal.close()
 tree=partition(source,ov)
 assert digest(source['cache_by_epsilon']['1/64'])==oldcache and digest(source['passed'])==oldpassed
 assert {k:digest(v) for k,v in source.items() if k not in ('cache_by_epsilon','passed')}==oldother and not source['cache_by_epsilon']['1/128']
 assert sha(SOURCE)==SOURCE_SHA and sha(HERE/'manifest.json')==PREP_SHA
 records=[]
 for ordinal in range(nextdispatch):
  task=tasks[ordinal];p=Path(task['point_checkpoint'])
  if p.exists():
   cp=json.load(open(p));assert cp['epsilon']=='1/128' and cp['kernel_sha256']==KERNEL_SHA and cp['point']==task['point'] and cp['pair_cache_key']==task['pair_cache_key']
   records.append({'ordinal':ordinal,'node':task['node'],'point':task['point'],'epsilon':'1/128','point_checkpoint':str(p),'point_checkpoint_sha256':sha(p),'stage':cp['stage'],'worker_pid':cp['pid'],'record':cp.get('record')})
 write(bdir/'point_manifest.json',{'source_sha256':SOURCE_SHA,'kernel_sha256':KERNEL_SHA,'execution_tasks_sha256':TASK_SHA,'records':records,'scope':'fresh H128 production values,independent audit pending'})
 validation={'status':'author_structural_checks_passed','source_old_cache_direct_unchanged':True,'source_old_passed_direct_unchanged':True,'source_other_fields_unchanged':True,'tree':tree,'new_corner_certificates':len(ov['new_passed'])+len(ov['new_splits']),'new_corners':16*(len(ov['new_passed'])+len(ov['new_splits']))}
 write(bdir/'author_validation.json',validation)
 fullpath=None
 if status=='finite194_complete' and not ov['inflight'] and not ov['faults'] and not ov['obstructions_H128']:
  final={**source,'task':'upper_mixed098','status':'mixed_G105_closed_candidate' if not ov['queue'] else 'mixed_G105_partial','passed':source['passed']+ov['new_passed'],
   'splits':source['splits']+ov['new_splits'],'queue':ov['queue'],'inflight':{},'resume_jobs':{},'cache_by_epsilon':{'1/64':source['cache_by_epsilon']['1/64'],'1/128':dict(ov['new_cache_H128'])},
   'accepted_leaf_epsilon':'per leaf,1/64 or1/128','g105_closed_cover':not ov['queue'],'new_H128_records':len(ov['new_cache_H128']),'production_authorized':True,
   'final098_source_sha256':SOURCE_SHA,'final098_overlay_sha256':sha(bdir/'state.json'),'final098_registration_sha256':sha(bdir/'registration.json'),
   'final098_scope':'same G mixed-epsilon candidate;await main and100;no I or C5 conclusion'}
  final['fingerprint098_final']=digest({'source':SOURCE_SHA,'overlay':sha(bdir/'state.json'),'journal':sha(bdir/'journal.jsonl'),'point_manifest':sha(bdir/'point_manifest.json')})
  fullpath=WORK/'final_state.json';write(fullpath,final)
 summary={'status':status,'master_pid':os.getpid(),'pool_evidence':evidence,'source_sha256':SOURCE_SHA,'selected':194,'counts':ov['counts'],'H64_passed_retained':100804,
  'new_H128_passed':len(ov['new_passed']),'new_H128_splits':len(ov['new_splits']),'queue':len(ov['queue']),'inflight':len(ov['inflight']),'H128_obstructions':len(ov['obstructions_H128']),
  'faults':len(ov['faults']),'tree':tree,'state_path':str(bdir/'state.json'),'state_sha256':sha(bdir/'state.json'),'point_manifest_sha256':sha(bdir/'point_manifest.json'),
  'full_state_path':str(fullpath) if fullpath else None,'full_state_sha256':sha(fullpath) if fullpath else None,'wall_seconds':time.monotonic()-began,'independent_numeric_audit_pending':True}
 write(bdir/'summary.json',summary);print(json.dumps(summary),flush=True)
 manifest={'task':'upper_mixed098_batch0001','source_sha256':SOURCE_SHA,'preparation_manifest_sha256':PREP_SHA,'code_sha256':code,
  'files':{p.name:sha(p) for p in sorted(bdir.iterdir()) if p.is_file() and p.name!='manifest.json'},'final_state_sha256':summary['full_state_sha256']}
 write(bdir/'manifest.json',manifest)
if __name__=='__main__':main()
