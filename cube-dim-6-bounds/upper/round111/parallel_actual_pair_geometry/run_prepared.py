"""111 parent: one release,12 disjoint workers, finite checkpoint continuation.

Latest user policy has no time/CPU cap. No geometry or cell executes before
the exact parent release. Only the parent writes snapshots and checkpoints.
"""
import sys
sys.dont_write_bytecode=True
import datetime,gzip,hashlib,json,os,resource,selectors,subprocess,time,traceback
from fractions import Fraction as Q
from pathlib import Path
from protocol import read_frame,send_frame
P=Path(__file__).resolve().parent
THREAD_ENV=['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS']
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def write_json(path,data):
 path=P/path;path.parent.mkdir(parents=True,exist_ok=True);tmp=path.with_name(path.name+'.tmp')
 tmp.write_text(json.dumps(data,indent=2,default=str)+'\n')
 if tmp.stat().st_size>=95*1024*1024:raise RuntimeError('95MiB per JSON archive file guard')
 tmp.replace(path)
def write_gzip(path,data):
 path=P/path;path.parent.mkdir(parents=True,exist_ok=True);tmp=path.with_name(path.name+'.tmp')
 with tmp.open('wb')as f:
  with gzip.GzipFile(fileobj=f,mode='wb',compresslevel=1,mtime=0)as z:
   import io
   with io.TextIOWrapper(z,encoding='utf-8')as text:json.dump(data,text,separators=(',',':'),default=str)
 if tmp.stat().st_size>=95*1024*1024:raise RuntimeError('95MiB per archive chunk guard')
 tmp.replace(path)
 return {'path':str(path.relative_to(P)),'bytes':path.stat().st_size,'sha256':sha(path)}
def snapshot(job_id,tag,state,transport_sha,stats):
 folder=f'checkpoints/{job_id}/{tag}';leaves=state['leaves'];parts=[];ids=sorted(leaves)
 for pos in range(0,len(ids),4096):
  block=[leaves[i]for i in ids[pos:pos+4096]]
  parts.append(write_gzip(f'{folder}/leaves-{pos//4096:03d}.json.gz',{'leaves':block}))
 header={k:v for k,v in state.items()if k!='leaves'}
 header.update({'transport_sha256':transport_sha,'producer_stats':stats,'leaf_chunks':parts,'leaf_count':len(ids),'complete_cover_via_initial_seeds_and_refinements':True,'geometry_sha256':sha(P/'geometry_certificate.json')})
 head=write_gzip(f'{folder}/header.json.gz',header)
 index={'job_id':job_id,'tag':tag,'header':head,'leaf_chunks':parts,'shape_lower':str(state['lower']),'shape_upper':str(state['upper']),'leaves':len(ids),'cell_bounds':state['cell_bounds'],'transport_sha256':transport_sha}
 write_json(f'{folder}/index.json',index)
 return {'index_path':f'{folder}/index.json','index_sha256':sha(P/f'{folder}/index.json'),'lower':state['lower'],'upper':state['upper'],'leaves':len(ids),'cell_bounds':state['cell_bounds'],'transport_sha256':transport_sha,'tag':tag}
def cleanup(children,normal):
 actions=[]
 if normal:
  for jid,p in children.items():
   if p.poll()is None:
    try:send_frame(p.stdin,{'command':'STOP'})
    except Exception:pass
 else:
  for jid,p in children.items():
   if p.poll()is None:p.terminate();actions.append({'job':jid,'action':'terminate_owned_child_after_failure_or_interrupt'})
 deadline=time.monotonic()+20
 for jid,p in children.items():
  if p.poll()is None:
   try:p.wait(timeout=max(0.01,deadline-time.monotonic()))
   except subprocess.TimeoutExpired:pass
 for jid,p in children.items():
  if p.poll()is None:p.kill();actions.append({'job':jid,'action':'kill_owned_child_after_bounded_cleanup_wait'})
 for jid,p in children.items():
  try:p.wait(timeout=10)
  except subprocess.TimeoutExpired:raise RuntimeError(f'owned child{jid}not reaped')
 return actions

def main():
 release_path=P/'parent_release.json'
 if not release_path.exists():raise SystemExit('111 PREPARED_NOT_RELEASED')
 rel=json.loads(release_path.read_text());reg=json.loads((P/'registration.json').read_text())
 assert rel['authorized']is True and rel['round']=='111'and rel['registration_sha256']==sha(P/'registration.json')
 for name,h in reg['code_sha256'].items():assert sha(P/name)==h
 assert rel['script_sha256']==sha(Path(__file__))
 assert sha(P/'source_manifest.json')==reg['source_manifest_sha256']and sha(P/'upstream_constants.json')==reg['upstream_constants_sha256']
 assert sha(P/'qualification_addendum.json')==reg['qualification_addendum_sha256']
 for source in json.loads((P/'source_manifest.json').read_text())['sources']:
  assert sha(Path(source['path']))==source['sha256']
 for name in ['started.json','checkpoint_index.json','results.json','run_receipt.json','geometry_certificate.json']:
  if(P/name).exists():raise SystemExit('111 one finite attempt only; no automatic restart')
 assert len(reg['jobs'])==12 and len({j['id']for j in reg['jobs']})==12 and reg['max_workers']==12
 for k in(2,3,4,5):
  spans=sorted([tuple(map(Q,j['u']))for j in reg['jobs']if j['k']==k])
  assert spans==[(Q(1,4),Q(5,12)),(Q(5,12),Q(7,12)),(Q(7,12),Q(3,4))]
 for e in THREAD_ENV:os.environ[e]='1'
 os.environ['PYTHONDONTWRITEBYTECODE']='1'
 start=time.monotonic();start_utc=datetime.datetime.now(datetime.timezone.utc).isoformat();before=resource.getrusage(resource.RUSAGE_CHILDREN)
 children={};selector=selectors.DefaultSelector();stderr={};checkpoints={};progress={};phase_records=[];phase1_lineage={};status='FIRST_FAILURE_STOP';failure=None;normal=False
 write_json('started.json',{'round':'111','pid':os.getpid(),'started_utc':start_utc,'release_sha256':sha(release_path),'runtime_limit':None,'CPU_limit':None,'mathematical_workers':12,'worker_BLAS_threads':1,'k1_read_only':True})
 try:
  from math_core import Arithmetic,generate_geometry
  ar=Arithmetic();A=ar.A
  geometry={'status':'COMPLETE_NEW_K2_TO5_GEOMETRY','by_k':{str(k):generate_geometry(k)for k in(2,3,4,5)},'new_pair_visits':7440,'k1_generated':False}
  write_json('geometry_certificate.json',geometry)
  up=json.loads((P/'upstream_constants.json').read_text());k1L=Q(up['k1_shape']['lower']);k1U=Q(up['k1_shape']['upper']);paid=Q(up['paid_CEH_I_upper'])
  pi=ar.arb.pi();q=Q(6,5);Kold=3*A(3).sqrt()/(64*pi**4);Kref=1/(8*pi**5)
  scale_old=Kold**A(q);scale_ref=Kref**A(q);scales={'old098_valid':(ar.end(scale_old,'lower'),ar.end(scale_old,'upper')),'ref108_frozen_reference':(ar.end(scale_ref,'lower'),ar.end(scale_ref,'upper'))}
  target=64*A(Q(1,65000)/6)**A(q);Tlo=ar.end(target,'lower');Thi=ar.end(target,'upper')
  def aggregate(stage):
   L=k1L;U=k1U;byk={}
   for k,n in[(2,5),(3,10),(4,10),(5,5)]:
    members=[j['id']for j in reg['jobs']if j['k']==k];assert len(members)==3 and all(j in checkpoints for j in members)
    kl=sum((checkpoints[j]['lower']for j in members),Q(0));ku=sum((checkpoints[j]['upper']for j in members),Q(0));L+=n*kl;U+=n*ku;byk[str(k)]={'multiplicity':n,'shape_lower':str(kl),'shape_upper':str(ku)}
   comparisons={}
   for name,(sl,su)in scales.items():
    il=sl*L;iu=su*U;Slo=Tlo-paid-iu;Shi=Thi-paid-il
    sep='STRICTLY_FITS_CONDITIONAL_PLAN'if Slo>0 else'STRICT_FIXED_EXPRESSION_OBSTRUCTION'if Shi<0 else'UNSEPARATED'
    comparisons[name]={'middle_I_lower_expression':str(il),'middle_I_upper_certificate':str(iu),'CEH_plus_middle_I_lower_expression':str(paid+il),'CEH_plus_middle_I_upper_certificate':str(paid+iu),'S_plan_lower':str(Slo),'S_plan_upper':str(Shi),'separation':sep,'g0':'unproved globally','far':'UNPAID, zero only in this plan diagnosis','reference_assumption':reg['reference_status']if name.startswith('ref')else'old098 source coefficient'}
   return {'phase':stage,'total_middle_shape_lower':str(L),'total_middle_shape_upper':str(U),'k1_shape_read_only':{'lower':str(k1L),'upper':str(k1U)},'new_k2_to5':byk,'comparisons':comparisons,'checkpoint_indices':{j:{k:v for k,v in d.items()if k not in('lower','upper')}for j,d in checkpoints.items()}}
  for job in reg['jobs']:
   proc=subprocess.Popen([sys.executable,str(P/'worker_entry.py')],cwd=P,stdin=subprocess.PIPE,stdout=subprocess.PIPE,stderr=subprocess.PIPE,start_new_session=True,env=os.environ.copy())
   children[job['id']]=proc;stderr[job['id']]=bytearray();selector.register(proc.stdout,selectors.EVENT_READ,(job['id'],'stdout'));selector.register(proc.stderr,selectors.EVENT_READ,(job['id'],'stderr'))
   send_frame(proc.stdin,{'command':'BOOT','job_id':job['id'],'release_sha256':sha(release_path),'geometry_sha256':sha(P/'geometry_certificate.json')})
  last_print=0
  for phase in(1,2):
   done=set()
   while len(done)<12:
    events=selector.select(timeout=1)
    for key,_ in events:
     jid,kind=key.data;proc=children[jid]
     if kind=='stderr':
      data=os.read(proc.stderr.fileno(),65536)
      if data:stderr[jid].extend(data[:max(0,262144-len(stderr[jid]))])
      else:selector.unregister(proc.stderr)
      continue
     message,transport_sha=read_frame(proc.stdout)
     assert message['job_id']==jid
     if message['kind']=='failure':
      if message.get('state')is not None:checkpoints[jid]=snapshot(jid,f'failure-phase{phase}',message['state'],transport_sha,{'failure':message['failure']})
      raise RuntimeError(message['failure'])
     assert message['phase']==phase
     progress[jid]={k:v for k,v in message.items()if k!='state'}
     if message['kind']in('partial_state','phase_done'):
      st=message['state']
      if phase==2:
       prior=phase1_lineage[jid];assert st['phase_entry_counts']['cell_bounds']>=prior['cell_bounds']
       prefix=hashlib.sha256(json.dumps(st['refinements'][:prior['refinements']],sort_keys=True,separators=(',',':')).encode()).hexdigest();assert prefix==prior['refinement_prefix_sha256']
      tag=f'phase{phase}-final'if message['kind']=='phase_done'else f'phase{phase}-partial-{st["cell_bounds"]}'
      checkpoints[jid]=snapshot(jid,tag,st,transport_sha,{k:v for k,v in message.items()if k!='state'})
      if message['kind']=='phase_done':
       assert jid not in done and message['pid']==proc.pid;done.add(jid)
       if phase==1:phase1_lineage[jid]={'pid':proc.pid,'transport_sha256':transport_sha,'cell_bounds':st['cell_bounds'],'refinements':len(st['refinements']),'refinement_prefix_sha256':hashlib.sha256(json.dumps(st['refinements'],sort_keys=True,separators=(',',':')).encode()).hexdigest()}
     write_json('checkpoint_index.json',{'phase':phase,'phase_done_jobs':sorted(done),'committed_snapshots':checkpoints,'progress':progress,'elapsed_seconds':time.monotonic()-start})
    for jid,proc in children.items():
     if proc.poll()is not None and jid not in done:raise RuntimeError(f'worker{jid}exited before phase completion: {proc.returncode};{bytes(stderr[jid]).decode(errors="replace")}')
    now=time.monotonic()
    if now-last_print>=5:
     print(json.dumps({'stage':'parallel_phase','phase':phase,'done_jobs':len(done),'leaves_seen':sum(int(x.get('leaves',0))for x in progress.values()),'elapsed_seconds':now-start}),flush=True);last_print=now
   if phase==1:write_json('phase1_lineage.json',phase1_lineage)
   summary=aggregate(phase);phase_records.append(summary);write_json(f'phase{phase}_summary.json',summary)
   refsep=summary['comparisons']['ref108_frozen_reference']['separation']
   if refsep!='UNSEPARATED' or phase==2:
    final=summary;normal=True;break
   for jid,proc in children.items():send_frame(proc.stdin,{'command':'CONTINUE','next_phase':2,'previous_snapshot_transport_sha256':phase1_lineage[jid]['transport_sha256']})
  actions=cleanup(children,True)
  assert all(proc.returncode==0 for proc in children.values())
  result={'status':'COMPLETE_FINITE_K2_TO5_PARALLEL_BATCH','round':'111','phases_completed':len(phase_records),'phase_records':phase_records,'final':final,'geometry_sha256':sha(P/'geometry_certificate.json'),'phase1_lineage_sha256':sha(P/'phase1_lineage.json'),'new_pair_visits':7440,'k1_recomputed':False,'mathematical_workers':12,'BLAS_threads_each':1,'worker_exit_codes':{j:p.returncode for j,p in children.items()},'cleanup_actions':actions,'scales':{'Kold':ar.enc(Kold),'Kref':ar.enc(Kref),'old_q':ar.enc(scale_old),'ref_q':ar.enc(scale_ref)},'reference_status':reg['reference_status'],'global_g_I_S_status':'UNPAID: globalg0/far; reference qualification recorded separately in reference_status','new_C6_bound':False,'code_sha256':reg['code_sha256'],'registration_sha256':sha(P/'registration.json'),'upstream_constants_sha256':sha(P/'upstream_constants.json')}
  write_json('results.json',result);status='COMPLETED'
 except KeyboardInterrupt:failure=traceback.format_exc();status='INTERRUPTED_STOP_NO_RETRY'
 except Exception:failure=traceback.format_exc()
 finally:
  cleanup_actions=[]
  if any(p.poll()is None for p in children.values()):cleanup_actions=cleanup(children,False)
  selector.close();after=resource.getrusage(resource.RUSAGE_CHILDREN)
  receipt={'status':status,'round':'111','pid':os.getpid(),'started_utc':start_utc,'receipt_constructed_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'wall_before_final_receipt_write':time.monotonic()-start,'parent_CPU_total_before_receipt_write':time.process_time(),'reaped_children_CPU_seconds':after.ru_utime+after.ru_stime-before.ru_utime-before.ru_stime,'owned_worker_pids':{j:p.pid for j,p in children.items()},'owned_worker_exit_codes':{j:p.returncode for j,p in children.items()},'all_owned_workers_reaped':all(p.poll()is not None for p in children.values()),'cleanup_actions_on_failure':cleanup_actions,'runtime_limit':None,'CPU_limit':None,'phases_completed':len(phase_records),'no_phase1_or_k1_recalculation':True,'first_failure':failure,'retries':0,'worker_stderr':{j:bytes(b).decode(errors='replace')for j,b in stderr.items()},'release_sha256':sha(release_path)}
  write_json('run_receipt.json',receipt);print(json.dumps({'status':status,'phases_completed':len(phase_records),'owned_workers':len(children),'all_reaped':receipt['all_owned_workers_reaped'],'wall_before_receipt':receipt['wall_before_final_receipt_write']}),flush=True)
 if status!='COMPLETED':raise SystemExit(1)
if __name__=='__main__':main()
