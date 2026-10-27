"""Main-release gated scheduler/observer, at most15 local single-core audit workers; no native oracle."""
import sys,os,json,hashlib,subprocess,time,resource,traceback
from pathlib import Path
P=Path(__file__).resolve().parent
rd=lambda p:json.loads(Path(p).read_text());sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
def save(n,x):(P/n).write_text(json.dumps(x,indent=2)+'\n')
def requirements():
 return{'authorized':True,'scope':'C129_TERMINAL_NEW_CERTIFICATES_ONLY','worker_script_sha256':sha(P/'check_new_certificates_worker.py'),'scheduler_sha256':sha(P/'run_registered_audit.py'),'registration_sha256':sha(P/'registration.json'),'prepared_jobs_sha256':sha(P/'prepared_new_jobs.json'),'structure_review_sha256':sha(P/'structure_review.json'),'new_certificates':18257,'new_native_records_to_validate':18257,'frozen_old_jets_only':4,'new_Arb768_corner_exp_max':584224,'workers_max':15,'old_baseline_Arb_or_validator_replay':0,'time_limit':None}
def pins():
 for r in rd(P/'terminal_source_pins.json')['files']:assert sha(r['path'])==r['sha256']and Path(r['path']).stat().st_size==r['bytes']
 for r in rd(P/'prepared_new_jobs.json')['shards']:assert sha(r['path'])==r['sha256']
def main():
 rel=rd(P/'parent_release.json')
 for k,v in requirements().items():assert rel.get(k)==v,'main release mismatch:'+k
 n=int(rel['workers']);assert 1<=n<=15
 assert rd(P/'structure_review.json')['status']=='PASS_TERMINAL_STRUCTURE_AND_SAVED_ARITHMETIC_PENDING_NEW_ARB768_ENDPOINTS'
 assert not(P/'STOP_REQUESTED.json').exists()and not(P/'actual_tool_receipt.json').exists()
 os.nice(max(0,10-os.getpriority(os.PRIO_PROCESS,0)))
 with(P/'started_audit.json').open('x')as f:json.dump({'pid':os.getpid(),'unix_time':time.time(),'release_sha256':sha(P/'parent_release.json'),'workers':n},f)
 pins();out=P/'worker_logs';out.mkdir(exist_ok=False);pending=list(range(15));active={};completed=[];maxactive=0;t=time.monotonic();before=resource.getrusage(resource.RUSAGE_CHILDREN);failure=None
 try:
  while pending or active:
   while pending and len(active)<n and not(P/'STOP_REQUESTED.json').exists():
    shard=pending.pop(0);log=(out/f'{shard:02d}.stdout').open('xb');err=(out/f'{shard:02d}.stderr').open('xb');child=subprocess.Popen(['/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/venv/bin/python','-B',str(P/'check_new_certificates_worker.py'),str(shard)],cwd=P,stdout=log,stderr=err,stdin=subprocess.DEVNULL);active[shard]=(child,log,err);maxactive=max(maxactive,len(active))
   for shard,(child,log,err)in list(active.items()):
    if child.poll()is None:continue
    code=child.wait();log.close();err.close();completed.append({'shard':shard,'pid':child.pid,'exit_code':code});del active[shard]
    if code and not(P/'STOP_REQUESTED.json').exists():save('STOP_REQUESTED.json',{'status':'WORKER_EXIT_FAILURE','shard':shard,'exit_code':code})
   save('scheduler_progress.json',{'pending':pending,'active':{k:v[0].pid for k,v in active.items()},'completed':completed,'max_active_workers':maxactive})
   if(P/'STOP_REQUESTED.json').exists()and not active:break
   if active:time.sleep(.1)
  assert not pending and len(completed)==15 and all(r['exit_code']==0 for r in completed),'finite audit stopped; retain checkpoints, no automatic retry'
  totals={k:0 for k in ['jobs','new_native_records','frozen_old_jets','certificates','corners']}
  for shard in range(15):
   result=rd(P/'worker_results'/f'result_{shard:02d}.json');assert result['status']=='PASS'
   for k in totals:totals[k]+=result[k]
  assert totals=={'jobs':18261,'new_native_records':18257,'frozen_old_jets':4,'certificates':18257,'corners':584224}
  pins();structure=rd(P/'structure_review.json');save('review_results.json',{'status':'PASS_COMPLETE_C129_NEW_SAVED_MATH_AND_TERMINAL_COVERAGE_BRIDGE',**totals,'continued_roots':55,'closed_continued_roots':20,'previously_closed_preserved':6,'newly_closed':14,'budget_stopped_roots':35,'unknown15_preserved_no_reissue':True,'old15193_validator_Arb_replay':0,'new_native_oracle_calls':0,'new_precision_bits':768,'global_g_or_new_C6_proved':False,'structure_review_sha256':sha(P/'structure_review.json')})
 except BaseException:failure=traceback.format_exc()
 finally:
  for shard,(child,log,err)in active.items():
   child.wait();log.close();err.close()
  after=resource.getrusage(resource.RUSAGE_CHILDREN);save('actual_tool_receipt.json',{'status':'PASS'if failure is None else'FAILED_NO_RETRY','failure':failure,'complete_worker_wall_seconds':time.monotonic()-t,'complete_worker_CPU_seconds':after.ru_utime+after.ru_stime-before.ru_utime-before.ru_stime,'children_reaped':True,'worker_count_requested':n,'max_active_workers':maxactive,'completed':completed,'scope':'all local audit workers launch through exit/reap and final source pins; scheduler receipt write excluded','new_native_oracle_calls':0,'old_baseline_Arb_or_validator_replay':0,'time_limit':None})
 if failure:raise SystemExit(1)
if __name__=='__main__':main()
