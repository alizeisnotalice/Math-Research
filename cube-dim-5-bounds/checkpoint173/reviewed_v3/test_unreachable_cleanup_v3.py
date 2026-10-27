"""One NEW integrated ZERO-H regression: fault report + truly stuck peer.
No repetition of the v2 three isolated tests. Uses real run_stage and cleanup.
"""
from pathlib import Path
import json,multiprocessing as mp,os,signal,subprocess,sys,time
from util173 import atomic,sha,read
H=Path(__file__).resolve().parent;T=Path('/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/upper_bandcurvature173/reviewed_v3_synthetic')
def synthetic_execute(phase,spec,outroot,manifest_sha,runtime,ready,result,start,stop,resume):
 sid=spec['shard_id'];out=Path(outroot)/phase/f'shard{sid:02d}';out.mkdir(parents=True,exist_ok=True)
 if sid==1:signal.signal(signal.SIGTERM,signal.SIG_IGN)
 atomic(out/'checkpoint.json',{'next_index':1,'last_hash':'synthetic_successful_prefix_only','status':'synthetic_not_math'});ready.put({'shard_id':sid,'pid':os.getpid(),'phase':phase,'synthetic_only':True});start.wait()
 if sid==0:
  # Deliberately report a fault WITHOUT setting stop: exercise the result loop.
  fault={'status':'synthetic_reported_fault','phase':phase,'shard_id':0,'pid':os.getpid(),'synthetic_only':True};atomic(out/'synthetic_fault.json',fault);result.put(fault);return
 if sid==1:
  while True:time.sleep(.05) # truly ignores shared stop and SIGTERM; helper must SIGKILL
 stop.wait();result.put({'status':'complete','phase':phase,'shard_id':sid,'pid':os.getpid(),'synthetic_only':True})

def driver(case):
 import run173 as M
 case=Path(case);meta=case/'meta';work=case/'work';meta.mkdir();work.mkdir();gate=meta/'synthetic_gate_not_qualification.json';atomic(gate,{'synthetic_only':True,'real_source_qualification':False});g={'path':str(gate),'sha256':sha(gate),'accepted':True};manifest=meta/'prepared_manifest.json';atomic(manifest,{'repository_files':{},'source_locks':{},'prepared_external_files':{},'known_qualified168':g});registry=work/'source_registry.json';atomic(registry,{});log=meta/'fake_log_not_mathematical.json';atomic(log,{'synthetic_only':True,'upper':'0'})
 plan={'weighted_band_terms':22295,'sources':455,'fixed_pairs':35872,'max_concurrent_workers':15,'source_registry_path':str(registry),'source_registry_sha256':sha(registry),'old_common_log_path':str(log),'old_common_log_sha256':sha(log),'time_endpoints':[],'phase1_mass_shards':[{'shard_id':i,'synthetic_only':True} for i in range(15)]};atomic(work/'finite_plan.json',plan);auth=meta/'synthetic_authorization_not_real.json';atomic(auth,{'reviewed_manifest_sha256':sha(manifest),'sources':455,'fixed_pairs':35872,'workers':15,'qualified168':g,'qualified172_174':g,'qualified_donor_endpoint_provenance':g,'synthetic_only':True})
 M.H=meta;M.W=work;M.execute=synthetic_execute
 def fake_g_metadata(*args,**kwargs):return [{'synthetic_only':True} for _ in range(15)]
 def forbidden(*args,**kwargs):raise AssertionError('must not enter any real H/band/q/corner function')
 M.original_jets_once=fake_g_metadata;M.assemble_full_domain=forbidden;sys.argv=['integrated_synthetic_only','--dispatch-reviewed','--reviewed-manifest-sha',sha(manifest),'--authorization',str(auth),'--authorization-sha',sha(auth)];raise SystemExit(M.main())

def main():
 assert not T.exists();T.mkdir();case=T/'fault_report_and_stuck_peer';case.mkdir();began=time.monotonic();p=subprocess.Popen([sys.executable,__file__,'--driver',str(case)],stdout=subprocess.PIPE,stderr=subprocess.PIPE,text=True,start_new_session=True)
 try:stdout,stderr=p.communicate(timeout=25) # synthetic watchdog only, never a normal numerical deadline
 except subprocess.TimeoutExpired:
  os.killpg(p.pid,signal.SIGTERM)
  try:stdout,stderr=p.communicate(timeout=2)
  except subprocess.TimeoutExpired:os.killpg(p.pid,signal.SIGKILL);stdout,stderr=p.communicate()
  atomic(T/'regression_fault.json',{'synthetic_watchdog_fired':True,'master_pid':p.pid,'returncode':p.returncode,'stdout':stdout,'stderr':stderr});raise AssertionError('integrated cleanup unreachable; synthetic group cleaned')
 O=case/'work/run_v3_001';summary=read(O/'actual_attempt0001_summary.json');stage=read(O/'actual_attempt0001_mass_exit_receipt.json');cleanup=read(O/'actual_attempt0001_mass_fault_cleanup.json');pool=read(O/'actual_attempt0001_mass_pool_registration.json');assert p.returncode==2 and summary['status']=='fault_prefix_preserved' and summary['candidate'] is None;assert stage['summaries']['0']['status']=='synthetic_reported_fault' and len(stage['summaries'])<15;assert 'enter_cleanup_now' in stage['fault']['exception'];assert cleanup['all_joined'] and len(cleanup['prefix_before_cleanup'])==15;stuck=pool['workers'][1]['pid'];assert any(a['action']=='kill' and a['pid']==stuck for a in cleanup['actions']);assert all(not w['alive'] and w['actual_join_called'] for w in cleanup['workers']);assert next(w['exitcode'] for w in cleanup['workers'] if w['pid']==stuck)==-signal.SIGKILL
 for w in cleanup['workers']:
  try:os.kill(w['pid'],0);raise AssertionError('synthetic worker remains alive')
  except ProcessLookupError:pass
 receipt={'status':'ONE_NEW_ZERO_H_INTEGRATED_REGRESSION_PASSED','test_pid':os.getpid(),'synthetic_master_pid':p.pid,'master_actual_exitcode':p.returncode,'fault_report_received_before_cleanup':True,'summaries_at_cleanup_entry':len(stage['summaries']),'did_not_wait_for_all15_summaries':True,'stuck_peer_pid':stuck,'stuck_peer_ignored_stop_and_SIGTERM':True,'SIGKILL_and_actual_join_confirmed':True,'all15_workers_exited_OS_confirmed':True,'successful_fake_prefixes_preserved':15,'production_cleanup_defaults_seconds':[5,2,2],'normal_numeric_wall_cap_added':False,'new_real_H_jet_donor_band_q_corner_calls':0,'old_v2_three_toys_repeated':False,'wall_seconds':time.monotonic()-began,'stage_receipt_sha256':sha(O/'actual_attempt0001_mass_exit_receipt.json'),'cleanup_receipt_sha256':sha(O/'actual_attempt0001_mass_fault_cleanup.json'),'master_summary_sha256':sha(O/'actual_attempt0001_summary.json'),'synthetic_metadata_is_not_qualification':True};atomic(T/'test_receipt.json',receipt);print(json.dumps(receipt,sort_keys=True),flush=True)
if __name__=='__main__':
 if len(sys.argv)>1 and sys.argv[1]=='--driver':driver(sys.argv[2])
 else:main()
