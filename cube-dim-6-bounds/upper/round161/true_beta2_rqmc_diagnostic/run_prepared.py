#!/usr/bin/env python3
"""Terra-only one finite NEW-sample batch, dual-gated; no automatic retries."""
import sys
sys.dont_write_bytecode=True
import concurrent.futures as cf,multiprocessing as mp,os,platform,resource,time,traceback,math
from fractions import Fraction as F
from common import P,need,read,sha,write_new,check_sources,check_freeze
STOP=None
THREAD_ENV=('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS')

def requirements():
 reg=read(P/'registration.json')
 return {'round':'C161_V2_TRUE_BETA2_RQMC_DIAGNOSTIC','authorized':True,'scope':'ONE_FIXED_NEW_SAMPLE_NUMERICAL_DIAGNOSTIC_ONLY','registration_sha256':sha(P/'registration.json'),'case_packet_sha256':sha(P/'case_packet.json'),'source_freeze':reg['source_freeze'],'jobs':168,'replicates':8,'max_workers':15,'nice_minimum':19,'time_limit':None,'automatic_retries':0,'density_compilations':0,'old_grid_or_MC_replays':0,'strict_mu2_bound_claim_allowed':False,'new_Cn_claim_allowed':False}
def preflight():
 reg=read(P/'registration.json');packet=read(P/'case_packet.json');pins=read(P/'source_pins.json');plan=read(P/'plan.json')
 check_freeze(reg['source_freeze']);check_sources(pins)
 need(sha(P/'case_packet.json')==reg['case_packet_sha256']and sha(P/'plan.json')==reg['plan_sha256']and sha(P/'source_pins.json')==reg['source_pins_sha256'],'FROZEN_PREPARED_INPUTS')
 need(reg['source_freeze']==packet['source_freeze']and packet['source_pins_sha256']==reg['source_pins_sha256'],'PREPARED_SOURCE_CHAIN')
 need(packet['cones']==read(pins['sources']['130_cone_table']['path'])['cones'],'UNCHANGED_ALL21_SOURCE_CONES')
 need(packet['h']=='1/128'and len(packet['jobs'])==168 and len({j['seed']for j in packet['jobs']})==168,'FIXED_H_AND_DISTINCT_SEEDS')
 need([(j['replicate'],j['cone_id'])for j in packet['jobs']]==[(r,c)for r in range(8)for c in range(21)],'EXACT168_JOB_PARTITION')
 need(reg['finite_counts']==plan['finite_counts']and reg['versions']==packet['versions']and platform.python_version()==reg['python_version'],'FINITE_REGISTRATION_AND_RUNTIME')
 return packet,reg

def init_worker(stop,frozen):
 global STOP
 STOP=stop
 for name in THREAD_ENV:os.environ[name]='1'
 os.nice(max(0,19-os.getpriority(os.PRIO_PROCESS,0)));check_freeze(frozen)

def worker(job,packet,release_hash):
 out=P/'job_outputs'/job['job_id'];out.mkdir(exist_ok=False)
 start=time.monotonic();cpu=time.process_time();counters={};sample_state={'phase':'PLANNED_NO_SOBOL_REQUEST','planned_points':16384,'sobol_constructor_attempted':False,'sobol_request_attempted':False,'sobol_request_completed':False,'actual_sample_points':0,'actual_sample_points_known':True}
 write_new(out/'intent.json',{'job':job,'pid':os.getpid(),'nice':os.getpriority(os.PRIO_PROCESS,0),'release_sha256':release_hash,'time_limit':None,'finite_points':16384,'all_body_forms_share_E':True})
 def checkpoint(name,value):write_new(out/name,value)
 if STOP.is_set():
  result={'status':'STOPPED_BEFORE_NEW_SAMPLES','job':job,'worker_pid':os.getpid(),'sample_accounting':dict(sample_state)};write_new(out/'receipt.json',result);return result
 try:
  check_freeze(packet['source_freeze']);write_new(out/'started.json',{'pid':os.getpid(),'time_limit':None})
  from rqmc_core import evaluate_job
  result=evaluate_job(job,packet['cones'][job['cone_id']],packet,checkpoint,counters,sample_state)
  check_freeze(packet['source_freeze']);result.update({'worker_pid':os.getpid(),'nice':os.getpriority(os.PRIO_PROCESS,0),'sample_accounting':dict(sample_state),'versions':packet['versions']})
  write_new(out/'result.json',result)
  record={'status':result['status'],'job':job,'result_path':str((out/'result.json').relative_to(P)),'result_sha256':sha(out/'result.json'),'worker_pid':os.getpid(),'sample_accounting':dict(sample_state),'actual_counts':counters,'wall_seconds':result['wall_seconds'],'CPU_seconds':result['CPU_seconds']}
  write_new(out/'receipt.json',record);return record
 except BaseException:
  if sample_state['sobol_request_attempted']and not sample_state['actual_sample_points_known']:
   counters['sample_points']=None;counters['uniform_coordinates']=None
  STOP.set();record={'status':'FAILED_STOP_NO_RETRY','job':job,'sample_accounting':dict(sample_state),'worker_pid':os.getpid(),'failure':traceback.format_exc(),'actual_counts':counters,'wall_seconds':time.monotonic()-start,'CPU_seconds':time.process_time()-cpu,'hard_interrupt_uncheckpointed_counts':'UNKNOWN','automatic_retry':False}
  write_new(out/'receipt.json',record);return record

def aggregate(packet,records):
 # All168 immutable results are required; do not substitute partial cone sets.
 rows=[]
 for record in records:
  path=P/record['result_path'];need(sha(path)==record['result_sha256'],'RESULT_IDENTITY');rows.append(read(path))
 rows.sort(key=lambda r:(r['job']['replicate'],r['job']['cone_id']))
 need([(r['job']['replicate'],r['job']['cone_id'])for r in rows]==[(r,c)for r in range(8)for c in range(21)],'COMPLETE168_RESULTS')
 kL,kU=map(F,packet['kappa_source_rational_enclosure']);kappa=float((kL+kU)/2)
 reps=[kappa*math.fsum(rows[r*21+c]['weighted_cone_mean_float']for c in range(21))for r in range(8)]
 estimate=math.fsum(reps)/8;sd=math.sqrt(math.fsum((x-estimate)**2 for x in reps)/7)
 counts={name:sum(row['actual_counts'][name]for row in rows)for name in rows[0]['actual_counts']}
 expected=read(P/'plan.json')['finite_counts']
 for key,value in expected.items():need(counts[key]==value,'EXACT_FINITE_COUNT:'+key)
 return {'status':'COMPLETE_FIXED_FINITE_RQMC_DIAGNOSTIC','mu2_estimate_float':estimate,'replicate_estimates_float':reps,'empirical_replicate_SD_float':sd,'replicate_min_float':min(reps),'replicate_max_float':max(reps),'empirical_replicate_SD_over_sqrt8_float':sd/math.sqrt(8),'dispersion_interpretation':'descriptive8-scramble dispersion only; no confidence level, no rigorous error bound','kappa_symbolic':packet['kappa_symbolic'],'kappa_source_rational_enclosure':packet['kappa_source_rational_enclosure'],'kappa_float_hex':kappa.hex(),'full168_job_records':records,'actual_counts':counts,'U_hashes':[{'job_id':row['job']['job_id'],'U_sha256':row['U_sha256']}for row in rows],'new_density_generations':0,'old_grid_MC_replays':0,'strict_mu2_lower':False,'strict_mu2_upper':False,'confidence_proof':False,'new_Cn_bound':False,'automatic_science_expansion_or_stop151':False,'scope':'fixed h1/128 original21cone32body shared Exp6 law; new scrambled Sobol samples only'}

def run():
 release=read(P/'parent_release.json')
 for key,value in requirements().items():need(release.get(key)==value,'SCIENCE_RELEASE:'+key)
 need(release.get('main_approved')is True and release.get('independent_review_approved')is True,'SCIENCE_DUAL_GATE')
 workers=release.get('workers');need(type(workers)is int and 1<=workers<=15,'MAX15_SINGLE_THREAD_WORKERS')
 need(not(P/'started.json').exists(),'NO_RESTART')
 packet,reg=preflight();need(all(os.environ.get(name)=='1'for name in THREAD_ENV),'ONE_THREAD_ENV_REQUIRED')
 os.nice(max(0,19-os.getpriority(os.PRIO_PROCESS,0)))
 start=time.monotonic();cpu=time.process_time();before=resource.getrusage(resource.RUSAGE_CHILDREN);records=[];issued=[];pending={};failure=None;pool_shutdown_observed=None;shutdown_failure=None
 write_new(P/'started.json',{'pid':os.getpid(),'nice':os.getpriority(os.PRIO_PROCESS,0),'workers':workers,'jobs':168,'release_sha256':sha(P/'parent_release.json'),'time_limit':None,'thread_environment':{name:os.environ[name]for name in THREAD_ENV},'resource_start':{'loadavg':list(os.getloadavg()),'physical_pages':os.sysconf('SC_PHYS_PAGES')if 'SC_PHYS_PAGES'in os.sysconf_names else None,'page_size':os.sysconf('SC_PAGE_SIZE')}})
 (P/'job_outputs').mkdir(exist_ok=False);context=mp.get_context('spawn');stop=context.Event();nextjob=0
 try:
  pool=cf.ProcessPoolExecutor(max_workers=workers,mp_context=context,initializer=init_worker,initargs=(stop,packet['source_freeze']))
  try:
   def fill():
    nonlocal nextjob
    while len(pending)<workers and nextjob<len(packet['jobs'])and not stop.is_set():
     job=packet['jobs'][nextjob]
     # Write-ahead issue is durable before submit; crash gaps remain UNKNOWN/no retry.
     write_new(P/('issued_'+str(nextjob).zfill(3)+'.json'),{'state':'ISSUED_PREPARE','job':job,'cannot_be_reclaimed_from_missing_worker_intent':True})
     issued.append(job);nextjob+=1
     future=pool.submit(worker,job,packet,sha(P/'parent_release.json'));pending[future]=job
   fill()
   while pending:
    ready,_=cf.wait(pending,timeout=10,return_when=cf.FIRST_COMPLETED)
    for future in ready:
     job=pending.pop(future)
     try:record=future.result()
     except BaseException:
      stop.set();record={'status':'WORKER_PROCESS_FAILURE_UNKNOWN_COUNTS','job':job,'failure':traceback.format_exc(),'uncheckpointed_counts':'UNKNOWN','sample_accounting':None,'actual_sample_points':'UNKNOWN_FOR_UNRETURNED_JOB'}
     records.append(record);write_new(P/('completion_'+str(len(records)).zfill(3)+'.json'),record)
     if record['status']!='COMPLETE_RQMC_CONE_DIAGNOSTIC':stop.set()
    fill()
  except BaseException:
   # Set stop BEFORE waiting for pool shutdown, including interruption in the body.
   stop.set();raise
  finally:
   try:
    pool.shutdown(wait=True,cancel_futures=False)
   except BaseException:
    stop.set();shutdown_failure=traceback.format_exc();pool_shutdown_observed=None;raise
   else:pool_shutdown_observed=True
  complete=len(records)==168 and all(r['status']=='COMPLETE_RQMC_CONE_DIAGNOSTIC'for r in records)
  check_freeze(reg['source_freeze']);check_sources(read(P/'source_pins.json'))
  if complete:
   result=aggregate(packet,records);write_new(P/'result.json',result)
  else:
   failure='partial or failed jobs; no retry, no substitute incomplete moment estimate'
   write_new(P/'partial_result.json',{'status':'PARTIAL_FAILED_NO_RETRY','issued_jobs':issued,'records':records,'mu2_estimate':None,'strict_bound':False})
 except BaseException:
  stop.set();failure=traceback.format_exc();raise
 finally:
  after=resource.getrusage(resource.RUSAGE_CHILDREN)
  write_new(P/'exit_receipt.json',{'status':'COMPLETE_FIXED_FINITE_RQMC_DIAGNOSTIC'if failure is None and len(records)==168 else'PARTIAL_FAILED_NO_RETRY','issued_count':len(issued),'completed_record_count':len(records),'failure':failure,'wall_seconds':time.monotonic()-start,'controller_CPU_seconds':time.process_time()-cpu,'observed_reaped_child_CPU_seconds':after.ru_utime+after.ru_stime-before.ru_utime-before.ru_stime,'child_CPU_scope':'complete_original_pool_only_if_shutdown_returned_wait_true; otherwise observed subset only','all_original_pool_workers_reaped_by_shutdown_return':pool_shutdown_observed,'pool_shutdown_observation':'RETURNED_WAIT_TRUE'if pool_shutdown_observed is True else'UNKNOWN_INTERRUPTED_OR_UNOBSERVED','shutdown_failure':shutdown_failure,'worker_pids_from_returned_records':sorted({r['worker_pid']for r in records if 'worker_pid'in r}),'submitted_job_ids_without_returned_record':[j['job_id']for j in issued if j['job_id']not in {r['job']['job_id']for r in records}],'unobserved_scope':'Worker PIDs may also exist only in per-job intent files; missing returns never imply zero samples. Independent operator verifies terminal/reap; no incomplete job is a mathematical false case','controller_maxrss':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,'loadavg_end':list(os.getloadavg()),'result_sha256':sha(P/'result.json')if(P/'result.json').exists()else None,'time_limit':None,'automatic_retries':0,'strict_bound':False})
 if failure:raise RuntimeError(failure)
 print('COMPLETE_FIXED_FINITE_RQMC_DIAGNOSTIC; descriptive only; no rigorous mu2 bound')
if __name__=='__main__':
 if sys.argv[1:]==['--release-requirements']:
  import json
  print(json.dumps(requirements(),indent=2))
 elif not sys.argv[1:]:run()
 else:raise SystemExit('usage: run_prepared.py [--release-requirements]')
