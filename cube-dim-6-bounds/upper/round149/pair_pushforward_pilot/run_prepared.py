#!/usr/bin/env python3
"""Dual-gated NEW canonical density pilot; executed only by the designated operator."""
from pathlib import Path
from fractions import Fraction as F
import concurrent.futures as cf,hashlib,json,multiprocessing as mp,os,resource,sys,time,traceback
sys.dont_write_bytecode=True
P=Path(__file__).resolve().parent
STOP=None
read=lambda p:json.loads(Path(p).read_text())
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
def need(v,msg):
 if not v:raise ArithmeticError(msg)
def write(path,x):
 with Path(path).open('x')as f:json.dump(x,f,indent=2);f.write('\n');f.flush();os.fsync(f.fileno())
def canonical(cols):
 rows=[tuple(map(F,v))for v in cols];return min(tuple(sorted(rows)),tuple(sorted((b,a)for a,b in rows)))
def key(cols):return hashlib.sha256((json.dumps([list(map(str,v))for v in canonical(cols)],separators=(',',':'))+'\n').encode()).hexdigest()
def release_requirements():
 return {'round':'149','authorized':True,'scope':'AT_MOST_100_NEW_CANONICAL_PUSHFORWARD_DENSITIES_NO_K_INTEGRAL','max_workers':15,'max_new_cases':100,'max_stage_checks_total':500,'max_literal_endpoint_substitutions_total':220000,'time_limit':None,'nice_minimum':10,'new_k_integral_special_function_native_kernel_calls_permitted':0,'old_density_keys_recompiled_permitted':0,'new_beta2_or_C6_claim_permitted':False,'run_sha256':sha(__file__),'core_sha256':sha(P/'density_core.py'),'prepare_sha256':sha(P/'prepare.py'),'proof_sha256':sha(P/'PROOF.md'),'case_packet_sha256':sha(P/'case_packet.json'),'source_manifest_sha256':sha(P/'source_manifest.json'),'registration_sha256':sha(P/'registration.json')}
def preflight():
 reg=read(P/'registration.json');packet=read(P/'case_packet.json');manifest=read(P/'source_manifest.json')
 for file,k in [('run_prepared.py','run_sha256'),('density_core.py','core_sha256'),('prepare.py','prepare_sha256'),('PROOF.md','proof_sha256'),('case_packet.json','case_packet_sha256'),('source_manifest.json','source_manifest_sha256')]:need(sha(P/file)==reg[k],'FROZEN_FILE:'+file)
 for row in manifest['sources']:need(sha(row['path'])==row['sha256']and Path(row['path']).stat().st_size==row['bytes'],'SOURCE_HASH')
 need(sha(packet['terra_pattern_path'])==packet['terra_pattern_sha256']and sha(packet['cone_table_path'])==packet['cone_table_sha256'],'SOURCE_TABLES')
 cases=packet['selected_cases'];excluded={r['canonical_key']for r in packet['excluded_passed_references']};need(len(cases)==100 and len(excluded)==3,'FIXED_COUNT_AND_REUSE_KEYS')
 need(len({r['canonical_key']for r in cases})==100 and not excluded&{r['canonical_key']for r in cases},'NO_RECOMPILATION_OR_DUPLICATE')
 patterns=read(packet['terra_pattern_path'])['patterns'];cones=read(packet['cone_table_path'])['cones'];sourcecount=0;weight=F(0)
 for case in cases:
  need(key(case['columns'])==case['canonical_key'],'CANONICAL_IDENTITY');original=patterns[case['terra_pattern_index']]
  need(canonical(original['columns'])==canonical(case['columns'])and original['sources']==case['sources']and original['source_count']==case['source_count']and F(original['weight'])==F(case['pattern_weight']),'EXACT_TERRA_PATTERN')
  w=F(0)
  for cid,a,b in case['sources']:
   need(key(zip(cones[cid]['nodes_by_body'][a],cones[cid]['nodes_by_body'][b]))==case['canonical_key'],'EXACT_ORIGINAL_PAIR');w+=F(cones[cid]['alpha'])
  need(w==F(case['pattern_weight']),'FULL_SOURCE_WEIGHT');weight+=w;sourcecount+=case['source_count']
 need(weight==F(packet['selected_pattern_weight'])and sourcecount==packet['selected_source_pair_count'],'SELECTED_WEIGHT_ACCOUNTING')
 need({r['geometry']['active_columns']for r in cases}==set(range(2,7))and max(r['geometry']['distinct_directions']for r in cases)==6,'COVER_M2_TO_M6_AND_MAX_DIRECTIONS')
 need(sum(r['geometry']['active_columns']-1 for r in cases)==reg['max_stages_total']<=500,'FINITE_STAGE_TOTAL')
 need(packet['caps']==reg['caps_per_case']and reg['max_literal_endpoint_substitutions_total']==220000,'CAPS_BINDING')
 need(packet['resonance_witness_key']in{r['canonical_key']for r in cases},'RESERVED_RESONANCE_WITNESS')
 return packet,reg

def init_worker(stop):
 global STOP
 STOP=stop
 for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS']:os.environ[k]='1'
 sys.dont_write_bytecode=True
 os.nice(max(0,10-os.getpriority(os.PRIO_PROCESS,0)))

def worker(case,caps,release_hash,frozen,resonance_key):
 outdir=P/'case_outputs'/case['case_id'];outdir.mkdir(exist_ok=False);start=time.monotonic();cpu=time.process_time();checkpoints=[];counter_state={};compilation_started=False
 intent={'case_id':case['case_id'],'canonical_key':case['canonical_key'],'columns':case['columns'],'pattern_weight':case['pattern_weight'],'sources':case['sources'],'caps':caps,'release_sha256':release_hash,'pid':os.getpid(),'nice':os.getpriority(os.PRIO_PROCESS,0),'time_limit':None}
 write(outdir/'intent.json',intent)
 if STOP.is_set():
  result={'status':'STOPPED_BEFORE_NEW_CASE','case_id':case['case_id'],'canonical_key':case['canonical_key'],'new_density_generations':0,'density_generation_started':False};write(outdir/'receipt.json',result);return result
 try:
  for n,h in frozen.items():need(sha(P/n)==h,'WORKER_SOURCE_FREEZE')
  write(outdir/'started.json',{'pid':os.getpid(),'nice':os.getpriority(os.PRIO_PROCESS,0),'intent_sha256':sha(outdir/'intent.json'),'time_limit':None})
  from density_core import compile_case
  def checkpoint(stage):
   name='stage_'+str(stage['active_columns'])+'.json';write(outdir/name,stage);checkpoints.append({'path':name,'sha256':sha(outdir/name),'bytes':(outdir/name).stat().st_size})
  compilation_started=True;result=compile_case(case,caps,checkpoint,counter_state)
  if case['canonical_key']==resonance_key:need(result['actual_counts']['resonant_branches']>0,'M3_DUPLICATE_COLUMN_RESONANCE_COVERAGE')
  result.update({'stage_checkpoints':checkpoints,'pattern_weight':case['pattern_weight'],'source_count':case['source_count'],'density_generation_started':True,'case_output_pid':os.getpid(),'wall_seconds':time.monotonic()-start,'CPU_seconds':time.process_time()-cpu})
  for n,h in frozen.items():need(sha(P/n)==h,'WORKER_SOURCE_FREEZE_AT_COMPLETION')
  write(outdir/'result.json',result);write(outdir/'receipt.json',{'status':result['status'],'canonical_key':case['canonical_key'],'result_sha256':sha(outdir/'result.json'),'actual_counts':result['actual_counts'],'wall_seconds':result['wall_seconds'],'CPU_seconds':result['CPU_seconds'],'new_k_integral_evaluations':0,'time_limit':None})
  return {'status':result['status'],'case_id':case['case_id'],'canonical_key':case['canonical_key'],'result_path':str((outdir/'result.json').relative_to(P)),'result_sha256':sha(outdir/'result.json'),'actual_counts':result['actual_counts'],'density_generation_started':True,'case_output_pid':os.getpid(),'CPU_seconds':result['CPU_seconds']}
 except BaseException:
  STOP.set();failed={'status':'FAILED_STOP_NO_RETRY','case_id':case['case_id'],'canonical_key':case['canonical_key'],'failure':traceback.format_exc(),'completed_stage_checkpoints':checkpoints,'actual_counts':counter_state,'density_generation_started':compilation_started,'case_output_pid':os.getpid(),'wall_seconds':time.monotonic()-start,'CPU_seconds':time.process_time()-cpu,'new_k_integral_evaluations':0,'automatic_retry':False};write(outdir/'receipt.json',failed);return failed

def run():
 rel=read(P/'parent_release.json')
 for k,v in release_requirements().items():need(rel.get(k)==v,'RELEASE:'+k)
 need(rel.get('main_approved')is True and rel.get('independent_review_approved')is True,'DUAL_RELEASE_REQUIRED');workers=rel.get('workers');need(type(workers)is int and 1<=workers<=15,'MAIN_WORKERS_1_TO_15')
 need(not(P/'started.json').exists(),'NO_AUTOMATIC_RESTART');packet,reg=preflight();os.nice(max(0,10-os.getpriority(os.PRIO_PROCESS,0)))
 need(all(os.environ.get(k)=='1'for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS']),'ONE_THREAD_ENV')
 releasehash=sha(P/'parent_release.json');frozen={n:sha(P/n)for n in ['run_prepared.py','density_core.py','PROOF.md','case_packet.json','registration.json','source_manifest.json']}
 write(P/'started.json',{'pid':os.getpid(),'nice':os.getpriority(os.PRIO_PROCESS,0),'workers':workers,'release_sha256':releasehash,'new_cases':100,'time_limit':None});(P/'case_outputs').mkdir(exist_ok=False)
 start=time.monotonic();cpu=time.process_time();childbefore=resource.getrusage(resource.RUSAGE_CHILDREN);done=[];launched=0;pending={};stop=None;failure=None
 try:
  context=mp.get_context('spawn');stop=context.Event();cases=packet['selected_cases'];nextindex=0
  with cf.ProcessPoolExecutor(max_workers=workers,mp_context=context,initializer=init_worker,initargs=(stop,))as pool:
   def fill():
    nonlocal nextindex,launched
    while len(pending)<workers and nextindex<len(cases)and not stop.is_set():
     case=cases[nextindex];nextindex+=1;pending[pool.submit(worker,case,packet['caps'],releasehash,frozen,packet['resonance_witness_key'])]=case['case_id'];launched+=1
   fill()
   while pending:
    ready,_=cf.wait(pending,return_when=cf.FIRST_COMPLETED)
    for future in ready:
     case_id=pending.pop(future)
     try:record=future.result()
     except BaseException:stop.set();record={'status':'WORKER_PROCESS_FAILURE','case_id':case_id,'failure':traceback.format_exc()}
     done.append(record)
     write(P/('completion_'+str(len(done)).zfill(3)+'.json'),record)
     if record['status']!='PASS_EXACT_DENSITY_CASE':stop.set()
    fill()
  # Exiting the executor reaps every original worker; in-flight cases preserve their independent checkpoint/result.
  success=len(done)==100 and all(r['status']=='PASS_EXACT_DENSITY_CASE'for r in done)
  counts_complete=all(v['status']!='WORKER_PROCESS_FAILURE'for v in done)
  stage_total=sum(r.get('actual_counts',{}).get('density_stage_checks',0)for r in done);endpoint_total=sum(r.get('actual_counts',{}).get('literal_endpoint_substitutions',0)for r in done)
  need(stage_total<=500 and endpoint_total<=220000,'GLOBAL_FINITE_CAPS')
  if success:need(stage_total==reg['max_stages_total'],'EXACT_TOTAL_STAGE_COUNT')
  for n,h in frozen.items():need(sha(P/n)==h,'SOURCE_FREEZE_AT_COMPLETION')
  out={'status':'COMPLETE_100_NEW_CANONICAL_DENSITIES'if success else'PARTIAL_STOP_NO_RETRY','new_cases_launched':launched,'density_generations_started':sum(bool(v.get('density_generation_started'))for v in done),'passed_new_cases':sum(r['status']=='PASS_EXACT_DENSITY_CASE'for r in done),'worker_limit':workers,'all_launched_workers_reaped':True,'case_records':done,'stage_checks_total':stage_total if counts_complete else None,'literal_endpoint_substitutions_total':endpoint_total if counts_complete else None,'observed_stage_checks':stage_total,'observed_literal_endpoint_substitutions':endpoint_total,'per_case_counts_complete':counts_complete,'unobserved_case_count_ids':[v['case_id']for v in done if v['status']=='WORKER_PROCESS_FAILURE'],'literal_endpoint_budget_upper_for_submitted_cases':launched*2200,'excluded_passed_keys_referenced_only':packet['excluded_passed_references'],'selected_pattern_weight':packet['selected_pattern_weight'],'selected_source_pair_count':packet['selected_source_pair_count'],'passed_pattern_weight':str(sum((F(v['pattern_weight'])for v in packet['selected_cases']if v['case_id']in{r['case_id']for r in done if r['status']=='PASS_EXACT_DENSITY_CASE'}),F(0))),'distinct_observed_worker_pids':sorted({v['case_output_pid']for v in done if 'case_output_pid'in v}),'new_k_integral_evaluations':0,'old_6D_grid_replays':0,'new_beta2_upper':False,'new_C6_bound':False,'scope':'density law compiler only; no k integral, beta2 value or new C6 result','run_sha256':sha(__file__),'core_sha256':sha(P/'density_core.py'),'release_sha256':releasehash,'wall_seconds':time.monotonic()-start,'controller_CPU_seconds':time.process_time()-cpu}
  write(P/'result.json',out)
  if not success:failure='Some new cases failed/stopped; no retry; independent successful checkpoints preserved'
 except BaseException:
  failure=traceback.format_exc()
  if stop is not None:stop.set()
  raise
 finally:
  childafter=resource.getrusage(resource.RUSAGE_CHILDREN);write(P/'exit_receipt.json',{'status':'COMPLETE_100_NEW_CANONICAL_DENSITIES'if failure is None else'FAILED_OR_PARTIAL_STOP_NO_RETRY','result_sha256':sha(P/'result.json')if(P/'result.json').exists()else None,'launched_case_count':launched,'completed_record_count':len(done),'failure':failure,'wall_seconds':time.monotonic()-start,'controller_CPU_seconds':time.process_time()-cpu,'complete_reaped_child_CPU_seconds':childafter.ru_utime+childafter.ru_stime-childbefore.ru_utime-childbefore.ru_stime,'time_limit':None,'automatic_restarts':0,'new_k_integral_evaluations':0,'new_beta2_upper':False,'new_C6_bound':False})
 if failure:raise RuntimeError(failure)
 print(json.dumps({'status':'COMPLETE_100_NEW_CANONICAL_DENSITIES','new_cases':100,'stage_checks_total':stage_total if counts_complete else None,'literal_endpoint_substitutions_total':endpoint_total if counts_complete else None,'observed_stage_checks':stage_total,'observed_literal_endpoint_substitutions':endpoint_total,'per_case_counts_complete':counts_complete,'unobserved_case_count_ids':[v['case_id']for v in done if v['status']=='WORKER_PROCESS_FAILURE'],'literal_endpoint_budget_upper_for_submitted_cases':launched*2200,'new_k_integral_evaluations':0}))
if __name__=='__main__':
 if sys.argv[1:]==['--saved-only-preflight']:
  packet,_=preflight();print(json.dumps({'status':'SAVED_ONLY_PREFLIGHT_PASS','new_cases':len(packet['selected_cases']),'excluded_passed_keys':3,'new_density_generations':0,'new_k_integral_evaluations':0}))
 elif not sys.argv[1:]:run()
 else:raise SystemExit('usage: run_prepared.py [--saved-only-preflight]')
