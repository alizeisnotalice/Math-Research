"""Terra-only v2 immutable-input audit; complete observed API state accounting."""
from pathlib import Path
import json,hashlib,os,sys,time,traceback,resource,copy
sys.dont_write_bytecode=True
P=Path(__file__).resolve().parent
from api_ledger import new_ledger
read=lambda p:json.loads(Path(p).read_text())
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
def need(v,m):
 if not v:raise ArithmeticError(m)
def write(n,v):
 with(P/n).open('x')as f:json.dump(v,f,indent=2);f.write('\n');f.flush();os.fsync(f.fileno())
def release_requirements():
 reg=read(P/'registration.json')
 return {'round':'155_V2','authorized':True,'scope':'SAME_FIRST_FIVE_IMMUTABLE_NEW149_CASES_EXACT_INTEGER_AUDIT_NO_RESELECTION','workers':1,'cases':5,'stages':15,'nice_minimum':15,'time_limit':None,'new_density_generations_permitted':0,'new_k_or_legacy_native_integral_calls_permitted':0,'primary_or_SymPy_baseline_replays_permitted':0,'registration_sha256':sha(P/'registration.json'),'case_packet_sha256':sha(P/'case_packet.json'),'source_manifest_sha256':sha(P/'source_manifest.json'),'provenance_sha256':sha(P/'provenance.json'),'source_sha256':reg['source_sha256'],'caps':reg['caps'],'counter_budget_semantics':reg['counter_budget_semantics'],'API_state_semantics':reg['API_state_semantics']}
def quiescent(state):
 return state['inflight_api']is None and state['totals']['attempted']==state['totals']['completed']+state['totals']['raised_exception']and all(v['attempted']==v['completed']+v['raised_exception']for v in state['by_category'].values())and all(state['totals'][key]==sum(v[key]for v in state['by_category'].values())for key in ['attempted','completed','raised_exception'])
def run():
 rel=read(P/'parent_release.json');reg=read(P/'registration.json');packet=read(P/'case_packet.json');manifest=read(P/'source_manifest.json')
 for k,v in release_requirements().items():need(rel.get(k)==v,'RELEASE:'+k)
 need(rel.get('main_approved')is True and rel.get('independent_review_approved')is True,'DUAL_GATE');need(not(P/'started.json').exists(),'NO_AUTOMATIC_RESTART')
 for n,h in reg['source_sha256'].items():need(sha(P/n)==h,'SOURCE_PIN')
 need(sha(P/'case_packet.json')=='6587f978228bc7478f25919425d0ed915ff24270786f4a9961027da17b3d495b'==reg['case_packet_sha256'],'FIRST_V1_PACKET_UNCHANGED')
 need(sha(P/'source_manifest.json')=='3012a7f9a22c838334158228ae2e44c8e2e4f84e10acdaf2fa6039b626e34d9b'==reg['source_manifest_sha256'],'FIRST_V1_SOURCE_MANIFEST_UNCHANGED')
 need(sha(P/'provenance.json')==reg['provenance_sha256'],'PROVENANCE_PIN')
 for row in manifest['sources']:need(sha(row['path'])==row['sha256']and Path(row['path']).stat().st_size==row['bytes'],'ORIGINAL_SAVED_INPUT_PIN')
 need(packet['caps']==reg['caps']and len(packet['selected_cases'])==5 and sum(len(c['stages'])for c in packet['selected_cases'])==15,'FIXED_SCOPE')
 need(all(os.environ.get(k)=='1'for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS']),'SINGLE_THREAD_ENV')
 os.nice(max(0,15-os.getpriority(os.PRIO_PROCESS,0)));write('started.json',{'pid':os.getpid(),'nice':os.getpriority(os.PRIO_PROCESS,0),'release_sha256':sha(P/'parent_release.json'),'case_packet_sha256':sha(P/'case_packet.json'),'cases':5,'stages':15,'process_crash_semantics':'uncheckpointed API totals and inflight call status UNKNOWN, never0'})
 start=time.monotonic();cpu=time.process_time();ledger=new_ledger();completed=[];current=None
 try:
  from verifier_core import verify_stage
  for case in packet['selected_cases']:
   rows=[]
   for entry in case['stages']:
    current={'case_id':case['case_id'],'active_columns':entry['active_columns'],'source_stage_path':entry['path'],'source_stage_sha256':entry['sha256']}
    name=case['case_id']+'_stage_'+str(entry['active_columns']);write(name+'_intent.json',{**current,'state':'STAGE_INFLIGHT_API_AND_PROOF_UNKNOWN','completed_stage_count_before':len(completed),'API_global_observed_before':copy.deepcopy(ledger)})
    need(sha(entry['path'])==entry['sha256'],'STAGE_IMMUTABLE')
    row=verify_stage(read(entry['path']),case['columns'],packet['caps'],ledger);need(quiescent(ledger),'QUIESCENT_GLOBAL_API_COUNTS');row.update(current);row['API_global_ledger_snapshot']=copy.deepcopy(ledger)
    write(name+'.json',row);completed.append({'path':name+'.json','sha256':sha(P/(name+'.json')),**current});rows.append(row);current=None
   write(case['case_id']+'_receipt.json',{'status':'PASS_ALL_SAVED_CASE_STAGES','case_id':case['case_id'],'canonical_key':case['canonical_key'],'saved_stage_checks':len(rows),'verifier_CPU_seconds':sum(r['CPU_seconds']for r in rows),'primary_saved_whole_case_CPU_seconds':case['primary_saved_case_CPU_seconds'],'timings_not_same_workload':True,'no_speedup_claim':True})
  need(len(completed)==15,'EXACT15_CHECKS')
  for row in manifest['sources']:need(sha(row['path'])==row['sha256'],'SAVED_INPUT_PIN_AT_END')
  for n,h in reg['source_sha256'].items():need(sha(P/n)==h,'SOURCE_PIN_AT_END')
  result={'status':'PASS15_V2_EXACT_INTEGER_CHECKS_SAME_FIRST_FIVE_CASES','completed':completed,'API_global_ledger':ledger,'API_totals_exact':True,'API_totals':ledger['totals'],'query_counts_separate':ledger['by_category']['query'],'budgeted_setup_construction_arithmetic_attempts':ledger['budgeted_setup_construction_arithmetic_attempts'],'wall_seconds':time.monotonic()-start,'CPU_seconds':time.process_time()-cpu,'maxrss_platform_units':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,'maxrss_units':'bytes'if sys.platform=='darwin'else'KiB','new_selection_runs':0,'new_density_generations':0,'new_k_integrals':0,'legacy_native_integral_kernel_calls':0,'SymPy_or_primary_baseline_replays':0,'new_beta2_Cn_bound':False,'performance_claim':'independent saved verifier only; no same-workload speedup or live checker replacement'}
  write('result.json',result);write('exit_receipt.json',{'status':result['status'],'result_sha256':sha(P/'result.json'),'API_totals_exact':True,'API_global_ledger':ledger,'completed_stage_count':15,'automatic_restarts':0});print(json.dumps({'status':result['status'],'API_totals':ledger['totals'],'queries':ledger['by_category']['query']}))
 except BaseException:
  exact=quiescent(ledger)and(ledger['observed_stage']is None or quiescent(ledger['observed_stage']['ledger']))
  write('exit_receipt.json',{'status':'FAILED_STOP_NO_RETRY','failure':traceback.format_exc(),'failed_or_inflight_stage':current,'completed':completed,'API_observed_ledger':ledger,'API_totals_exact':exact,'API_totals':ledger['totals']if exact else None,'inflight_or_interrupted_API_status':'NONE_ALL_INVOCATIONS_QUIESCENT'if exact else'UNKNOWN','native_invocations_raised_exception':ledger['totals']['raised_exception'],'new_density_generations':0,'new_k_integrals':0,'legacy_native_integral_kernel_calls':0,'automatic_restarts':0});raise
if __name__=='__main__':run()
