"""C158 Terra-only single-process audit:95 fixed remaining saved cases,5 actual references."""
from pathlib import Path
from fractions import Fraction as F
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
def quiescent(state):
 return state['inflight_api']is None and state['totals']['attempted']==state['totals']['completed']+state['totals']['raised_exception']and all(v['attempted']==v['completed']+v['raised_exception']for v in state['by_category'].values())and all(state['totals'][key]==sum(v[key]for v in state['by_category'].values())for key in ['attempted','completed','raised_exception'])
def release_requirements():
 reg=read(P/'registration.json')
 return {'round':'158','authorized':True,'scope':'FIXED95_REMAINING_IMMUTABLE149_SAVED_CASES_REUSE5_ACTUAL155V2_NO_REPLAY','workers':1,'NEW_audit_cases':95,'NEW_audit_stages':365,'reused_actual_cases':5,'reused_actual_stages':15,'nice_minimum':19,'time_limit':None,'new_density_generations_permitted':0,'new_k_or_legacy_native_integral_calls_permitted':0,'primary_or_SymPy_baseline_replays_permitted':0,'registration_sha256':sha(P/'registration.json'),'case_packet_sha256':sha(P/'case_packet.json'),'coverage_ledger_sha256':sha(P/'coverage_ledger.json'),'source_manifest_sha256':sha(P/'source_manifest.json'),'source_sha256':reg['source_sha256'],'caps':reg['caps'],'scope_failure_is_not_mathematical_family_exclusion':True}
def run():
 rel=read(P/'parent_release.json');reg=read(P/'registration.json');packet=read(P/'case_packet.json');coverage=read(P/'coverage_ledger.json');manifest=read(P/'source_manifest.json')
 for k,v in release_requirements().items():need(rel.get(k)==v,'RELEASE:'+k)
 need(rel.get('main_approved')is True and rel.get('independent_review_approved')is True,'DUAL_GATE');need(not(P/'started.json').exists(),'NO_AUTOMATIC_RESTART')
 for n,h in reg['source_sha256'].items():need(sha(P/n)==h,'SOURCE_PIN')
 for n,k in [('case_packet.json','case_packet_sha256'),('coverage_ledger.json','coverage_ledger_sha256'),('source_manifest.json','source_manifest_sha256')]:need(sha(P/n)==reg[k],'FROZEN_PACKET:'+n)
 for row in manifest['sources']:need(sha(row['path'])==row['sha256']and Path(row['path']).stat().st_size==row['bytes'],'IMMUTABLE_SOURCE_PIN')
 catalog={c['case_id']:c for c in packet['case_catalog']};newids=packet['NEW_audit_case_ids'];reused=packet['reused_case_ids']
 need(len(packet['case_catalog'])==len(catalog)==100 and len({c['canonical_key']for c in catalog.values()})==100,'EXACT100_CASE_KEY_UNIVERSE')
 need(newids==coverage['NEW_audit_case_ids']and reused==coverage['reused_case_ids']and len(newids)==len(set(newids))==95 and len(reused)==len(set(reused))==5 and not set(newids)&set(reused)and set(newids)|set(reused)==set(catalog),'FIXED95_PLUS5_PARTITION')
 need(sum(len(catalog[c]['stages'])for c in newids)==365 and sum(len(catalog[c]['stages'])for c in reused)==15 and packet['caps']==reg['caps']and packet['caps']['native_arithmetic_total']==73000000,'FIXED365_API_SCOPE')
 need(all(os.environ.get(k)=='1'for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS']),'SINGLE_THREAD_ENV')
 os.nice(max(0,19-os.getpriority(os.PRIO_PROCESS,0)));write('started.json',{'pid':os.getpid(),'nice':os.getpriority(os.PRIO_PROCESS,0),'release_sha256':sha(P/'parent_release.json'),'case_packet_sha256':sha(P/'case_packet.json'),'coverage_ledger_sha256':sha(P/'coverage_ledger.json'),'NEW_audit_cases':95,'NEW_audit_stages':365,'reused_actual_cases':5,'process_crash_semantics':'uncheckpointed API counts and inflight proof status UNKNOWN, never0'})
 start=time.monotonic();cpu=time.process_time();ledger=new_ledger();completed=[];passed_cases=[];current=None
 def summary(status,complete=False):
  exact=quiescent(ledger)and(ledger['observed_stage']is None or quiescent(ledger['observed_stage']['ledger']))
  return {'status':status,'full149100_independent_saved_audit_complete':complete,'planned_NEW_audit_cases':95,'planned_NEW_audit_stages':365,'passed_NEW_audit_cases':len(passed_cases),'passed_NEW_audit_case_ids':list(passed_cases),'completed_NEW_stage_checks':len(completed),'completed':list(completed),'reused_actual155_cases':coverage['reused_actual_certificates'],'reused_cases':5,'reused_stage_checks':15,'full149_case_count_if_complete':100 if complete else None,'full149_stage_count_if_complete':380 if complete else None,'failed_or_inflight_stage':current,'API_observed_global_ledger':copy.deepcopy(ledger),'API_totals_exact':exact,'API_totals':copy.deepcopy(ledger['totals'])if exact else None,'NEW_audit_passed_pattern_weight':str(sum((F(catalog[c]['pattern_weight'])for c in passed_cases),F(0))),'planned_NEW_audit_pattern_weight':coverage['NEW_audit_pattern_weight'],'reused_pattern_weight':coverage['reused_pattern_weight'],'full149_pattern_weight':coverage['full_pattern_weight'],'source_case_packet_sha256':sha(P/'case_packet.json'),'coverage_ledger_sha256':sha(P/'coverage_ledger.json'),'wall_seconds':time.monotonic()-start,'CPU_seconds':time.process_time()-cpu,'maxrss_platform_units':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,'maxrss_units':'bytes'if sys.platform=='darwin'else'KiB','scope_failure_is_not_mathematical_family_exclusion':True,'new_density_generations':0,'new_k_integrals':0,'legacy_native_integral_kernel_calls':0,'primary_or_SymPy_baseline_replays':0,'new_beta2_Cn_bound':False,'global9961_density_coverage_claim':False,'performance_claim':'independent exact saved verification only; no same-workload producer speedup'}
 try:
  from verifier_core import verify_stage
  for cid in newids:
   case=catalog[cid];rows=[]
   for entry in case['stages']:
    current={'case_id':cid,'canonical_key':case['canonical_key'],'active_columns':entry['active_columns'],'source_stage_path':entry['path'],'source_stage_sha256':entry['sha256']}
    name=cid+'_stage_'+str(entry['active_columns']);write(name+'_intent.json',{**current,'state':'STAGE_INFLIGHT_API_AND_PROOF_UNKNOWN','completed_stage_count_before':len(completed),'API_observed_before':copy.deepcopy(ledger)})
    need(sha(entry['path'])==entry['sha256'],'SOURCE_STAGE_STILL_IMMUTABLE');row=verify_stage(read(entry['path']),case['columns'],packet['caps'],ledger);need(quiescent(ledger),'QUIESCENT_GLOBAL_API_COUNTS');row.update(current);row['API_global_ledger_snapshot']=copy.deepcopy(ledger)
    write(name+'.json',row);completed.append({'path':name+'.json','sha256':sha(P/(name+'.json')),**current});rows.append(row);current=None
   write(cid+'_receipt.json',{'status':'PASS_ALL_FIXED_SAVED_CASE_STAGES','case_id':cid,'canonical_key':case['canonical_key'],'source_result_sha256':case['source_result_sha256'],'saved_stage_checks':len(rows),'verifier_CPU_seconds':sum(r['CPU_seconds']for r in rows),'primary_saved_whole_case_CPU_seconds':case['primary_saved_whole_case_CPU_seconds'],'pattern_weight':case['pattern_weight'],'source_count':case['source_count'],'no_same_workload_speedup_claim':True});passed_cases.append(cid)
  need(passed_cases==newids and len(completed)==365,'EXACT95_CASE365_STAGE_SUCCESS')
  for row in manifest['sources']:need(sha(row['path'])==row['sha256'],'ALL_ORIGINAL_AND_REUSED_CERT_PINS_AT_END')
  for n,h in reg['source_sha256'].items():need(sha(P/n)==h,'SOURCE_PIN_AT_END')
  for n,k in [('case_packet.json','case_packet_sha256'),('coverage_ledger.json','coverage_ledger_sha256'),('source_manifest.json','source_manifest_sha256')]:need(sha(P/n)==reg[k],'PACKET_PIN_AT_END')
  result=summary('PASS95_REMAINING_PLUS5_REUSED_FULL149100_SAVED_AUDIT',True);write('result.json',result);write('exit_receipt.json',{'status':result['status'],'result_sha256':sha(P/'result.json'),'API_totals_exact':result['API_totals_exact'],'API_observed_global_ledger':ledger,'passed_NEW_audit_cases':95,'completed_NEW_stage_checks':365,'reused_actual_cases':5,'full149100_independent_saved_audit_complete':True,'automatic_restarts':0});print(json.dumps({'status':result['status'],'API_totals':result['API_totals'],'replayed155_cases':0}))
 except BaseException:
  failure=traceback.format_exc();partial=summary('FAILED_SCOPE_AUDIT_INCOMPLETE_NO_FAMILY_EXCLUSION');partial['failure']=failure;write('partial_result.json',partial);write('exit_receipt.json',{'status':partial['status'],'partial_result_sha256':sha(P/'partial_result.json'),'API_totals_exact':partial['API_totals_exact'],'API_totals':partial['API_totals'],'API_observed_global_ledger':ledger,'failed_or_inflight_stage':current,'passed_NEW_audit_cases':len(passed_cases),'completed_NEW_stage_checks':len(completed),'full149100_independent_saved_audit_complete':False,'scope_failure_is_not_mathematical_family_exclusion':True,'automatic_restarts':0,'failure':failure});raise
if __name__=='__main__':run()
