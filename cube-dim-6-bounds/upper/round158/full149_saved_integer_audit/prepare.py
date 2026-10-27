"""Terra-only C158 input/accounting preparation; no numerical verifier or density import."""
from pathlib import Path
from fractions import Fraction as F
import json,hashlib,os,sys
sys.dont_write_bytecode=True
P=Path(__file__).resolve().parent
S=Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube6-sol-brascamp108/cube-dim-6-bounds/upper/round149/pair_pushforward_pilot')
V=Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube6-sol-brascamp108/cube-dim-6-bounds/upper/round155/exact_laplace_verifier_v2')
PINS149={'source_manifest.json':'fc05fd3a2febd455c07417ebae55c7dbf2b0191b99df8da6d2752eb4e444924a','case_packet.json':'a91d0dad30efdad4290dd2f66376d55474dad805464758a9c4e44b476ba128d1','registration.json':'170c179e59de60b86bbb943afe13e7dd25240872efed9210dc3633b9433f25f6','parent_release.json':'2de7a93029cd1e92531056938bd50c5386bd530c329eb51bcc1bcf40b57de6db','result.json':'e2012d1a0d70912a3e7bd5bf154246f3e2a73da645ba5e2358d588f8a0feae81','exit_receipt.json':'9df86adb3492ef2367f79890c3419bbd60904632739a26e4b93759bdd4657dde','operator_terminal_receipt.json':'b44ac1feb56095e6876388da0eff2c8b92c1b15be8496382ee8d549537b0ae7a'}
PINS155={'case_packet.json':'6587f978228bc7478f25919425d0ed915ff24270786f4a9961027da17b3d495b','result.json':'d65f9ab5b03127334762340d53ed67729ca69ac847e886b06ef9604f6d342eab','verifier_core.py':'ed67062f1bc5d1cfc6fea1ea67a3e92d120e3a08970fd1f224b5a7b104cfd0c6','api_ledger.py':'4ea3048d36ad128130fbe261b6d07c8996c445c98602247f99b1d7b95ba16661'}
read=lambda p:json.loads(Path(p).read_text())
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
def need(v,m):
 if not v:raise ArithmeticError(m)
def write(n,v):
 with(P/n).open('x')as f:json.dump(v,f,indent=2);f.write('\n');f.flush();os.fsync(f.fileno())
def contained(root,rel):
 p=(root/rel).resolve();need(p.is_relative_to(root.resolve()),'CONTAINED_REFERENCE');return p
def main():
 need(not any((P/n).exists()for n in ['case_packet.json','coverage_ledger.json','source_manifest.json','registration.json','prepare_receipt.json','parent_release.json','started.json','result.json']),'NO_PREPARE_RETRY_OR_OVERWRITE')
 pins={}
 def pin(p,expected=None):
  h=sha(p);need(expected is None or h==expected,'FROZEN_SOURCE_PIN');pins[str(p)]={'path':str(p),'sha256':h,'bytes':Path(p).stat().st_size};return h
 for n,h in PINS149.items():pin(S/n,h)
 for n,h in PINS155.items():pin(V/n,h)
 for n in ['verifier_core.py','api_ledger.py']:need(sha(P/n)==PINS155[n],'UNCHANGED_ACTUAL155V2_ALGEBRA_AND_API_LEDGER')
 packet=read(S/'case_packet.json');source_result=read(S/'result.json');source_exit=read(S/'exit_receipt.json');terminal=read(S/'operator_terminal_receipt.json');start=read(S/'started.json');pin(S/'started.json')
 need(source_result['status']==source_exit['status']=='COMPLETE_100_NEW_CANONICAL_DENSITIES'and source_exit['failure']is None and source_exit['result_sha256']==PINS149['result.json'],'REAL149_TERMINAL_SUCCESS')
 need(source_result['passed_new_cases']==source_result['density_generations_started']==source_result['new_cases_launched']==100 and source_result['per_case_counts_complete']is True and source_result['unobserved_case_count_ids']==[]and source_result['all_launched_workers_reaped']is True,'REAL149_COMPLETE100')
 need(source_result['stage_checks_total']==380 and source_result['literal_endpoint_substitutions_total']==3634,'REAL149_STAGE_AND_LITERAL_COUNTS')
 need(terminal['exit_code']==0 and terminal['automatic_restarts']==0 and terminal['per_case_counts_complete']is True and terminal['exact_density_generations_started']==100 and terminal['unobserved_case_ids']==[]and terminal['controller_pid']==start['pid']and terminal['release_sha256']==PINS149['parent_release.json'],'REAL149_OPERATOR_TERMINAL')
 for n,k in [('operator_science_stdout.txt','stdout_sha256'),('operator_science_stderr.txt','stderr_sha256')]:pin(S/n,terminal[k])
 for row in read(S/'source_manifest.json')['sources']:pin(Path(row['path']),row['sha256'])
 pin(S/'source_manifest.json')
 cases=packet['selected_cases'];need(len(cases)==100 and len({c['case_id']for c in cases})==len({c['canonical_key']for c in cases})==100,'EXACT149100_UNIQUE_KEYS')
 source_records={r['case_id']:r for r in source_result['case_records']};need(set(source_records)=={c['case_id']for c in cases}and len(source_result['case_records'])==100,'COMPLETE149_CASE_RECORD_SET')
 catalog=[];stage_index={};literal_total=0
 for case in cases:
  cid=case['case_id'];record=source_records[cid];need(record['status']=='PASS_EXACT_DENSITY_CASE'and record['canonical_key']==case['canonical_key'],'149_CASE_RECORD_IDENTITY')
  resultpath=contained(S,record['result_path']);pin(resultpath,record['result_sha256']);cr=read(resultpath);root=resultpath.parent;receiptpath=root/'receipt.json';pin(receiptpath);receipt=read(receiptpath)
  need(cr['status']==receipt['status']=='PASS_EXACT_DENSITY_CASE'and cr['case_id']==cid and cr['canonical_key']==receipt['canonical_key']==case['canonical_key']and receipt['result_sha256']==sha(resultpath)and cr['actual_counts']==receipt['actual_counts']==record['actual_counts'],'149_CASE_COMPLETE_RESULT_RECEIPT')
  need(F(cr['pattern_weight'])==F(case['pattern_weight'])and cr['source_count']==case['source_count'],'FULL_PATTERN_WEIGHT_AND_SOURCE_COUNT')
  m=case['geometry']['active_columns'];need(len(cr['stage_checkpoints'])==m-1 and cr['actual_counts']['density_stage_checks']==m-1,'ALL_SOURCE_STAGES_PASSED');literal_total+=cr['actual_counts']['literal_endpoint_substitutions'];stages=[]
  for k,entry in zip(range(2,m+1),cr['stage_checkpoints']):
   path=contained(root,entry['path']);pin(path,entry['sha256']);need(path.stat().st_size==entry['bytes'],'STAGE_BYTES');st=read(path);need(st['active_columns']==k and st['Laplace_difference_numerator']=='0'and st['normalization']==1,'SOURCE_COMPLETE_STAGE_ASSERTION')
   row={'path':str(path),'sha256':sha(path),'bytes':path.stat().st_size,'active_columns':k};stages.append(row);stage_index[(cid,k)]=row
  catalog.append({'case_id':cid,'canonical_key':case['canonical_key'],'columns':case['columns'],'geometry':case['geometry'],'pattern_weight':case['pattern_weight'],'source_count':case['source_count'],'sources':case['sources'],'source_result_path':str(resultpath),'source_result_sha256':sha(resultpath),'source_receipt_path':str(receiptpath),'source_receipt_sha256':sha(receiptpath),'primary_saved_whole_case_CPU_seconds':cr['CPU_seconds'],'stages':stages})
 need(len(stage_index)==380 and literal_total==3634,'ENTIRE_SOURCE_STAGE_AND_LITERAL_ACCOUNTING')
 prior=read(V/'result.json');priorpacket=read(V/'case_packet.json');priorexit=read(V/'exit_receipt.json');priorterm=read(V/'operator_terminal_receipt.json');priorstart=read(V/'started.json')
 for n in ['exit_receipt.json','operator_terminal_receipt.json','started.json','parent_release.json','registration.json','provenance.json','source_manifest.json']:pin(V/n)
 need(prior['status']==priorexit['status']=='PASS15_V2_EXACT_INTEGER_CHECKS_SAME_FIRST_FIVE_CASES'and priorexit['result_sha256']==PINS155['result.json']and prior['API_totals_exact']is True and priorexit['API_totals_exact']is True,'ACTUAL155V2_COMPLETE15')
 need(priorterm['exit_code']==0 and priorterm['automatic_retries']==0 and priorterm['pid']==priorstart['pid']and priorterm['release_sha256']==priorstart['release_sha256']==sha(V/'parent_release.json'),'ACTUAL155V2_TERMINAL')
 for n,k in [('operator_science_stdout.txt','stdout_sha256'),('operator_science_stderr.txt','stderr_sha256')]:pin(V/n,priorterm[k])
 need([c['case_id']for c in priorpacket['selected_cases']]==['case_024','case_015','case_008','case_003','case_000'],'FIRST_FIVE_UNCHANGED')
 byid={c['case_id']:c for c in catalog};reused_ids=[c['case_id']for c in priorpacket['selected_cases']];reused=[]
 for case in priorpacket['selected_cases']:
  cid=case['case_id'];need(case['canonical_key']==byid[cid]['canonical_key']and case['columns']==byid[cid]['columns'],'REUSED_CASE_IDENTICAL_SOURCE_LAW');rp=V/(cid+'_receipt.json');pin(rp);rec=read(rp);need(rec['status']=='PASS_ALL_SAVED_CASE_STAGES'and rec['canonical_key']==case['canonical_key'],'REUSED_COMPLETE_CASE_RECEIPT');reused.append({'case_id':cid,'canonical_key':case['canonical_key'],'actual155_case_receipt_path':str(rp),'actual155_case_receipt_sha256':sha(rp),'stages':[]})
 reusedby={c['case_id']:c for c in reused};seen=set()
 for entry in prior['completed']:
  key=(entry['case_id'],entry['active_columns']);need(key in stage_index and entry['case_id']in reusedby and key not in seen,'EXACT_REUSED_STAGE_SET');seen.add(key);original=stage_index[key]
  need(entry['source_stage_path']==original['path']and entry['source_stage_sha256']==original['sha256'],'REUSED_CERT_SAME_IMMUTABLE_SOURCE_STAGE')
  p=contained(V,entry['path']);pin(p,entry['sha256']);r=read(p);need(r['status']=='PASS_DIRECT_EXACT_INTEGER_LAPLACE_NUMERATOR_ZERO'and r['case_id']==entry['case_id']and r['active_columns']==entry['active_columns']and r['source_stage_path']==original['path']and r['source_stage_sha256']==original['sha256']and r['exact_integer_numerator_is_zero']is True and r['normalization']==1 and r['API_stage_ledger']['inflight_api']is None,'REUSED_ACTUAL_EXACT_CERTIFICATE')
  reusedby[entry['case_id']]['stages'].append({'active_columns':entry['active_columns'],'certificate_path':str(p),'certificate_sha256':sha(p),'source_stage_path':original['path'],'source_stage_sha256':original['sha256']})
 expected_reused={(cid,s['active_columns'])for cid in reused_ids for s in byid[cid]['stages']};need(seen==expected_reused and len(seen)==15,'COMPLETE_REUSED15_NO_REPLAY')
 new_ids=[c['case_id']for c in catalog if c['case_id']not in reused_ids];new_stage_count=sum(len(byid[cid]['stages'])for cid in new_ids)
 need(len(new_ids)==95 and new_stage_count==365 and set(new_ids).isdisjoint(reused_ids)and set(new_ids)|set(reused_ids)==set(byid),'EXACT_FIXED95_PLUS_REUSED5_PARTITION')
 totalweight=sum((F(c['pattern_weight'])for c in catalog),F(0));reusedweight=sum((F(byid[c]['pattern_weight'])for c in reused_ids),F(0));newweight=sum((F(byid[c]['pattern_weight'])for c in new_ids),F(0));need(newweight+reusedweight==totalweight==F(packet['selected_pattern_weight']),'FULL100_WEIGHT_PARTITION');need(sum(c['source_count']for c in catalog)==packet['selected_source_pair_count'],'ALL_SOURCE_PAIR_COUNT')
 caps=dict(priorpacket['caps']);caps['native_arithmetic_total']=365*caps['native_arithmetic_per_stage'];need(caps['native_arithmetic_total']==73000000,'FINITE_NEW_API_ATTEMPT_SCOPE')
 coverage={'status':'FIXED149100_EQUALS_REUSED155V2_5_PLUS_NEW95','all149_case_ids':[c['case_id']for c in catalog],'reused_case_ids':reused_ids,'NEW_audit_case_ids':new_ids,'all149_cases':100,'reused_cases':5,'NEW_audit_cases':95,'all149_stages':380,'reused_stages':15,'NEW_audit_stages':365,'full_pattern_weight':str(totalweight),'reused_pattern_weight':str(reusedweight),'NEW_audit_pattern_weight':str(newweight),'source_pair_count':packet['selected_source_pair_count'],'reused_actual_certificates':reused}
 write('coverage_ledger.json',coverage);write('case_packet.json',{'status':'PREPARED_UNRUN_FIXED95_SAVED_CASE_AUDIT_REUSE_ACTUAL155V2_FIVE','case_catalog':catalog,'NEW_audit_case_ids':new_ids,'reused_case_ids':reused_ids,'NEW_audit_cases':95,'NEW_audit_stages':365,'source149_pins':PINS149,'source155_pins':PINS155,'coverage_ledger_sha256':sha(P/'coverage_ledger.json'),'caps':caps,'workers':1,'nice_minimum':19,'counter_policy':'setup/construction/arithmetic ATTEMPTS capped73000000; every explicit query separately counted','no_case_selection_or_scope_expansion':True,'new_density_generations':0,'new_k_integrals':0,'legacy_native_integral_kernel_calls':0,'primary_baseline_replays':0})
 write('source_manifest.json',{'sources':[pins[p]for p in sorted(pins)]})
 source={n:sha(P/n)for n in ['prepare.py','api_ledger.py','verifier_core.py','run_prepared.py','PROOF.md','OPERATOR.md']}
 write('registration.json',{'status':'PREPARED_UNRUN_PENDING_REVIEW128_AND_MAIN_DUAL_GATE','source_sha256':source,'case_packet_sha256':sha(P/'case_packet.json'),'coverage_ledger_sha256':sha(P/'coverage_ledger.json'),'source_manifest_sha256':sha(P/'source_manifest.json'),'NEW_audit_cases':95,'NEW_audit_stages':365,'reused_actual_cases':5,'reused_actual_stages':15,'workers':1,'nice_minimum':19,'caps':caps,'time_limit':None,'automatic_restarts':0,'mathematical_core_is_actual155v2':PINS155['verifier_core.py'],'API_ledger_is_actual155v2':PINS155['api_ledger.py'],'scope_failure_is_not_mathematical_family_exclusion':True,'new_density_generations':0,'new_beta2_Cn_bound':False})
 write('prepare_receipt.json',{'status':'PREPARED_ONLY95_FIXED_SAVED_CASES_REUSE5_ACTUAL','registration_sha256':sha(P/'registration.json'),'case_packet_sha256':sha(P/'case_packet.json'),'coverage_ledger_sha256':sha(P/'coverage_ledger.json'),'new_verifier_or_FLINT_calls':0,'new_density_generations':0,'replayed155_cases':0})
 print(json.dumps({'status':'PREPARED_UNRUN','NEW_audit_cases':95,'NEW_audit_stages':365,'reused_cases':5,'setup_construction_arithmetic_attempt_budget':73000000,'registration_sha256':sha(P/'registration.json')}))
if __name__=='__main__':main()
