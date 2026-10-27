"""C160B actual audit-package metadata/byte gate only; no root record decode, math or producer imports."""
from pathlib import Path
import hashlib,json,os,time,collections
P=Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube6-sol-brascamp108/cube-dim-6-bounds/audit/upper160_saved');A=Path(__file__).resolve().parent
os.nice(max(0,19-os.getpriority(os.PRIO_PROCESS,0)))
def read(p):return json.loads(Path(p).read_text())
def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb')as f:
  for b in iter(lambda:f.read(1048576),b''):h.update(b)
 return h.hexdigest()
counts=collections.Counter()
def ck(v,n):
 counts[n]+=1
 if not v:raise RuntimeError(n)
t=time.monotonic();expected={'inputs.json':'006397a805460d4f52575b7656cd0eb5718291016d3b23c08665c27efd61e252','source_manifest.json':'184bb765e108b9de88053d01de6f1ce65e1ff3d2b379c5558595f2b6c7b891e9','registration.json':'97b0252223d250b9f647321415504f2e180ecf69664a59c492564c613e70e4bb','prepare_receipt.json':'7e1352e774d3b6273dd10cca25e589f0625e8afd7dd9de9fbccd2ec2484c373c'}
for n,h in expected.items():ck(sha(P/n)==h,'actual_main_package_hash')
ck(sha(P/'code_manifest.json')=='6343adc441221cc329c2668bdae178c84da665a7eac6179cf732f2a88588f7b0','already_reviewed_code_manifest')
cm=read(P/'code_manifest.json');reg=read(P/'registration.json');inp=read(P/'inputs.json');manifest=read(P/'source_manifest.json');receipt=read(P/'prepare_receipt.json')
for row in cm['files']:ck(sha(P/row['path'])==row['sha256']and(P/row['path']).stat().st_size==row['bytes'],'reviewed_code_copy_unchanged')
for n,h in reg['code_sha256'].items():ck(sha(P/n)==h,'registered_five_code_members')
ck(reg['inputs_sha256']==expected['inputs.json']and reg['source_manifest_sha256']==expected['source_manifest.json']and reg['limits']==inp['limits'],'actual_registration_input_source_scope')
ck(reg['status']=='PREPARED_UNRUN_PENDING_DUAL_GATE'and reg['round']=='160B_SAVED_INDEPENDENT'and receipt['status']=='METADATA_FREEZE_ONLY_NO_MATH'and receipt['root_inequalities_or_math_calls']==0 and receipt['registration_sha256']==expected['registration.json'],'actual_metadata_only_prepare_receipt')
for row in manifest['sources']:ck(sha(row['path'])==row['sha256']and Path(row['path']).stat().st_size==row['bytes'],'all_frozen_source_bytes')
pinned={row['path']:row['sha256']for row in manifest['sources']};ck(len(pinned)==len(manifest['sources'])==104,'unique104_source_paths')
limit=inp['limits'];ck(limit=={'workers':15,'fixed_prefix_jobs':15,'prefix_states':735471,'tail_records':16,'saved_root_inequalities':735487,'new_fifthroots':0,'new_Arb_exp_erfc_sqrt':0,'optimizer_native_RNG_density_old_baseline_calls':0,'arithmetic_bits':192,'max_saved_root_file_bytes':100000000,'max_saved_root_rows_per_job':49032,'nice_minimum':10,'time_or_CPU_limit':None,'automatic_retries':0,'configurations':1},'exact_fixed_finite_resources')
ck(inp['python_executable']=='/opt/homebrew/Cellar/python@3.15/3.15.0_1/Frameworks/Python.framework/Versions/3.15/bin/python3.15'and inp['python_executable']in pinned,'actual_registered_python_byte_identity')
src=Path(inp['producer_root']);d=src/'phaseB';cert=read(d/'candidate_certificate.json');jobs=read(d/'jobs.json')
ck(sha(d/'candidate_certificate.json')==inp['producer_certificate_sha256']=='d573b8796ecfbdff133c7444228d89ae4bee99e2736c52988ce266713df8c519','exact_one_actual160B_candidate')
ck(len(jobs)==len(cert['chunks'])==15 and jobs[0]['start']==0 and jobs[-1]['end_exclusive']==735471,'fixed_complete15_job_partition')
rootnames=set();total_declared_states=0
for i,(job,row)in enumerate(zip(jobs,cert['chunks'])):
 header=read(d/'chunks'/('header_'+str(i).zfill(2)+'.json'));f=d/row['root_file'];rootnames.add(f.name)
 ck(header==row and all(header[k]==v for k,v in job.items())and job['job']==i and job['states']==job['end_exclusive']-job['start']and 0<job['states']<=49032,'exact_original_header_job_identity')
 ck(i==0 or jobs[i-1]['end_exclusive']==job['start'],'no_job_gap_or_overlap')
 ck(row['root_file']=='chunks/roots_'+str(i).zfill(2)+'.json.gz'and sha(f)==row['root_file_sha256']==pinned[str(f)]and f.stat().st_size<100000000,'saved_root_file_byte_identity_no_decode')
 total_declared_states+=job['states']
ck(rootnames=={p.name for p in(d/'chunks').glob('*.gz')}and total_declared_states==735471,'all_and_only15_saved_prefix_files')
ck(len(cert['tail_records'])==16 and[row['K']for row in cert['tail_records']]==list(range(17,33))and cert['prefix_states']==735471 and cert['logical_root_requests']==735487 and cert['arithmetic_bits']==192,'original16_tail_and192_metadata')
ck(cert['status']=='COMPLETE_ONE_NEW_NORMALIZED_FULL_H_CERTIFICATE','actual_producer_terminal_scope')
release=read(d/'parent_release.json');accept=read(inp['source156_acceptance']);aug=release['source156_saved_audit'];audit156=read(P.parent/'upper156_saved/certificate.json');roles={row['role']:Path(row['path'])for row in read(src/'STATIC_SOURCE_PINS.json')['sources']};r156=read(roles['156_result']);q=read(src/'quantization_certificate.json')
ck(sha(inp['source156_acceptance'])=='b4743d20f78e1f4320aeeb602e84d803e6a2fac16c899633219012b9327f9381'and accept['main_accepted']is True and accept['independent_saved_review_passed']is True,'actual_original156_main_acceptance')
ck(set(aug)==set(accept)and all(aug[k]==v for k,v in accept.items()if k!='evidence_files')and aug['evidence_files'][:-1]==accept['evidence_files']and aug['evidence_files'][-1]['sha256']=='18eca9a3a6b33c987952e42244221286efb143c038c54335a5bb02125e9713e4','original_acceptance_plus_one_seal_not_replaced')
for row in aug['evidence_files']:ck(pinned[row['path']]==row['sha256']==sha(row['path']),'all156_actual_audit_evidence_pinned')
ck(audit156['producer_result_sha256']==accept['result_sha256']==sha(roles['156_result'])==q['source_result_sha256']==cert['source156_result_sha256'],'actual156_result_value_source_identity')
ck(audit156['certified_true_beta2_upper']==accept['accepted_true_beta2_upper_exact']==r156['budgets']['finite_best']['valid_true_beta2_upper']==q['complete_source_upper_exact'],'exact_saved_cap_strings_no_comparison_or_ceiling')
ck(sha(inp['source128_independent_review'])=='6e5ee7283d8e426926808ee91d33c4fe2673261b1540e0bcd0cc1795f62ba67e'and read(inp['source128_independent_review'])['status']=='PASS_COMPLETE_INDEPENDENT_SAVED_JOINT_CANDIDATE','accepted128_dependency_not_replayed')
g1=read(src/'operator_phaseB_prepare_receipt.json');g2=read(src/'main_operator_phaseB_prepare_terminal.json');ck(g1['pid']==1971 and g1['exit_code']==1 and g2['child_pid']==2329 and g2['exit_code']==1 and release['known_guard_failures_preserved']==2,'two_distinct_B_guard_failures_not_merged')
ck(all(str(src/n)in pinned for n in('PARTIAL_FAILURE_OBSERVATION.json','operator_prepare_failure_observation.json','V2_FAILURE_AND_DELTA.md')),'distinct_earlier_partial_A_prepare_evidence_preserved')
op=read(P/'operator_prepare_terminal_receipt.json');ck(op['exit_code']==0 and op['pid']==7426,'actual_unique_metadata_prepare_exit')
for k,n in [('stdout_sha256','operator_prepare_stdout.txt'),('stderr_sha256','operator_prepare_stderr.txt')]:ck(op[k]==sha(P/n),'actual_prepare_observer_log_chain')
for name in('started.json','observer_started.json','workers','audit_model.json','certificate.json','exit_receipt.json'):ck(not(P/name).exists(),'saved_math_unstarted')
for n,h in expected.items():ck(sha(P/n)==h,'prepared_hashes_unchanged_after_gate')
out={'status':'PASS_CONCRETE_160B_SAVED_PACKAGE_READY_FOR_MAIN_UNIQUE_OBSERVER_RELEASE','scope':'metadata counts/strings/header partition and file-byte identities only; no root decode or arithmetic verification','actual_hashes':expected,'code_manifest_sha256':sha(P/'code_manifest.json'),'checks':dict(counts),'actual_prepare_PID':7426,'source_rows':104,'fixed_prefix_jobs':15,'declared_prefix_states':735471,'tail_records':16,'logical_saved_root_checks_for_Terra':735487,'root_record_decodes_or_math_checks_by_Sol':0,'source_imports_fixture_prepare_root_Arb_optimizer_or_baseline_runs':0,'currentbest_or_C6_update_authorized_by_this_gate':False,'execution_entry':'registered Python -B observe.py only, then its original checker and15 owned workers','new_main_exact_release_required':True,'metadata_gate_wall_seconds':time.monotonic()-t}
(A/'CONCRETE_METADATA_GATE.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out))
