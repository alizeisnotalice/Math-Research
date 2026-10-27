"""Terra-only metadata freeze. Never imports a producer/evaluator or does numerical validation."""
from pathlib import Path
import hashlib,json,os,sys
sys.dont_write_bytecode=True
P=Path(__file__).resolve().parent
BASE=Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube6-sol-brascamp108/cube-dim-6-bounds')
TARGET=BASE/'audit/upper160_saved'
SRC=BASE/'upper/round160/normalized_joint_with_true_beta2'
PRIOR=Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube6-sol-cover090-audit/cube-dim-6-bounds/audit/upper128_129_preflight/round128_saved_review')
CODE=['prepare.py','check_saved.py','observe.py','PROOF_REVIEW.md','OPERATOR.md']
CERT='d573b8796ecfbdff133c7444228d89ae4bee99e2736c52988ce266713df8c519'
ACCEPT='b4743d20f78e1f4320aeeb602e84d803e6a2fac16c899633219012b9327f9381'
def read(p):return json.loads(Path(p).read_text())
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def need(b,m):
 if not b:raise RuntimeError(m)
def save(n,v):
 with(P/n).open('x')as f:json.dump(v,f,indent=2);f.write('\n');f.flush();os.fsync(f.fileno())
def main():
 need(P==TARGET,'TERRA_ROOT_ONLY')
 need(not any((P/n).exists()for n in ['inputs.json','source_manifest.json','registration.json','parent_release.json','started.json','observer_started.json','certificate.json']),'ONE_FRESH_PREPARE')
 d=SRC/'phaseB';need(sha(d/'candidate_certificate.json')==CERT,'ACTUAL160B_CERTIFICATE')
 need(sha(SRC/'inputs.json')=='b74e867ab92449f8604a5f0dc07a2ce9896163410f2b62444b943ed5a25250aa','ORIGINAL_INPUT_BYTES')
 need(sha(SRC/'quantization_certificate.json')=='737479ae78c4a841a77c9af23ec9b9f85a70d19d976ba0c6c65b6029d26159ad','ORIGINAL_CAP_WITNESS')
 cert=read(d/'candidate_certificate.json');need(len(cert['chunks'])==15 and len(cert['tail_records'])==16,'FINITE_SAVED_SCOPE')
 files=[SRC/n for n in ['inputs.json','quantization_certificate.json','STATIC_SOURCE_PINS.json','PACKAGE_MANIFEST.json','V2_PACKAGE_MANIFEST.json','common.py','run_phaseA.py','run_phaseA_v2.py','run_phaseB.py','prepare_for_terra.py','repair_prepare_v2.py','operator_prepare_v2.py','PROOF.md','PARTIAL_FAILURE_OBSERVATION.json','operator_prepare_failure_observation.json','V2_FAILURE_AND_DELTA.md','main_operator_phaseB_prepare_terminal.json','main_operator_phaseB_prepare_stderr.txt','main_operator_phaseB_prepare_v2_terminal.json','main_operator_phaseB_prepare_v2_stdout.txt','main_operator_phaseB_prepare_v2_stderr.txt']]
 files.extend(d/n for n in ['candidate_certificate.json','candidate.json','coefficients.json','integer_tables.json','jobs.json','selection.json','registration.json','parent_release.json','started.json','exit_receipt.json','operator_terminal_receipt.json','operator_science_stdout.txt','operator_science_stderr.txt'])
 files.extend(SRC/'phaseA'/n for n in ['result.json','exit_receipt.json','started.json','parent_release.json','operator_terminal_receipt.json'])
 for row in read(SRC/'STATIC_SOURCE_PINS.json')['sources']:
  path=Path(row['path']);need(sha(path)==row['sha256'],'STATIC_SOURCE_IDENTITY');files.append(path)
 for i,row in enumerate(cert['chunks']):
  need(row['job']==i and row['root_file']=='chunks/roots_'+str(i).zfill(2)+'.json.gz','FIXED_15_ROOT_FILES')
  path=d/row['root_file'];need(sha(path)==row['root_file_sha256'],'SAVED_ROOT_HASH');files.extend([path,d/'chunks'/('header_'+str(i).zfill(2)+'.json')])
 need({p.name for p in(d/'chunks').glob('*.gz')}=={Path(r['root_file']).name for r in cert['chunks']},'NO_UNLISTED_PAID_ROOT_FILE')
 release=read(d/'parent_release.json');audit=release['source156_saved_audit']
 for row in audit['evidence_files']:
  need(sha(row['path'])==row['sha256'],'156_ACCEPTED_EVIDENCE');files.append(Path(row['path']))
 acceptance=BASE/'upper/round156/pair_integral_without_full_density/main_acceptance.json'
 need(sha(acceptance)==ACCEPT,'ACTUAL156_MAIN_ACCEPTANCE');files.append(acceptance)
 original=read(acceptance)
 need(set(audit)==set(original) and all(audit[k]==v for k,v in original.items()if k!='evidence_files'),'ACCEPTANCE_ORIGINAL_FIELDS_PRESERVED')
 need(audit['evidence_files'][:len(original['evidence_files'])]==original['evidence_files']and len(audit['evidence_files'])==len(original['evidence_files'])+1,'ORIGINAL_EVIDENCE_PLUS_SINGLE_REAL_SEAL')
 need(audit['evidence_files'][-1]['sha256']=='18eca9a3a6b33c987952e42244221286efb143c038c54335a5bb02125e9713e4','ADDITIONAL_REAL156_SEAL')
 files.extend(SRC/n for n in ['operator_phaseB_prepare_receipt.json','operator_phaseB_prepare_stdout.txt','operator_phaseB_prepare_stderr.txt','main_operator_phaseB_prepare_stdout.txt','main_operator_phaseB_prepare_handle.json','main_operator_phaseB_prepare_reservation.json','main_operator_phaseB_prepare_v2_handle.json','main_operator_phaseB_prepare_v2_reservation.json'])
 files.append(PRIOR/'review_results.json')
 need(sha(PRIOR/'review_results.json')=='6e5ee7283d8e426926808ee91d33c4fe2673261b1540e0bcd0cc1795f62ba67e','ACTUAL128_INDEPENDENT_REVIEW')
 files.append(Path(sys.executable).resolve())
 rows={str(p):{'path':str(p),'sha256':sha(p),'bytes':p.stat().st_size}for p in files}
 limits={'workers':15,'fixed_prefix_jobs':15,'prefix_states':735471,'tail_records':16,'saved_root_inequalities':735487,'new_fifthroots':0,'new_Arb_exp_erfc_sqrt':0,'optimizer_native_RNG_density_old_baseline_calls':0,'arithmetic_bits':192,'max_saved_root_file_bytes':100000000,'max_saved_root_rows_per_job':49032,'nice_minimum':10,'time_or_CPU_limit':None,'automatic_retries':0,'configurations':1}
 save('inputs.json',{'producer_root':str(SRC),'producer_certificate_sha256':CERT,'source156_acceptance':str(acceptance),'source128_independent_review':str(PRIOR/'review_results.json'),'python_executable':str(Path(sys.executable).resolve()),'limits':limits})
 save('source_manifest.json',{'sources':[rows[k]for k in sorted(rows)]})
 save('registration.json',{'round':'160B_SAVED_INDEPENDENT','status':'PREPARED_UNRUN_PENDING_DUAL_GATE','code_sha256':{n:sha(P/n)for n in CODE},'inputs_sha256':sha(P/'inputs.json'),'source_manifest_sha256':sha(P/'source_manifest.json'),'limits':limits,'new_C6_bound':False})
 save('prepare_receipt.json',{'status':'METADATA_FREEZE_ONLY_NO_MATH','registration_sha256':sha(P/'registration.json'),'source_manifest_sha256':sha(P/'source_manifest.json'),'root_inequalities_or_math_calls':0})
 print(json.dumps({'status':'PREPARED_UNRUN','registration_sha256':sha(P/'registration.json')}))
if __name__=='__main__':main()
