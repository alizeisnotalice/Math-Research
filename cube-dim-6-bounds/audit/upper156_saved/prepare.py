"""Terra-only source/receipt freeze for one independent C156 saved audit; no evaluator imports."""
from pathlib import Path
import json,hashlib,os,sys
sys.dont_write_bytecode=True
P=Path(__file__).resolve().parent
TARGET=Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube6-sol-brascamp108/cube-dim-6-bounds/audit/upper156_saved')
SRC=TARGET.parents[1]/'upper/round156/pair_integral_without_full_density'
AUDIT130=Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube6-sol-cover090-audit/cube-dim-6-bounds/audit/upper130_saved')
CODE=['prepare.py','check_saved.py','observe.py','PROOF_REVIEW.md','OPERATOR.md']
def read(p):return json.loads(Path(p).read_text())
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def need(x,m):
 if not x:raise RuntimeError(m)
def write(n,v):
 with(P/n).open('x')as f:json.dump(v,f,indent=2);f.write('\n');f.flush();os.fsync(f.fileno())
def main():
 need(P==TARGET,'TERRA_ONLY_ROOT');need(not any((P/n).exists()for n in ['inputs.json','registration.json','source_manifest.json','parent_release.json','started.json','observer_started.json','certificate.json']),'FRESH_ONE_PREPARE_ONLY')
 need(sha(SRC/'result.json')=='f9a0689da14fb34d48c3ffb65b6abaef06193906ff71cc3cab77420074bee21f','ACTUAL156_RESULT_PIN')
 result=read(SRC/'result.json');ex=read(SRC/'exit_receipt.json');term=read(SRC/'operator_terminal_receipt.json');started=read(SRC/'started.json')
 need(ex['result_sha256']==sha(SRC/'result.json')and ex['status']==result['status']=='COMPLETE_ONE_TRUE_JOINT_LAW_BOUND_FAMILY_DIAGNOSTIC'and term['exit_code']==0 and term['automatic_retries']==0 and term['pid']==started['pid'],'ACTUAL156_TERMINAL')
 need(sha(SRC/'parent_release.json')==term['release_sha256']==started['release_sha256'],'ORIGINAL_RELEASE')
 files=[SRC/n for n in ['result.json','inputs.json','source_manifest.json','registration.json','diagnose.py','PROOF.md','parent_release.json','started.json','exit_receipt.json','operator_terminal_receipt.json','operator_science_stdout.txt','operator_science_stderr.txt']]
 for name,key in [('operator_science_stdout.txt','stdout_sha256'),('operator_science_stderr.txt','stderr_sha256')]:need(sha(SRC/name)==term[key],'ORIGINAL_STDIO_HASH')
 original=read(SRC/'source_manifest.json')['sources']
 for row in original:need(sha(row['path'])==row['sha256']and Path(row['path']).stat().st_size==row['bytes'],'SOURCE156_UPSTREAM');files.append(Path(row['path']))
 headers=result['pair_files']+result['cauchy_root_certificate_files'];need(len(result['pair_files'])==21 and len(result['cauchy_root_certificate_files'])==1,'FIXED_ACTUAL_FILE_SCOPE')
 for row in headers:
  path=(SRC/row['file']).resolve();need(path.is_relative_to((SRC/'pairs').resolve())and sha(path)==row['sha256'],'PAIR_OR_ROOT_FILE_PIN');files.append(path)
 need({str(x.resolve())for x in(SRC/'pairs').glob('*.gz')}=={str((SRC/r['file']).resolve())for r in headers},'ALL_PAID_GZIP_FILES_LISTED')
 for n in ['manifest.json','source_pins.json','review_results.json']:
  need((AUDIT130/n).is_file(),'PRIOR130_CERTIFICATE_SOURCE');files.append(AUDIT130/n)
 need(read(AUDIT130/'review_results.json')['status']=='PASS_SAVED_BOX_COVERAGE_POSITIVE_SUMS_AND_ARB768_KERNELS','PRIOR130_SAVED_AUDIT_PASS')
 qualified130=[x['sha256']for x in read(AUDIT130/'source_pins.json')['files']if Path(x['path']).name=='results.json']
 current130=next(x['sha256']for x in original if Path(x['path']).name=='results.json')
 need(qualified130==[current130],'CURRENT130_RESULTS_ARE_SAME_INDEPENDENTLY_CERTIFIED_BYTES')
 prior_manifest={x['path']:x['sha256']for x in read(AUDIT130/'manifest.json')['files']}
 for n in ['source_pins.json','review_results.json']:need(sha(AUDIT130/n)==prior_manifest[n],'PRIOR130_CERTIFICATE_MANIFEST_LINK')
 runtime=Path('/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/venv/lib/python3.15/site-packages/flint')
 files.extend([Path(sys.executable).resolve(),runtime/'__init__.py',*[runtime/'.dylibs'/n for n in ['libflint.24.0.dylib','libgmp.10.dylib','libmpfr.6.dylib']]])
 # Pin installed native bindings too, without importing or testing their APIs.
 bindings=list(runtime.rglob('*.so'));need(bindings,'INSTALLED_FLINT_NATIVE_BINDINGS');files.extend(bindings)
 rows={str(x):{'path':str(x),'sha256':sha(x),'bytes':x.stat().st_size}for x in files}
 paths={Path(x['path']).name:x['path']for x in original}
 limits={'cones':21,'pairs_per_cone':496,'raw_pairs':10416,'max_caps_per_pair':35,'max_total_caps':364560,'actual_expected_caps':229314,'determinants':156240,'four_pairing_products':468720,'root_entries':147,'new_Cauchy_root_calls':0,'fresh_scalar_Arb_bits':768,'new_exp':3,'new_erfc':2,'new_erf':0,'new_sqrt_h':1,'new_sqrt_pi':1,'new_pi_constant':1,'workers':1,'nice_minimum':15,'time_or_CPU_limit':None}
 write('inputs.json',{'producer_root':str(SRC),'producer_result_sha256':sha(SRC/'result.json'),'upstream_paths':paths,'prior130_audit_manifest':str(AUDIT130/'manifest.json'),'python_executable':str(Path(sys.executable).resolve()),'flint_init_path':str((runtime/'__init__.py').resolve()),'limits':limits,'no_producer_import':True,'old_baseline_replays':0})
 write('source_manifest.json',{'sources':[rows[k]for k in sorted(rows)]})
 write('registration.json',{'round':'156_SAVED_INDEPENDENT','status':'PREPARED_UNRUN_PENDING_DUAL_GATE','source_sha256':{n:sha(P/n)for n in CODE},'inputs_sha256':sha(P/'inputs.json'),'source_manifest_sha256':sha(P/'source_manifest.json'),'limits':limits,'python_flint_version':'0.9.0','FLINT_version':'3.6.0','automatic_retries':0,'producer_or_old_baseline_calls':0,'new_true_beta2_candidate_only':True,'new_C6_bound':False})
 write('prepare_receipt.json',{'status':'SOURCE_FREEZE_ONLY_NO_MATH','registration_sha256':sha(P/'registration.json'),'inputs_sha256':sha(P/'inputs.json'),'source_manifest_sha256':sha(P/'source_manifest.json'),'checker_or_Arb_calls':0})
 print(json.dumps({'status':'PREPARED_UNRUN','registration_sha256':sha(P/'registration.json')}))
if __name__=='__main__':main()
