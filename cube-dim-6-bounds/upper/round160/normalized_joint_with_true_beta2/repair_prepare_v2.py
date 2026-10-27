#!/usr/bin/env python3
"""Resume only missing metadata; preserve the failed v1 source and partial bytes."""
from pathlib import Path
from fractions import Fraction as F
from importlib.metadata import version as metadata_version_v2
from importlib.util import spec_from_file_location as spec_from_file_location_v2
from importlib.util import module_from_spec as module_from_spec_v2
import argparse,hashlib,json,os,shutil,sys
sys.dont_write_bytecode=True
sys.set_int_max_str_digits(0)
SOURCE_V2=Path(__file__).resolve().parent
def read_v2(p):return json.loads(Path(p).read_text())
def sha_v2(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def require_v2(ok,msg):
    if not ok:raise ArithmeticError(msg)
def new_v2(root,name,value):
    with(Path(root)/name).open('x')as f:json.dump(value,f,indent=2,allow_nan=False);f.write('\n');f.flush();os.fsync(f.fileno())
def version_metadata_v2():
    return {'numpy_version':metadata_version_v2('numpy'),'scipy_version':metadata_version_v2('scipy')}
def load_common_v2(root):
    spec=spec_from_file_location_v2('c160_preserved_common_v2',Path(root)/'common.py')
    mod=module_from_spec_v2(spec);spec.loader.exec_module(mod);return mod
def assert_snapshot_v2(root,source):
    observation=read_v2(Path(source)/'PARTIAL_FAILURE_OBSERVATION.json')
    require_v2(str(Path(root).resolve())==observation['target_root'],'ONE_FIXED_PARTIAL_TARGET')
    for r in observation['files']:
        require_v2(sha_v2(Path(root)/r['name'])==r['sha256'],'PRESERVE_PARTIAL:'+r['name'])
    require_v2(sha_v2(Path(root)/'PACKAGE_MANIFEST.json')==observation['original_package_sha256'],'PRESERVED_V1_PACKAGE')
    return observation
def assert_v2_files_v2(root):
    manifest=read_v2(Path(root)/'V2_PACKAGE_MANIFEST.json')
    for r in manifest['files']:
        require_v2(set(r)=={'name','bytes','sha256'}and Path(r['name']).name==r['name'],'V2_MANIFEST_NAME_SCHEMA')
        require_v2(sha_v2(Path(root)/r['name'])==r['sha256'],'FIXED_V2_SOURCE:'+r['name'])
def verify_existing_witness_v2(root):
    # Consumer validation of the already saved witness. No ceil and no new cap selection.
    q=read_v2(Path(root)/'quantization_certificate.json');inp=read_v2(Path(root)/'inputs.json')
    source=Path(q['source_result_path']);require_v2(sha_v2(source)==q['source_result_sha256']=='f9a0689da14fb34d48c3ffb65b6abaef06193906ff71cc3cab77420074bee21f','FULL_SOURCE156_RESULT')
    original=read_v2(source)['budgets']['finite_best']['valid_true_beta2_upper']
    require_v2(q['complete_source_upper_exact']==original and q['source_field_selector']==['budgets','finite_best','valid_true_beta2_upper'],'PRESERVED_FULL_UPPER_SOURCE_STRING')
    b=F(q['complete_source_upper_exact']);N=int(q['N']);D=int(q['D']);require_v2(D==1<<96 and b.denominator>0,'EXISTING96GRID_WITNESS')
    # Exactly three strict/weak integer cross-product inequalities:
    require_v2(0 < b.numerator*D,'WITNESS_1_POSITIVE_SOURCE')
    require_v2(b.numerator*D <= N*b.denominator,'WITNESS_2_UPPER_CONTAINMENT')
    require_v2(N*b.denominator < b.numerator*D+b.denominator,'WITNESS_3_LESS_THAN_ONE_GRID_STEP')
    cap=F(q['betaD_exact']);stored=F(inp['new_beta2_cap_exact']);model=F(inp['model_input']['beta_by_h'][5][1])
    require_v2(cap.numerator*D==N*cap.denominator and stored==cap==model,'UNCHANGED_SAVED_CAP_IDENTITY')
    return {'kind':'STRICT_CONSUMER_CHECK_OF_EXISTING_N_WITNESS','integer_cross_product_inequalities':3,'all_passed':True,'new_ceil_evaluations':0,'new_cap_selections':0,'old_M_H_or_scalar_evaluations':0,'input_bytes_unchanged':True,'quantization_bytes_unchanged':True}
def contract_v2(root):
    root=Path(root);base=load_common_v2(root).phase_requirements('phaseA');obs=read_v2(root/'PARTIAL_FAILURE_OBSERVATION.json')
    partial={r['name']:r['sha256']for r in obs['files']}
    base.update({'preparation_version':'V2_RESUME_MISSING_METADATA_ONLY','v2_package_manifest_sha256':sha_v2(root/'V2_PACKAGE_MANIFEST.json'),'v2_entrypoint_sha256':sha_v2(root/'run_phaseA_v2.py'),'partial_failure_observation_sha256':sha_v2(root/'PARTIAL_FAILURE_OBSERVATION.json'),'operator_failure_observation_sha256':partial['operator_prepare_failure_observation.json'],'v2_prepare_receipt_required':True,'engineering_scope_fixture_required':True,'original_failure_capture_incomplete_acknowledged':True,'old_partial_inputs_sha256':partial['inputs.json'],'old_quantization_certificate_sha256':partial['quantization_certificate.json'],'v1_raw_input_and_output_paths_unchanged':True})
    return base
def resume_v2(root):
    root=Path(root).resolve()
    require_v2(root!=SOURCE_V2 and root.parents[3]!=SOURCE_V2.parents[3]and(root.parents[3]/'.git').exists(),'INDEPENDENT_EXISTING_CLONE_ROOT')
    require_v2(root.name=='normalized_joint_with_true_beta2'and root.parent.name=='round160','FIXED_C160_ROOT')
    assert_v2_files_v2(SOURCE_V2);obs=assert_snapshot_v2(root,SOURCE_V2)
    require_v2((root/'phaseA').is_dir()and list((root/'phaseA').iterdir())==[],'OBSERVED_EMPTY_PHASE_A_STILL_EMPTY')
    require_v2(not(root/'prepare_v2').exists(),'NO_AUTOMATIC_V2_RESTART')
    # Copy only new versioned source names; never overwrite any original file.
    v2manifest=read_v2(SOURCE_V2/'V2_PACKAGE_MANIFEST.json')
    for n in [r['name']for r in v2manifest['files']]+['V2_PACKAGE_MANIFEST.json']:
        if(root/n).exists():require_v2(sha_v2(root/n)==sha_v2(SOURCE_V2/n),'EXISTING_V2_SOURCE_IDENTITY')
        else:shutil.copy2(SOURCE_V2/n,root/n)
    fixture=read_v2(root/'operator_v2_scope_fixture_receipt.json')
    require_v2(fixture['status']=='ENGINEERING_SCOPE_FIXTURE_PASS'and fixture['repair_prepare_v2_sha256']==sha_v2(root/'repair_prepare_v2.py')and fixture['operator_wrapper_v2_sha256']==sha_v2(root/'operator_prepare_v2.py'),'ACTUAL_ENGINEERING_SCOPE_FIXTURE_REQUIRED')
    d=root/'prepare_v2';d.mkdir(exist_ok=False)
    new_v2(d,'started.json',{'scope':'ONE_MISSING_METADATA_RECOVERY','pid':os.getpid(),'original_package_sha256':obs['original_package_sha256'],'v2_package_manifest_sha256':sha_v2(root/'V2_PACKAGE_MANIFEST.json'),'no_automatic_retry':True})
    try:
        witness=verify_existing_witness_v2(root)
        versions=version_metadata_v2()
        reg={'round':'160','phase':'phaseA','status':'TERRA_V2_METADATA_PREPARED_UNRUN','optimizer_invocations':1,'max_objective_gradient_calls':150,'max_optimizer_calls_before_final':149,'max_iterations':120,'ftol':1e-12,'workers':1,'nice_minimum':10,'time_or_CPU_limit':None,'automatic_restarts':0,'prefix_states':735471,'tail_rows':16,'features':49,'x0_numerical_floor':'1/16777216','candidate_denominator':'17592186044416','probability_renormalization':False,**versions,'package_manifest_sha256':sha_v2(root/'PACKAGE_MANIFEST.json'),'inputs_sha256':sha_v2(root/'inputs.json'),'quantization_certificate_sha256':sha_v2(root/'quantization_certificate.json'),'source156_saved_audit_required':True,'new_kernel_erfc_native_density_old_functional_calls':0,'new_M_upper':False,'new_C6_bound':False,'preparation_version':'V2_RESUME_MISSING_METADATA_ONLY','v2_package_manifest_sha256':sha_v2(root/'V2_PACKAGE_MANIFEST.json'),'original_failed_prepare_preserved':True,'original_failure_capture_complete':False}
        new_v2(root/'phaseA','registration.json',reg)
        requirements=contract_v2(root);new_v2(root/'phaseA','release_requirements.json',requirements)
        assert_snapshot_v2(root,root)
        receipt={'status':'V2_METADATA_RECOVERY_COMPLETE_SCIENCE_UNRUN','witness_consumer_validation':witness,'original_failure_capture_complete':False,'operator_failure_observation_sha256':requirements['operator_failure_observation_sha256'],'preserved_original_package_sha256':obs['original_package_sha256'],'preserved_inputs_sha256':sha_v2(root/'inputs.json'),'preserved_quantization_sha256':sha_v2(root/'quantization_certificate.json'),'registration_sha256':sha_v2(root/'phaseA/registration.json'),'release_requirements_sha256':sha_v2(root/'phaseA/release_requirements.json'),'scope_fixture_receipt_sha256':sha_v2(root/'operator_v2_scope_fixture_receipt.json'),'new_optimization_H_root_scalar_calls':0,'new_C6_bound':False}
        new_v2(d,'exit_receipt.json',receipt)
        print(json.dumps({'status':receipt['status'],'registration_sha256':receipt['registration_sha256'],'release_requirements_sha256':receipt['release_requirements_sha256'],'new_math_science_calls':0,'saved_witness_consumer_integer_inequalities':3}))
    except BaseException as err:
        new_v2(d,'exit_receipt.json',{'status':'V2_METADATA_RECOVERY_FAILED_STOP_NO_RETRY','error_type':type(err).__name__,'error':str(err),'new_optimization_H_root_scalar_calls':0});raise
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--root',required=True);args=ap.parse_args();resume_v2(Path(args.root))
