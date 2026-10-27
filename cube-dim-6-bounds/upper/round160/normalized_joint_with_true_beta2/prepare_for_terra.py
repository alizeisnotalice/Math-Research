#!/usr/bin/env python3
"""Terra preparation only. No optimizer, law enumeration, Joint import or root."""
from pathlib import Path
from fractions import Fraction as F
import argparse,hashlib,importlib.metadata,json,shutil,sys
sys.dont_write_bytecode=True
sys.set_int_max_str_digits(0)
SOURCE=Path(__file__).resolve().parent
from common import RESULT156_SHA,FIELD156,CODE_FILES,read,sha,req,new,rows,coeffs
def phaseA(root):
    req(not any((root/n).exists()for n in ['inputs.json','phaseA/registration.json','phaseA/started.json','phaseA/result.json']),'FRESH_PHASE_A_PREPARATION')
    root.mkdir(parents=True,exist_ok=True)
    package=read(SOURCE/'PACKAGE_MANIFEST.json')
    for r in package['files']:req(sha(SOURCE/r['name'])==r['sha256'],'FIXED_SOURCE_PACKAGE')
    for n in CODE_FILES+['PACKAGE_MANIFEST.json']:
        if(root/n).exists():req(sha(root/n)==sha(SOURCE/n),'EXISTING_SOURCE_BYTES')
        else:shutil.copy2(SOURCE/n,root/n)
    pins=read(root/'STATIC_SOURCE_PINS.json')
    for r in pins['sources']:req(sha(r['path'])==r['sha256'],'FROZEN_EXTERNAL_SOURCE')
    s={r['role']:Path(r['path'])for r in pins['sources']};old=read(s['128_input']);r156=read(s['156_result']);candidate=read(s['128_candidate'])
    req(sha(s['156_result'])==RESULT156_SHA,'FULL156RESULT')
    b=F(r156['budgets']['finite_best']['valid_true_beta2_upper']);D=1<<96;N=(b.numerator*D+b.denominator-1)//b.denominator;betaD=F(N,D)
    req(r156['h']=='1/128'and read(s['156_inputs'])['h']=='1/128'and b==F(r156['budgets']['finite_best']['budget_expression_interval'][1]),'FULL156_H_AND_UPPER_ENDPOINT_IDENTITY')
    req(0<b<=betaD<b+F(1,D)and betaD<F(old['beta_by_h'][5][1]),'UPWARD96_AND_STRICT_CAP_IMPROVEMENT')
    req(old['h_values'][5]=='1/128'and all(old['Phi_table'][n][2]==n*(n-1)//2 for n in range(33)),'TRUE_H1_128_RANK2_IDENTITY')
    req((N-1)*b.denominator<D*b.numerator<=N*b.denominator,'INTEGER_CEILING_CERTIFICATE')
    new(root,'quantization_certificate.json',{'source_result_path':str(s['156_result']),'source_result_sha256':RESULT156_SHA,'source_field_selector':FIELD156,'complete_source_upper_exact':str(b),'D':str(D),'N':str(N),'betaD_exact':str(betaD),'strict_integer_lower_ceiling_inequality':True,'weak_integer_upper_ceiling_inequality':True,'BU_le_betaD_lt_BU_plus_2minus96':True,'source156_saved_audit_required_before_science':True,'from_display_decimal':False})
    m={k:old[k]for k in ['model','h_values','beta0','p_intervals','Phi_table','threshold_M_power6_exact','prefix_Kmax','prefix_state_count','tail_K_values']};m['beta_by_h']=[row.copy()for row in old['beta_by_h']];m['beta_by_h'][5][1]=str(betaD)
    coeffs(m,candidate)
    cert128=read(s['128_certificate'])
    new(root,'inputs.json',{'round':'160','model_input':m,'warm_start_49':candidate,'changed_flat_index':32,'changed_table_entry':[5,1],'h':'1/128','rank':2,'new_beta2_cap_exact':str(betaD),'threshold_gstar_power6_exact':str(F(8,9)**6*24576**5/63500**6),'expected_state_hash':cert128['state_bytes_sha256'],'saved128_certificate_path':str(s['128_certificate']),'qualification':'genuine upper cap; normalized cost<=1, not optimistic bL; no science yet'})
    a=root/'phaseA';a.mkdir(exist_ok=False)
    new(a,'registration.json',{'round':'160','phase':'phaseA','status':'TERRA_PREPARED_UNRUN','optimizer_invocations':1,'max_objective_gradient_calls':150,'max_optimizer_calls_before_final':149,'max_iterations':120,'ftol':1e-12,'workers':1,'nice_minimum':10,'time_or_CPU_limit':None,'automatic_restarts':0,'prefix_states':735471,'tail_rows':16,'features':49,'x0_numerical_floor':'1/16777216','candidate_denominator':str(1<<44),'probability_renormalization':False,'numpy_version':importlib.metadata.version('numpy'),'scipy_version':importlib.metadata.version('scipy'),'package_manifest_sha256':sha(root/'PACKAGE_MANIFEST.json'),'inputs_sha256':sha(root/'inputs.json'),'quantization_certificate_sha256':sha(root/'quantization_certificate.json'),'source156_saved_audit_required':True,'new_kernel_erfc_native_density_old_functional_calls':0,'new_M_upper':False,'new_C6_bound':False})
    # Importing the copied helpers does not invoke a mathematical preflight.
    import importlib.util
    spec=importlib.util.spec_from_file_location('c160_target_common_A',root/'common.py');mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
    new(a,'release_requirements.json',mod.phase_requirements('phaseA'))
    print(json.dumps({'status':'PHASE_A_PREPARED_UNRUN','root':str(root),'source156_saved_audit_required_before_science':True,'new_objective_H_root_calls':0,'phaseB_prepared':False}))
def phaseB(root,selection_path):
    req((root/'phaseA/result.json').exists(),'ACTUAL_PHASE_A_REQUIRED')
    req(not(root/'phaseB').exists(),'FRESH_PHASE_B_PREPARATION')
    import importlib.util
    spec=importlib.util.spec_from_file_location('c160_target_common_B',root/'common.py');mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
    inp,s=mod.model_preflight();result=read(root/'phaseA/result.json');receipt=read(root/'phaseA/exit_receipt.json');sel=read(selection_path)
    req(receipt['status']=='COMPLETE_ONE_GENUINE_CAP_NUMERICAL_SELECTOR'and receipt['result_sha256']==sha(root/'phaseA/result.json'),'COMPLETE_BOUND_PHASE_A_RESULT')
    req(result['status']=='COMPLETE_ONE_GENUINE_CAP_NUMERICAL_SELECTOR'and result['new_M_upper']is False and result['H_power5_below_saved128_float_diagnostic_only']is True,'EVIDENCED_NEW_CANDIDATE_ONLY')
    req(sel.get('round')=='160'and sel.get('main_selected')is True and sel.get('phaseB_prepare_authorized')is True and isinstance(sel.get('evidence'),str)and len(sel['evidence'])>0,'EXPLICIT_PRIMARY_SELECTION_REQUIRED')
    req(sel.get('phaseA_result_sha256')==sha(root/'phaseA/result.json')and sel.get('candidate_canonical_sha256')==hashlib.sha256(mod.canon(result['candidate'])).hexdigest(),'SELECTION_BINDS_ACTUAL_NEW_CANDIDATE')
    mod.check_source156_acceptance(read(root/'phaseA/parent_release.json'))
    x,caps,a=coeffs(inp['model_input'],result['candidate']);d=root/'phaseB';d.mkdir(exist_ok=False)
    new(d,'selection.json',sel);new(d,'candidate.json',result['candidate']);new(d,'coefficients.json',{'raw_a_quantized':False,'a_exact':list(map(str,a)),'cap_vector_exact':list(map(str,caps)),'sum_a_cap_exact':'1','canonical_a_sha256':hashlib.sha256(mod.canon(list(map(str,a)))).hexdigest()})
    jobs=[];start=0;q,r=divmod(735471,15)
    for j in range(15):
        end=start+q+(j<r);jobs.append({'job':j,'start':start,'end_exclusive':end,'states':end-start});start=end
    req(start==735471,'FINITE_COMPLETE_JOBS');new(d,'jobs.json',jobs)
    new(d,'registration.json',{'round':'160','phase':'phaseB','status':'TERRA_PREPARED_UNRUN','new_candidates':1,'workers':15,'prefix_states':735471,'tail_rows':16,'logical_root_requests':735487,'arithmetic_bits':192,'raw_a_quantized':False,'root_cache_per_process_cap':65536,'nice_minimum':10,'time_or_CPU_limit':None,'automatic_restarts':0,'new_optimizer_kernel_erfc_native_density_old_functional_calls':0,'phaseA_result_sha256':sha(root/'phaseA/result.json'),'candidate_sha256':sha(d/'candidate.json'),'coefficients_sha256':sha(d/'coefficients.json'),'selection_sha256':sha(d/'selection.json'),'jobs_sha256':sha(d/'jobs.json'),'package_manifest_sha256':sha(root/'PACKAGE_MANIFEST.json'),'inputs_sha256':sha(root/'inputs.json'),'quantization_certificate_sha256':sha(root/'quantization_certificate.json'),'new_C6_bound':False})
    new(d,'release_requirements.json',mod.phase_requirements('phaseB'))
    print(json.dumps({'status':'PHASE_B_PREPARED_UNRUN','root':str(root),'prefix_states':735471,'workers':15,'new_H_root_calls':0,'awaiting_separate_dual_gate':True}))
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--root',required=True);ap.add_argument('--phase',choices=['A','B'],required=True);ap.add_argument('--selection');args=ap.parse_args();root=Path(args.root).resolve()
    req(root!=SOURCE and root.parents[3]!=SOURCE.parents[3]and(root.parents[3]/'.git').exists(),'INDEPENDENT_EXISTING_CLONE_ROOT')
    req(root.name=='normalized_joint_with_true_beta2'and root.parent.name=='round160'and root.parents[1].name=='upper'and root.parents[2].name=='cube-dim-6-bounds','EXACT_C160_ROOT')
    if args.phase=='A':req(args.selection is None,'NO_SELECTION_IN_PHASE_A');phaseA(root)
    else:req(args.selection is not None,'PHASE_B_NEEDS_PRIMARY_SELECTION');phaseB(root,Path(args.selection).resolve())
if __name__=='__main__':main()
