#!/usr/bin/env python3
"""Read-only byte/field audit of primary's one completed135 scalar.

No numerical consumer import, arithmetic expression evaluation, root check,
scalar rerun, seven-link rerun, H/R/gradient/log/exp/Gamma or optimization.
All copied numerical fields remain the primary's recorded evidence.
"""
from datetime import datetime,timezone
from pathlib import Path
import hashlib
import json
import os
import time

HERE=Path(__file__).resolve().parent
OWN=HERE.parent
MAIN=Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube-dim-5/cube-dim-5-bounds')
MAIN15=MAIN/'round015'
P132=Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube5-audit123-sol/cube-dim-5-bounds/upper_join132')
P134=Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube5-auditcurvature129-sol/cube-dim-5-bounds/upper_facepeak134')
LOCKS={
    'actual_inputs':(MAIN15/'join135_inputs.json','c6eb027efee798f140f6cb7362721b246de257a722b58cf65ba57ea1a6ad426c'),
    'actual_registration':(MAIN15/'join135_actual/registration.json','b1ca8c818d5cc0928141f83a72fc21f4f199792e0bc28f8b5d35d1fccf4b1bab'),
    'actual_results':(MAIN15/'join135_actual/results.json','8b0623db708b493cf7d196e5f8dc0f88724d401154a7c888da95987568502686'),
    'frozen_consumer':(OWN/'upper_join135/certify_one_scalar.py','80c162f2d1f90d6cc88140c700ce94fe83c776a0157ef279632b0fce17947512'),
    'frozen_plan':(OWN/'upper_join135/proposed_plan.json','98a200a9645f2532252af43ed03ac4871346aa95b2cc028aa3a835578c0df50d'),
    'frozen135':(OWN/'upper_join135/manifest.json','bc05ffe0e93315c764928f7e9b34dde450cccdc47777a5b223f8cffa41772d54'),
    'manifest134':(P134/'manifest.json','1d13b5eb6f7252818e2df203630fdead4a18efbc3e9915f91ff0241c36d645ea'),
    'condition134':(P134/'condition_results.json','1cf1aa32b44eb580226692bfe37ab40d2ef96deeda4edc155f16a2f0a48cb024'),
    'manifest132':(P132/'manifest.json','fbc30d041ad0d48480501633274d8883e4ee79868f4e19973a33eb48920745da'),
    'components132':(P132/'portable_components.json','2ec64713d80f2ecc791a7207a706d035b4c6723c24868cc44383c4d91ff9855b'),
    'IH_primary':(MAIN/'round014/primary_integrals_acceptance.json','14374097c6d164fcbaca56a4b6d003069e92e5e16ffa990a2ebb6cb470d2c219'),
}

def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def read(path):
    return json.loads(Path(path).read_bytes())

def field(obj,path):
    for key in path.split('.'):
        obj=obj[key]
    return obj

def write(name,value):
    with (HERE/name).open('x') as f:
        json.dump(value,f,sort_keys=True,indent=2);f.write('\n')

def main():
    assert not (HERE/'results.json').exists(),'Preserve completed source audit; no blind overwrite'
    began=time.monotonic()
    before={}
    for key,(path,expected) in LOCKS.items():
        assert sha(path)==expected,(key,'byte hash mismatch')
        before[key]={'path':str(path),'sha256':expected}
    authorization_path=MAIN15/'primary_join135_authorization.json'
    before['primary_authorization']={'path':str(authorization_path),'sha256':sha(authorization_path)}
    inp=read(LOCKS['actual_inputs'][0]);reg=read(LOCKS['actual_registration'][0]);result=read(LOCKS['actual_results'][0])
    plan=read(LOCKS['frozen_plan'][0]);auth=read(authorization_path)
    for name,row in plan['sources'].items():
        assert sha(row['path'])==row['sha256']
        before['plan/'+name]={'path':row['path'],'sha256':row['sha256']}
    auditreg={'task':'execution135_read_only_source_binding_final_check','pid':os.getpid(),
        'registered_utc':datetime.now(timezone.utc).isoformat(),'code_sha256':sha(__file__),
        'source_locks':before,'actual_scalar_recomputations':0,'actual_seven_chain_recomputations':0,
        'new_H_R_gradient_kernel_log_exp_Gamma_calls':0,'new_M_optimizations':0,
        'scope':'byte hashes and exact JSON field equality only; numerical formulas not re-evaluated'}
    write('registration.json',auditreg)
    assert auth['executor']=='primary uses readonly frozen author consumer; writes own main round015 output only'
    assert auth['code_path']==str(LOCKS['frozen_consumer'][0])
    for key,lock in [('code_sha256','frozen_consumer'),('plan_sha256','frozen_plan'),('inputs_sha256','actual_inputs')]:
        assert auth[key]==reg[key]==LOCKS[lock][1]
    assert auth['finite_actual_scalar_cases']==reg['case_count']==result['case_count']==1
    assert auth['root_bits']==reg['root_bits']==plan['root_bits']==192
    assert reg['pid']==result['pid']==51083 and reg['pool'] is None and result['pool'] is None
    assert inp['epsilon']==reg['epsilon']==result['epsilon']==plan['epsilon']=='1/512'
    assert inp['g']==reg['g']==plan['g']=='21/200000'
    assert reg['sources']==plan['sources'] and reg['qualifications']==inp['qualifications']
    assert result['registration_sha256']==LOCKS['actual_registration'][1]
    assert result['status']=='one_conditional_join_cap_scalar_complete_pending_primary'
    assert result['all_sources_unchanged_after'] is True
    assert result['global_I_g_Cn_claim'] is False and result['actual_J_R_or_I_lower_claim'] is False
    assert result['conditional_cost_upper_below_32'] is True
    # This reads the stored sign; it does not recompute the residual expression.
    numerator,denominator=result['conditional_cost_fourth_power_residual_32'].split('/')
    assert not numerator.startswith('-') and numerator!='0' and not denominator.startswith('-')
    matches=[]
    for quantity,input_key in [('Cbar','Cbar_upper'),('M','M_upper'),('IH','IH_upper')]:
        proof=inp['qualifications'][quantity]
        assert proof['epsilon']=='1/512' and proof['domain']=='complete_F'
        assert sha(proof['source_path'])==proof['source_sha256']
        source=read(proof['source_path'])
        source_value=field(source,proof['field'])
        assert inp[input_key]==source_value==result[input_key]
        matches.append({'quantity':quantity,'input_field':input_key,'source_path':proof['source_path'],
            'source_sha256':proof['source_sha256'],'source_field':proof['field'],
            'exact_recorded_value':source_value,'matches_input_registration_and_result':True,
            'qualification_remains_primary_and136_responsibility':True})
    c=read(LOCKS['condition134'][0]);m134=read(LOCKS['manifest134'][0])
    assert c['epsilon']==m134['epsilon']=='1/512'
    assert c['status']=='single_condition_and_new_A_R_bounds_certified'
    assert m134['condition_results_sha256']==LOCKS['condition134'][1]
    assert c['Cbar_upper_256dyadic']==m134['Cbar_upper_256dyadic']==inp['Cbar_upper']
    p132=read(LOCKS['components132'][0]);m132=read(LOCKS['manifest132'][0])
    assert p132['epsilon']=='1/512' and p132['q']=='5/4'
    assert p132['M_and_mu_reusable_without_H_or_T_optimization'] is True
    assert m132['components_sha256']==LOCKS['components132'][1]
    assert m132['M_analyticstationary_upper_cap']==inp['M_upper']
    assert p132['UH']==inp['IH_upper']
    ih=read(LOCKS['IH_primary'][0])
    assert ih['status']=='primary_accepts_complete_alpha1_I103_then_I112'
    assert ih['qualified_g']==inp['g'] and ih['full_face_all_positive_time'] is True
    assert ih['complete_H512_integral_upper']==inp['IH_upper']
    gproof=read(plan['sources']['qualified_g105']['path'])
    assert gproof['qualified_g']==inp['g']
    assert gproof['same029_general_L1_and_all_positive_common_scales'] is True
    stored=read(plan['sources']['stored076']['path'])['derived_constants']['pi_outward_384_dyadic']
    assert result['pi_interval']==[stored['lower'],stored['upper']]
    for row in before.values():
        assert sha(row['path'])==row['sha256'],('source changed during binding audit',row['path'])
    out={'task':'execution135','status':'PASS_primary_once_only_scalar_all_input_source_fields_bound',
        'audit_pid':os.getpid(),'audit_seconds':time.monotonic()-began,
        'actual_primary_scalar_pid':result['pid'],'actual_primary_seconds_as_recorded':result['seconds'],
        'primary_reported_exitcode':0,'exitcode_evidence':'trusted primary task message; not inferred from absence of a process',
        'actual_case_count':1,'actual_pool':None,'source_field_matches':matches,
        'selected_M':'132 analytic_stationary M_upper_portable_cap exact rational;not raw tighter M.upper or display decimal',
        'pi_source_field_match':True,'code_plan_authorization_input_registration_result_chain_match':True,
        'actual_registration_sha256':LOCKS['actual_registration'][1],
        'actual_results_sha256':LOCKS['actual_results'][1],
        'actual_inputs_sha256':LOCKS['actual_inputs'][1],
        'recorded_root_endpoints':[result['root_lower'],result['root_upper']],
        'recorded_complete_I_upper':result['conditional_complete_I_upper'],
        'recorded_residual_32_sign':'strictly positive as stored;expression not recalculated',
        'global_I_g_Cn_claim':False,'pending_gates':['independent136','primary029 and final analytic/source acceptance'],
        'actual_scalar_recomputations':0,'actual_seven_chain_recomputations':0,
        'new_H_R_gradient_kernel_log_exp_Gamma_calls':0,'new_M_optimizations':0,
        'all_sources_unchanged_after_audit':True,'frozen135_not_modified':True,
        'registration_sha256':sha(HERE/'registration.json'),'committed':False,'pushed':False}
    write('results.json',out)
    print(json.dumps({k:v for k,v in out.items() if k not in ['source_field_matches','recorded_root_endpoints','recorded_complete_I_upper']},sort_keys=True))

if __name__=='__main__':main()
