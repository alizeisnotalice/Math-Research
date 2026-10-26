#!/usr/bin/env python3
"""Read-only final bindings of primary's completed136 run. No math-case rerun."""
from datetime import datetime, timezone
from pathlib import Path
import hashlib
import json
import os
import subprocess
import time

HERE=Path(__file__).resolve().parent
MAIN=Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube-dim-5/cube-dim-5-bounds')
ROUND=MAIN/'round015'

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def read(p):return json.loads(Path(p).read_bytes())
def select(data,keys):
    for key in keys:data=data[key]
    return data

def main():
    output=HERE/'execution136_checks.json'
    assert not output.exists(),'Preserve completed read-only binding receipt'
    start=time.monotonic()
    named={
        'audit136_inputs':(ROUND/'audit136_inputs.json','98668f14db23ac5f3d824e2e6570752dffcf09e618656fc39b653193abcc7c33'),
        'audit136_registration':(ROUND/'audit136_actual/registration.json','c0db2e449c1f5f77b1cd62b168fb88eb4a3d3be0e8bb67a066cf8d65edb2f3d6'),
        'audit136_results':(ROUND/'audit136_actual/results.json','bcc9caa56f076a2f0a20256304b8d48b8d4596b189bffa26573dba703f53b92c'),
        'join135_results':(ROUND/'join135_actual/results.json','8b0623db708b493cf7d196e5f8dc0f88724d401154a7c888da95987568502686'),
        'own_audit_code':(HERE/'audit_candidate.py','90f1f77221ace67f32adac70f338f9bbdabebc40753004c3b45c71b4fedd2b07'),
        'own_plan':(HERE/'proposed_plan.json','fc077f35031d8d2b021b5e423d9b9f3e80d1cc6b2852cb051c49a22f1e10b3db'),
        'own_prep_manifest':(HERE/'manifest.json','f614285f1e561dbed528557a38173a220f9aaa2864c0e6813b27ff2c4480e083'),
    }
    before={str(p):sha(p) for p,_ in named.values()}
    for p,expected in named.values():assert before[str(p)]==expected,str(p)
    inp=read(named['audit136_inputs'][0]);reg=read(named['audit136_registration'][0]);res=read(named['audit136_results'][0])
    preg=read(ROUND/'join135_actual/registration.json');pres=read(ROUND/'join135_actual/results.json')
    pinput=read(ROUND/'join135_inputs.json')
    auth=read(ROUND/'primary_audit136_authorization.json');pauth=read(ROUND/'primary_join135_authorization.json')
    assert reg['code_sha256']==auth['code_sha256']==named['own_audit_code'][1]
    assert reg['plan_sha256']==auth['plan_sha256']==named['own_plan'][1]
    assert reg['inputs_sha256']==auth['input_sha256']==named['audit136_inputs'][1]
    assert reg['pid']==res['pid']==51896
    assert res['registration_sha256']==named['audit136_registration'][1]
    assert reg['sources']==inp['sources'] and reg['M_pointer']==inp['M_pointer']
    assert reg['maximum_physical_cases']==1 and reg['workers']==1 and reg['precision_bits']==384
    source_data={}
    for name,row in inp['sources'].items():
        p=Path(row['path']);before[str(p)]=sha(p)
        assert before[str(p)]==row['sha256'],name
        source_data[name]=read(p)
    for row in preg['sources'].values():
        p=Path(row['path']);before[str(p)]=sha(p)
        assert before[str(p)]==row['sha256'],str(p)
    pcode=Path(pauth['code_path'])
    pplan=pcode.parent/'proposed_plan.json'
    assert sha(pcode)==preg['code_sha256']==pauth['code_sha256']
    assert sha(pplan)==preg['plan_sha256']==pauth['plan_sha256']
    assert sha(ROUND/'join135_inputs.json')==preg['inputs_sha256']==pauth['inputs_sha256']
    assert sha(ROUND/'join135_actual/registration.json')==pres['registration_sha256']
    assert preg['pid']==pres['pid']==51083
    assert preg['case_count']==pres['case_count']==1
    assert pinput['qualifications']==preg['qualifications']
    for quantity,key in [('Cbar','Cbar_upper'),('M','M_upper'),('IH','IH_upper')]:
        qual=pinput['qualifications'][quantity]
        assert qual['epsilon']==inp['epsilon']=='1/512' and qual['domain']=='complete_F'
        p=Path(qual['source_path']);before[str(p)]=sha(p)
        assert before[str(p)]==qual['source_sha256']
        value=select(read(p),qual['field'].split('.'))
        assert value==pinput[key]==pres[key],quantity
    assert inp['g']==pinput['g']==source_data['g105']['qualified_g']=='21/200000'
    assert inp['eta']==source_data['condition134']['eta']=='1/10'
    assert inp['epsilon']==source_data['condition134']['epsilon']=='1/512'
    assert inp['M_qualification']['epsilon']=='1/512'
    assert inp['M_qualification']['domain']=='complete_F'
    assert inp['M_qualification']['bound']=='H_epsilon<=M'
    assert select(source_data['M132'],inp['M_pointer'])==res['M_consumed']==pres['M_upper']
    assert source_data['integrals014']['complete_H512_integral_upper']==res['IH_consumed']==pres['IH_upper']
    assert source_data['condition134']['Cbar_upper_256dyadic']==res['Cbar_134_dyadic_consumed']==pres['Cbar_upper']
    pp=source_data['stored076']['derived_constants']['pi_outward_384_dyadic']
    assert pres['pi_interval']==[pp['lower'],pp['upper']]
    assert res['producer_complete_I_upper']==pres['conditional_complete_I_upper']
    assert res['producer_exact32_gate']==pres['conditional_cost_fourth_power_residual_32']
    assert res['producer_cap_endpoints_pass_exact_fourth_powers'] is True
    assert res['cost_strictly_below32'] is True and pres['conditional_cost_upper_below_32'] is True
    assert res['actual_J_R_or_I_lower_claim'] is False and pres['actual_J_R_or_I_lower_claim'] is False
    assert res['new_Cn_claim_without_primary_acceptance'] is False
    assert res['all_source_hashes_unchanged'] is True and pres['all_sources_unchanged_after'] is True
    assert res['short_rational_chain_complete_I_upper']==inp['short_rational_chain']['I_upper']
    assert res['short_rational_chain_cost_upper']==inp['short_rational_chain']['cost_upper']
    # Validate unchanged preparation artifacts, without rewriting the earlier manifest.
    prep=read(HERE/'manifest.json')
    for name,row in prep['files'].items():assert sha(HERE/name)==row['sha256'],name
    for p in [pcode,pplan,ROUND/'join135_inputs.json',ROUND/'join135_actual/registration.json',
        ROUND/'primary_audit136_authorization.json',ROUND/'primary_join135_authorization.json']:
        before[str(p)]=sha(p)
    process=[]
    for pid in (51896,51083):
        done=subprocess.run(['ps','-p',str(pid),'-o','pid=,stat=,command='],text=True,capture_output=True)
        assert done.returncode==1 and not done.stdout.strip(),f'PID {pid} present; inspect before qualification'
        process.append({'pid':pid,'ps_returncode':done.returncode,'stdout':done.stdout,'stderr':done.stderr,'present':False})
    after={path:sha(path) for path in before}
    assert before==after,'A read-only source changed during final binding review'
    result={
        'task':'audit_upper136_readonly_execution_binding',
        'status':'PASS_source_code_input_qualification_bindings_no_case_rerun',
        'review_pid':os.getpid(),'review_utc':datetime.now(timezone.utc).isoformat(),
        'review_seconds':time.monotonic()-start,
        'review_code_sha256':sha(__file__),
        'mathematical_code_author':'independent Sol136',
        'actual_complete_candidate_executor':'primary',
        'actual_complete_candidate_pid':51896,
        'actual_complete_candidate_exitcode':{'value':0,'evidence':'trusted primary task message; not inferred from JSON status'},
        'actual_complete_candidate_seconds_recorded':res['seconds'],
        'actual_producer135_executor':'primary','actual_producer135_pid':51083,
        'actual_candidates_reexecuted_in_this_review':0,
        'seven_rational_links_reexecuted_in_this_review':0,
        'math_or_root_functions_called_in_this_review':0,
        'all_current_source_hashes_match_runtime_locks':True,
        'binding_review_before_hashes':before,'binding_review_after_hashes':after,
        'historical_runtime_before_after_statement':'Readonly136 code checked source hashes on entry and exit; frozen result records true. This review independently rechecks current bytes twice; it does not invent a historical separate snapshot.',
        'runtime_pid_checks':process,
        'independent_complete_I_upper_copied':res['independent_complete_I_upper'],
        'independent_cost_interval_copied':res['independent_cost_interval'],
        'independent_cost_display_only_copied':res['cost_display_only'],
        'short_rational_I_cap_copied':res['short_rational_chain_complete_I_upper'],
        'short_rational_cost_cap_copied':res['short_rational_chain_cost_upper'],
        'cost_interval_lower_is_Cn_lower':False,
        'shared_dependencies':['076 saved C3/C1/pi','057 P3 normalization via132','132 complete-face M','111/112/121 complete UH','g105 global cover','029 finite-source and all-scale transfer'],
        'formal_Cn_ledger_acceptance_owner':'primary',
        'new_H_R_gradient_kernel_log_exp_Gamma_calls':0,
        'old_baseline_recomputations':0,'shared_MATLAB_calls':0,
        'agents_dispatched':0,'lower_work':0,'committed':False,'pushed':False,
    }
    with output.open('x') as stream:json.dump(result,stream,ensure_ascii=False,indent=2,sort_keys=True);stream.write('\n')
    print(json.dumps({'status':result['status'],'source_files_checked':len(before),'actual_candidates_reexecuted':0,'pid_checks':process},ensure_ascii=False))

if __name__=='__main__':main()
