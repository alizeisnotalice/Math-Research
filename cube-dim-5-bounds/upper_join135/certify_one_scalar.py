#!/usr/bin/env python3
"""One prepared join-cap scalar; real physical budget NOT YET RUN.

Consumes qualified Cbar/M/IH and saved pi endpoints. No H/R/kernel,
log/exp/Gamma/Arb calls or H-sup optimization. All root endpoints are saved.
"""
from datetime import datetime, timezone
from fractions import Fraction as F
from math import isqrt
from pathlib import Path
import argparse
import hashlib
import json
import os
import time

HERE=Path(__file__).resolve().parent

def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def fourth_root_interval(value,bits=192):
    value=F(value)
    assert value>=0 and type(bits) is int and bits>=1
    den=1<<bits
    n=value.numerator*den**4//value.denominator
    k=isqrt(isqrt(n))
    lower=F(k,den)
    assert lower**4<=value<F(k+1,den)**4
    upper=lower if lower**4==value else F(k+1,den)
    return lower,upper

def join_cap_interval(epsilon,Cbar,M,pi_interval,bits=192):
    """Enclose Jcap at chosen conservative Cbar,M; not an actual J lower."""
    epsilon,Cbar,M=map(F,(epsilon,Cbar,M))
    pil,piu=map(F,pi_interval)
    assert epsilon>0 and Cbar>0 and M>0 and 0<pil<=piu
    D=Cbar/(2*epsilon**2)
    root_argument=(4*D+3*M)/7
    root_lower,root_upper=fourth_root_interval(root_argument,bits)
    factor=F(5,6)*epsilon**2
    cap_lower=factor*root_lower/piu
    cap_upper=factor*root_upper/pil
    return {'D':D,'root_argument':root_argument,'root_lower':root_lower,
        'root_upper':root_upper,'cap_lower':cap_lower,'cap_upper':cap_upper}

def write(path,obj):
    with Path(path).open('x') as f:
        json.dump(obj,f,sort_keys=True,indent=2)
        f.write('\n')

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--reviewed-plan-sha',required=True)
    ap.add_argument('--inputs',type=Path,required=True)
    ap.add_argument('--reviewed-inputs-sha',required=True)
    ap.add_argument('--output',type=Path,required=True)
    args=ap.parse_args()
    planpath=HERE/'proposed_plan.json'
    assert sha(planpath)==args.reviewed_plan_sha
    assert sha(args.inputs)==args.reviewed_inputs_sha
    plan=json.loads(planpath.read_bytes())
    assert plan['maximum_actual_scalar_cases']==1 and sha(__file__)==plan['consumer_sha256']
    for row in plan['sources'].values():
        assert sha(row['path'])==row['sha256']
    inputs=json.loads(args.inputs.read_bytes())
    eps=F(inputs['epsilon']);g=F(inputs['g'])
    assert str(eps)==plan['epsilon']=='1/512' and str(g)==plan['g']=='21/200000'
    Cbar,M,IH=map(F,(inputs['Cbar_upper'],inputs['M_upper'],inputs['IH_upper']))
    assert Cbar>0 and M>0 and IH>=0
    for quantity in ('Cbar','M','IH'):
        row=inputs['qualifications'][quantity]
        assert row['epsilon']==str(eps) and row['domain']=='complete_F'
        assert sha(row['source_path'])==row['source_sha256']
    assert inputs['qualifications']['Cbar']['bound']=='p_t<=Cbar*t^-3 for0<t<=epsilon'
    assert inputs['qualifications']['M']['bound']=='H_epsilon<=M'
    assert inputs['qualifications']['IH']['bound']=='integral_F H_epsilon^(5/4)<=IH'
    args.output.mkdir(exist_ok=False)
    reg={'task':'upper_join135_one_join_cap_scalar','pid':os.getpid(),
        'registered_utc':datetime.now(timezone.utc).isoformat(),'code_sha256':sha(__file__),
        'plan_sha256':args.reviewed_plan_sha,'inputs_sha256':args.reviewed_inputs_sha,
        'case_count':1,'epsilon':str(eps),'g':str(g),'root_bits':plan['root_bits'],
        'sources':plan['sources'],'qualifications':inputs['qualifications'],
        'workers':1,'pool':None,'wall_deadline':None,
        'new_H_R_gradient_kernel_log_exp_Gamma_calls':0,'old_baseline_recomputations':0,
        'primary_must_validate_source_theorems_not_only_hashes':True}
    write(args.output/'registration.json',reg)
    began=time.monotonic()
    pi_source=json.loads(Path(plan['sources']['stored076']['path']).read_bytes())['derived_constants']['pi_outward_384_dyadic']
    pi_pair=(F(pi_source['lower']),F(pi_source['upper']))
    values=join_cap_interval(eps,Cbar,M,pi_pair,plan['root_bits'])
    Iupper=IH+values['cap_upper']
    residual=F(32)**4*g**5-F(3125)*Iupper**4
    result={'task':'upper_join135','status':'one_conditional_join_cap_scalar_complete_pending_primary',
        'pid':os.getpid(),'epsilon':str(eps),'Cbar_upper':str(Cbar),'M_upper':str(M),'IH_upper':str(IH),
        'pi_interval':list(map(str,pi_pair)),**{k:str(v) for k,v in values.items()},
        'root_endpoint_power_checks':{'lower_fourth_power':str(values['root_lower']**4),
            'upper_fourth_power':str(values['root_upper']**4),'direction':'lower^4<=argument<=upper^4'},
        'conditional_complete_I_upper':str(Iupper),'conditional_cost_fourth_power_residual_32':str(residual),
        'conditional_cost_upper_below_32':residual>0,
        'actual_J_R_or_I_lower_claim':False,'global_I_g_Cn_claim':False,
        'requires_primary_acceptance':['Cbar complete shorttime peak qualification',
            'same-epsilon complete-face M and IH qualifications','same-model g105 and029 transfer'],
        'case_count':1,'new_H_R_gradient_kernel_log_exp_Gamma_calls':0,
        'old_baseline_recomputations':0,'pool':None,'registration_sha256':sha(args.output/'registration.json')}
    for row in plan['sources'].values():
        assert sha(row['path'])==row['sha256']
    for row in inputs['qualifications'].values():
        assert sha(row['source_path'])==row['source_sha256']
    result['all_sources_unchanged_after']=True
    result['seconds']=time.monotonic()-began
    write(args.output/'results.json',result)
    print(json.dumps(result,sort_keys=True))

if __name__=='__main__':
    main()
