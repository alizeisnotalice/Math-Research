#!/usr/bin/env python3
"""Prepared bounded scalar consumer; actual physical cases NOT YET RUN.

Requires primary-reviewed plan/input hashes and a fresh output directory.
Consumes stored076 constants only. Never evaluates H/R, optimizes H sup,
imports a numerical kernel, or calls log/exp/Gamma/Arb.
"""
from datetime import datetime, timezone
from fractions import Fraction as F
from pathlib import Path
import argparse
import hashlib
import json
import os
import time
from scalar_interface import norm_fifth_power, mixed_increment_interval

HERE = Path(__file__).resolve().parent

def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def write(p, value):
    with Path(p).open('x') as f:
        json.dump(value, f, sort_keys=True, indent=2)
        f.write('\n')

def fifth_root_interval(value, bits):
    value=F(value)
    assert value>=0
    den=1<<bits
    n=value.numerator*den**5//value.denominator
    lo,hi=0,1<<((n.bit_length()+4)//5)
    while hi-lo>1:
        mid=(lo+hi)//2
        if mid**5<=n:
            lo=mid
        else:
            hi=mid
    a=F(lo,den)
    assert a**5<=value<F(lo+1,den)**5
    return a,a if a**5==value else F(lo+1,den)

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--reviewed-plan-sha',required=True)
    ap.add_argument('--inputs',type=Path,required=True)
    ap.add_argument('--reviewed-inputs-sha',required=True)
    ap.add_argument('--output',type=Path,required=True)
    args=ap.parse_args()
    planpath=HERE/'proposed_scalar_plan.json'
    assert sha(planpath)==args.reviewed_plan_sha
    assert sha(args.inputs)==args.reviewed_inputs_sha
    plan=json.loads(planpath.read_bytes())
    assert sha(__file__)==plan['scalar_consumer_sha256']
    assert sha(HERE/'scalar_interface.py')==plan['scalar_interface_sha256']
    locked=plan['sources']
    for source in locked.values():
        assert sha(source['path'])==source['sha256']
    inputs=json.loads(args.inputs.read_bytes())
    eps=F(inputs['epsilon'])
    assert str(eps)==plan['epsilon']=='1/512'
    assert inputs['norm'] is True
    mixed=inputs.get('mixed')
    assert mixed is None or isinstance(mixed,dict)
    args.output.mkdir(exist_ok=False)
    reg={'task':'upper_rearrangement133_finite_scalar_budget','pid':os.getpid(),
        'registered_utc':datetime.now(timezone.utc).isoformat(),'plan_sha256':args.reviewed_plan_sha,
        'inputs_sha256':args.reviewed_inputs_sha,'code_sha256':sha(__file__),
        'norm_cases':1,'mixed_cases':int(mixed is not None),'max_total_cases':2,
        'z_intervals_per_mixed_case':plan['z_intervals'],'root_bits':plan['root_bits'],
        'workers':1,'pool':None,'wall_deadline':None,'sources':locked,
        'new_H_R_gradient_kernel_log_exp_Gamma_calls':0,'old_baseline_recomputations':0}
    write(args.output/'registration.json',reg)
    began=time.monotonic()
    old=json.loads(Path(locked['stored076']['path']).read_bytes())['derived_constants']
    Cbar=F(old['C3']['upper'])+F(old['C1']['upper'])*eps**2
    pi_interval=(F(old['pi_outward_384_dyadic']['lower']),F(old['pi_outward_384_dyadic']['upper']))
    fifth=norm_fifth_power(eps,Cbar)
    norm=fifth_root_interval(fifth,plan['root_bits'])
    result={'task':'upper_rearrangement133','status':'finite_scalar_budget_complete_pending_primary_qualification',
        'pid':os.getpid(),'epsilon':str(eps),'Cbar':str(Cbar),
        'new_norm_fifth_power_upper':str(fifth),
        'new_norm_certificate_interval':list(map(str,norm)),
        'new_norm_upper_display':float(norm[1]),'actual_R_norm_lower_claim':False,
        'mixed':None,'registration_sha256':sha(args.output/'registration.json'),
        'new_H_R_gradient_kernel_log_exp_Gamma_calls':0,'old_baseline_recomputations':0,
        'global_I_g_or_C5_claim':False}
    write(args.output/'norm_complete.json',result)
    if mixed is not None:
        assert mixed['epsilon']==str(eps) and mixed['domain']=='complete_F'
        assert mixed['component']=='H_epsilon' and F(mixed['M_upper'])>=0
        assert sha(mixed['source_path'])==mixed['source_sha256']
        lower,upper=mixed_increment_interval(eps,Cbar,pi_interval,F(mixed['M_upper']),plan['z_intervals'],plan['root_bits'])
        result['mixed']={'M_upper':mixed['M_upper'],'M_source_path':mixed['source_path'],
            'M_source_sha256':mixed['source_sha256'],'comparator_increment_interval':[str(lower),str(upper)],
            'actual_R_or_I_lower_claim':False,'upper_display':float(upper),
            'requires_M_global_sup_qualification':True}
        if 'IH_upper' in mixed:
            assert F(mixed['IH_upper'])>=0
            assert sha(mixed['IH_source_path'])==mixed['IH_source_sha256']
            result['mixed']['conditional_complete_I_upper']=str(F(mixed['IH_upper'])+upper)
            result['mixed']['IH_upper']=mixed['IH_upper']
            result['mixed']['IH_source_path']=mixed['IH_source_path']
            result['mixed']['IH_source_sha256']=mixed['IH_source_sha256']
            result['mixed']['requires_same_epsilon_complete_IH_qualification']=True
    for source in locked.values():
        assert sha(source['path'])==source['sha256']
    if mixed is not None:
        assert sha(mixed['source_path'])==mixed['source_sha256']
        if 'IH_upper' in mixed:
            assert sha(mixed['IH_source_path'])==mixed['IH_source_sha256']
    result['seconds']=time.monotonic()-began
    result['all_sources_unchanged_after']=True
    write(args.output/'results.json',result)
    print(json.dumps(result,sort_keys=True))

if __name__=='__main__':
    main()
