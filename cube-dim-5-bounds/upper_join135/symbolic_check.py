#!/usr/bin/env python3
"""Finite new135 symbolic and artificial rational checks; no133 toy rerun."""
from datetime import datetime,timezone
from fractions import Fraction as F
from pathlib import Path
import hashlib
import json
import os
import time
import sympy as s
from certify_one_scalar import join_cap_interval

HERE=Path(__file__).resolve().parent

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def write(name,value):
    with (HERE/name).open('x') as f:
        json.dump(value,f,sort_keys=True,indent=2);f.write('\n')

def main():
    assert not (HERE/'symbolic_results.json').exists(),'Preserve completed new135 toy'
    reg={'task':'upper_join135_new_symbolic_and_artificial_checks',
        'pid':os.getpid(),'registered_utc':datetime.now(timezone.utc).isoformat(),
        'symbolic_code_sha256':sha(__file__),'consumer_sha256':sha(HERE/'certify_one_scalar.py'),
        'finite_symbolic_objects':['probability normalization','fourth moment3/7',
            'mean affine argument','Jensen coefficient cancellation','root concavity derivative'],
        'finite_artificial_toys':['two-point Jensen D<M','two-point Jensen M<D','constant M=D exact cap'],
        'actual_physical_scalar_cases':0,'new_H_R_gradient_kernel_numeric_log_exp_Gamma_calls':0,
        'old133_toy_recomputations':0,'old_baseline_recomputations':0}
    write('symbolic_registration.json',reg);began=time.monotonic()
    z=s.symbols('z',nonnegative=True)
    D,M,C,eps=s.symbols('D M C eps',positive=True)
    assert s.integrate(3*z*z,(z,0,1))==1
    assert s.integrate(3*z*z*z**4,(z,0,1))==s.Rational(3,7)
    assert s.simplify(s.integrate(3*z*z*(D+(M-D)*z**4),(z,0,1))-(4*D+3*M)/7)==0
    x=s.symbols('x',positive=True)
    assert s.simplify(s.diff(x**s.Rational(1,4),x,2)+s.Rational(3,16)*x**(-s.Rational(7,4)))==0
    m=eps**4/(s.pi*C);d=C/(2*eps**2)
    assert s.simplify(s.Rational(5,3)*m*d-5*eps**2/(6*s.pi))==0
    mean_roots=(F(1)+F(2))/2;mean_powers=(F(1)+F(16))/2
    assert mean_roots**4<mean_powers
    assert ((F(2)+F(1))/2)**4<(F(16)+F(1))/2
    # Artificial P=1 in place of pi, epsilon1,Cbar32 =>D16,M16.
    # Jensen is equality, root=2, cap=(5/6)*2=5/3.
    toy=join_cap_interval(F(1),F(32),F(16),(F(1),F(1)),24)
    assert toy['D']==16 and toy['root_argument']==16
    assert toy['root_lower']==toy['root_upper']==2
    assert toy['cap_lower']==toy['cap_upper']==F(5,3)
    out={'task':'upper_join135','status':'PASS_new_Jensen_identities_and_artificial_toys',
        'pid':os.getpid(),'seconds':time.monotonic()-began,
        'registration_sha256':sha(HERE/'symbolic_registration.json'),
        'normalized_weight':'3*z^2 dz','E_z4':'3/7','E_argument':'(4D+3M)/7',
        'cap_formula':'5epsilon^2/(6pi)*((4D+3M)/7)^(1/4)',
        'two_point_root_mean_fourth_power':str(mean_roots**4),
        'two_point_argument_mean':str(mean_powers),
        'artificial_constant_case':{k:str(v) for k,v in toy.items()},
        'artificial_constant_case_is_physical_budget':False,
        'actual_physical_scalar_cases':0,'new_H_R_gradient_kernel_numeric_log_exp_Gamma_calls':0,
        'old133_toy_recomputations':0,'old_baseline_recomputations':0,'new_I_g_Cn_claim':False}
    write('symbolic_results.json',out);print(json.dumps(out,sort_keys=True))

if __name__=='__main__':main()
