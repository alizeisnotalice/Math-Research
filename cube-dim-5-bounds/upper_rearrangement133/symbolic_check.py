#!/usr/bin/env python3
"""Finite new symbolic identities and artificial rational toys only."""
from datetime import datetime, timezone
from fractions import Fraction as F
from pathlib import Path
import hashlib
import json
import os
import time
import sympy as s
from scalar_interface import fourth_root_interval, mixed_increment_interval

HERE = Path(__file__).resolve().parent

def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def write(name, value):
    with (HERE/name).open('x') as f:
        json.dump(value, f, sort_keys=True, indent=2)
        f.write('\n')

def main():
    assert not (HERE/'symbolic_results.json').exists(), 'Preserve frozen finite result'
    reg = {'task':'upper_rearrangement133_new_symbolic_and_artificial_toys',
        'registered_utc':datetime.now(timezone.utc).isoformat(),'pid':os.getpid(),
        'code_sha256':sha(__file__),'interface_sha256':sha(HERE/'scalar_interface.py'),
        'symbolic_objects':['profile derivative','profile endpoint mass','beta/gamma product',
            'quartic rational decomposition','new norm fifth power','positive mixed integral M=D'],
        'finite_artificial_objects':{'stoploss_levels':[0,1,2,4,5],
            'fourth_root_inputs':['0','1','16','81/16','2'],
            'mixed_constant_case':{'epsilon':'1','Cbar':'2','P_interval':['1','1'],'M':'1','z_intervals':8}},
        'actual_physical_scalar_budget_cases':0,'new_H_R_gradient_kernel_numeric_log_Gamma_calls':0,
        'old_baseline_recomputations':0}
    write('symbolic_registration.json',reg)
    began = time.monotonic()
    C, eps, x = s.symbols('C eps x', positive=True)
    m = eps**4/(s.pi*C)
    D = C/(2*eps**2)
    Fprofile = s.sqrt(C*x/s.pi)-C*x/(2*eps**2)
    Rprofile = C/s.Integer(2)*((s.pi*C*x)**(-s.Rational(1,2))-eps**-2)
    assert s.simplify(s.diff(Fprofile,x)-Rprofile)==0
    assert s.simplify(Fprofile.subs(x,m)-eps**2/(2*s.pi))==0
    gamma_product = s.gammasimp(s.gamma(s.Rational(3,4))*s.gamma(s.Rational(9,4)))
    assert s.simplify(gamma_product-5*s.pi*s.sqrt(2)/16)==0
    z = s.symbols('z',positive=True)
    assert s.simplify((1+z**2)/(1+z**4)-s.Rational(1,2)*(1/(z*z+s.sqrt(2)*z+1)+1/(z*z-s.sqrt(2)*z+1)))==0
    q=s.Rational(5,4)
    Q=m*D**q*5*s.pi*s.sqrt(2)/16
    assert s.simplify(Q**4-s.Rational(625,524288)*C*eps**6)==0
    # For the abstract M=D comparator the positive integral is elementary.
    assert s.simplify(5*m*D*s.integrate(z*z*D**s.Rational(1,4),(z,0,1))-s.Rational(5,3)*m*D**q)==0
    actual = list(map(F,[4,2,1,1]))
    envelope = list(map(F,[5,5,1,1]))
    assert all(sum(sorted(actual,reverse=True)[:k]) <= sum(envelope[:k]) for k in range(1,5))
    stoploss=[]
    for level in map(F,[0,1,2,4,5]):
        a=sum(max(v-level,F(0)) for v in actual)/4
        b=sum(max(v-level,F(0)) for v in envelope)/4
        assert a<=b
        stoploss.append({'lambda':str(level),'actual_toy':str(a),'envelope_toy':str(b)})
    square_actual=sum(v*v for v in actual)/4
    square_envelope=sum(v*v for v in envelope)/4
    assert square_actual<=square_envelope
    # Phi(r)=-r is convex with Phi(0)=0 but not increasing: it reverses this.
    assert -sum(actual)/4 > -sum(envelope)/4
    for value in map(F,[0,1,16,F(81,16),2]):
        a,b=fourth_root_interval(value,24)
        assert a**4<=value<=b**4
    toy_lower,toy_upper=mixed_increment_interval(F(1),F(2),(F(1),F(1)),F(1),8,24)
    assert toy_lower<=F(5,6)<=toy_upper
    result={'task':'upper_rearrangement133','status':'PASS_new_symbolic_identities_and_artificial_toys',
        'pid':os.getpid(),'seconds':time.monotonic()-began,'registration_sha256':sha(HERE/'symbolic_registration.json'),
        'new_E_fifth_power_identity':'625/524288*C_epsilon*epsilon^6',
        'Gamma3over4_Gamma9over4':'5*pi*sqrt(2)/16','stoploss_toy':stoploss,
        'square_convex_toy':{'actual':str(square_actual),'envelope':str(square_envelope)},
        'negative_linear_convex_counterexample_verified':True,
        'artificial_mixed_constant_case':{'exact':'5/6','enclosure':[str(toy_lower),str(toy_upper)],'is_physical_budget':False},
        'actual_physical_scalar_budget_cases':0,'new_H_R_gradient_kernel_numeric_log_Gamma_calls':0,
        'old_baseline_recomputations':0,'new_global_I_g_or_C5':False}
    write('symbolic_results.json',result)
    print(json.dumps(result,sort_keys=True))

if __name__=='__main__':
    main()
