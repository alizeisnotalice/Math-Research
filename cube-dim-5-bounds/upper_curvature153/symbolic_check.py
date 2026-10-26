#!/usr/bin/env python3
"""New153 moment identities/artificial ratios only; no150 symbols rerun."""
from datetime import datetime,timezone
from fractions import Fraction as F
from pathlib import Path
import hashlib,json,os,time
import sympy as s
from penalty import curvature_caps

HERE=Path(__file__).resolve().parent
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def write(n,x):
    with (HERE/n).open('x') as f:json.dump(x,f,sort_keys=True,indent=2);f.write('\n')
def main():
    assert not (HERE/'symbolic_results.json').exists(),'Preserve finite153 result'
    reg={'task':'upper_curvature153_new_symbolic_and_artificial_ratios','pid':os.getpid(),'utc':datetime.now(timezone.utc).isoformat(),
        'code_sha256':sha(__file__),'interface_sha256':sha(HERE/'penalty.py'),
        'finite_new_identities':['global second-inverse-time moment2/3','finite doubled-cutoff moment7/10','tail doubled-cutoff moment1/6',
            'mixture coefficient8/15','tail ratio1/16','same-input129-versus153 gap','q128 below256 threshold27/32'],
        'artificial_ratio_examples':3,'invalid_input_examples':4,
        'actual_point_q_evaluations':0,'actual_H_gradient_log_corner_calls':0,'old150_symbol_recomputations':0}
    write('symbolic_registration.json',reg);began=time.monotonic()
    e,t=s.symbols('e t',positive=True);r,A=s.symbols('r A',real=True)
    mass=lambda a,b:s.integrate(t**-5,(t,a,b))
    moment=lambda a,b:s.integrate(t**-7,(t,a,b))/mass(a,b)
    assert s.simplify(moment(e,s.oo)-s.Rational(2,3)/e**2)==0
    assert s.simplify(moment(e,2*e)-s.Rational(7,10)/e**2)==0
    assert s.simplify(moment(2*e,s.oo)-s.Rational(1,6)/e**2)==0
    assert s.simplify((1-r)*s.Rational(7,10)+r*s.Rational(1,6)-(s.Rational(7,10)-s.Rational(8,15)*r))==0
    assert s.simplify(mass(2*e,s.oo)/mass(e,s.oo)-s.Rational(1,16))==0
    assert s.simplify((s.Rational(7,10)-s.Rational(8,15)/16)-s.Rational(2,3))==0
    assert s.simplify(A*(1-s.Rational(3,4)*r)-A*(s.Rational(7,10)-s.Rational(8,15)*r)-A*(18-13*r)/60)==0
    assert s.solve(s.Eq(s.Rational(7,10)-s.Rational(8,15)*r,s.Rational(1,4)),r)==[s.Rational(27,32)]
    a=curvature_caps(F(1,64),0,1);b=curvature_caps(F(1,128),F(9,10),1);c=curvature_caps(F(1,128),0,1,F(1,100))
    assert a['new_q_upper']==F(512,3)
    assert b['new_q_upper']==F(5632,25)<256
    assert c['new_q_upper']==F(1,100)<c['A']/6
    invalid=[(0,0,1,None),(F(1,64),2,1,None),(F(1,64),1,1,None),(F(1,64),0,1,0)]
    for args in invalid:
        try:curvature_caps(*args)
        except ValueError:pass
        else:raise AssertionError('Invalid source/parameter shape accepted')
    out={'task':'upper_curvature153','status':'PASS_new153_identities_and_artificial_ratios','pid':os.getpid(),'seconds':time.monotonic()-began,
        'registration_sha256':sha(HERE/'symbolic_registration.json'),
        'constants':{'global':'2/3','first_band':'7/10','second_band':'1/6','mixture_subtraction':'8/15','true_r_min':'1/16'},
        'artificial_examples':{'eps64_L0_q':str(a['new_q_upper']),'eps128_rLower9over10_q':str(b['new_q_upper']),'qualified_old_min_interface_toy':str(c['new_q_upper'])},
        'examples_are_not_actual_H_or_q_evidence':True,'invalid_inputs_rejected':len(invalid),
        'actual_point_q_evaluations':0,'actual_H_gradient_log_corner_calls':0,'old150_symbol_recomputations':0,'new_global_g_or_Cn_claim':False}
    write('symbolic_results.json',out);print(json.dumps(out,sort_keys=True))
if __name__=='__main__':main()
