#!/usr/bin/env python3
"""Fraction checks for the p=3 binary reduction and its endpoint derivative."""
from fractions import Fraction as Q
from pathlib import Path
import json

ROOT=Path(__file__).resolve().parent

def direct(q,delta):
    a=1-delta/q; b=(1-q*a)/(1-q)
    ef3=q*a**3+(1-q)*b**3
    return (ef3-1)/(4*delta*delta)

def formula(q,delta):
    return Q(3,1)/(4*q*(1-q))+delta*(2*q-1)/(4*q*q*(1-q)*(1-q))

def endpoint(q): return (2-q)/(4*q*(1-q)*(1-q))

def main():
    checks=[]
    for q in (Q(1,7),Q(1,4),Q(2,5),Q(3,7),Q(3,5),Q(5,6)):
        for delta in (q/Q(4),q/Q(2),q*Q(3,4)):
            left=direct(q,delta);right=formula(q,delta);assert left==right
            checks.append({'q':str(q),'delta':str(delta),'R':str(left)})
    # Endpoint derivative has sign -sign(q^2-3q+1) on (0,1).
    signs=[]
    for q in (Q(1,4),Q(3,8),Q(19,50),Q(2,5),Q(1,2),Q(3,4)):
        poly=q*q-3*q+1
        # A positive common denominator leaves exactly this sign in -R'(q).
        signs.append({'q':str(q),'q2_minus_3q_plus_1':str(poly),'endpoint_R':str(endpoint(q))})
    assert Q(19,50)**2-3*Q(19,50)+1>0 and Q(2,5)**2-3*Q(2,5)+1<0
    # The root is (3-sqrt(5))/2; brackets 19/50 < root < 2/5 follow from signs.
    out={'scope':'exact rational checks only','affine_formula':'3/[4q(1-q)]+delta(2q-1)/[4q^2(1-q)^2]',
         'endpoint':'(2-q)/[4q(1-q)^2] for q<1/2 and delta=q (a downarrow 0)',
         'checks':checks,'derivative_sign_samples':signs,
         'root_bracket':'19/50 < (3-sqrt(5))/2 < 2/5'}
    (ROOT/'p3_exact_audit.json').write_text(json.dumps(out,indent=2))
    print('p3_fraction_identity_checks=%d derivative_bracket PASS' % len(checks))
if __name__=='__main__':main()
