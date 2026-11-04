#!/usr/bin/env python3
"""Exact endpoint windows for homogeneous norm-infinity singularities.

For t^(1/n)<=1/8 and 0<tau<=0.9, (1+2t^(1/n))^tau<2.  The matching
lower estimate is therefore (1/2)t^(tau/n), while a background B/a=R
contributes at most 2^n R t/(1-t).  This script certifies the exact dyadic
comparison after (not in place of) the analytic geometric estimates.
"""
from fractions import Fraction as F
from math import lcm
from pathlib import Path
import json

ROOT=Path(__file__).resolve().parent
TAUS=(F(1,10),F(1,4),F(2,5),F(1,2),F(3,4),F(9,10))
RATIOS=(F(4),F(16),F(64),F(130))
Q_AUDIT=(F(1),F(2))

def dyadic(exponent):
    assert exponent.denominator==1 and exponent>=0
    return F(1,2**exponent.numerator)

def endpoint_window(n,tau,R):
    theta=tau/F(n); assert 0<theta<1
    # L makes both t^(1/n) and t^theta exact dyadic quantities.
    L=lcm(n,theta.denominator)
    for k in range(1,2001):
        exponent=L*k
        if exponent<3*n: continue  # t^(1/n)<=1/8
        t=dyadic(F(exponent))
        lower=F(1,2)*dyadic(F(exponent)*theta)
        upper=(2**n)*R*t/(1-t)
        if lower>upper:
            return dict(n=n,tau=str(tau),theta=str(theta),R=str(R),L=L,k0=k,
                        t0=str(t),lower=str(lower),background_upper=str(upper),
                        strict=True,geometry_condition='t^(1/n)<=1/8')
    raise AssertionError('window not found')

def second_order_audit(n,tau):
    theta=tau/F(n); twice_theta=2*theta; tests=[]
    for q in Q_AUDIT:
        polynomial=tau*tau-(n+2*q)*tau+n*q
        valid=(tau<1 and tau<F(n,2) and polynomial>0)
        if valid: assert twice_theta<1
        tests.append({'q':str(q),'quadratic_positive':polynomial>0,'quadratic_value':str(polynomial),
                      'valid_second_order':valid})
    return {'n':n,'tau':str(tau),'theta':str(theta),'twice_theta':str(twice_theta),
            'tau_lt_1':tau<1,'tau_lt_n_over_2':tau<F(n,2),
            'twice_theta_lt_1':twice_theta<1,'second_order_tests':tests}

def main():
    rows=[];audit=[]
    for n in (1,2,3):
        for tau in TAUS:
            # First-order endpoint window needs theta<1, not twice_theta<1.
            for R in RATIOS: rows.append(endpoint_window(n,tau,R))
            audit.append(second_order_audit(n,tau))
    out={'scope':'exact verification of stated lower/upper endpoint comparison',
         'analytic_connection':('For 0<t<=t0, the quotient of (1/2)t^theta and '
             '2^n R t/(1-t) increases as t decreases because theta<1. Hence the '
             'strict t0 inequality holds throughout (0,t0].'),
         'count':len(rows),'rows':rows,'second_order_audit':audit,
         'limits':'Exact algebra only after the analytic lower/upper geometric estimates.'}
    (ROOT/'family_windows.json').write_text(json.dumps(out,indent=2))
    valid=sum(t['valid_second_order'] for x in audit for t in x['second_order_tests'])
    print('family_windows=%d first_order_all=%s valid_second_order=%d/%d' %
          (len(rows),all(r['strict'] for r in rows),valid,len(audit)*len(Q_AUDIT)))

if __name__=='__main__': main()
