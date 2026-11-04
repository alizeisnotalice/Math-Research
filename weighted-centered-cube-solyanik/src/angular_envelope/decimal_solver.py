#!/usr/bin/env python3
"""Numerical-only Decimal candidate roots, independent of Fraction verifier."""
from decimal import Decimal,getcontext
from fractions import Fraction as Q
from pathlib import Path
import json
getcontext().prec=90
ROOT=Path(__file__).resolve().parent

def ds():
 return [Q(1,1000),Q(1,10),Q(1,2),Q(13,20),Q(2,3)-Q(1,2**20),Q(2,3),Q(7,10),Q(3,2),Q(2)-Q(1,2**30)]
def dec(x):return Decimal(x.numerator)/Decimal(x.denominator)
def dq(q):return 3*q*(1-q)*(1-2*q)/(1-3*q+3*q*q)
def G(d,q):return d*d*(3*q*(1-q)+d*(2*q-1)/2)/(4*q*q*(1-q)*(1-q))
def solve(d):
 if d>=Decimal(2)/3:return d/2,d*(4-d)/(2-d)**2
 l=Decimal(1)/3;u=Decimal(1)/2
 for _ in range(300):
  m=(l+u)/2
  if dq(m)>d:l=m
  else:u=m
 q=(l+u)/2;return q,G(d,q)
def main():
 rows=[]
 for x in ds():
  d=dec(x);q,H=solve(d);rows.append({'d':str(x),'q_decimal':str(q),'H_decimal':str(H),'branch':'inner' if d<Decimal(2)/3 else 'outer'})
 (ROOT/'decimal_candidates.json').write_text(json.dumps({'precision':getcontext().prec,'scope':'numerical candidates only','rows':rows},indent=2))
 print('decimal_candidates=%d precision=%d' % (len(rows),getcontext().prec))
if __name__=='__main__':main()
