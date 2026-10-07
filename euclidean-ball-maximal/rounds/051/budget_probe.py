#!/usr/bin/env python3
"""Actual-budget constants: exact dyadic bounds, exploratory high-precision counts."""
from fractions import Fraction as F
from decimal import Decimal as D,localcontext,ROUND_CEILING
import json

def ceil_log2_integer(q):return (q-1).bit_length()
def counts(n,precision):
 with localcontext() as ctx:
  ctx.prec=precision;ln2=D(2).ln();answer=[]
  for a in range(1,n):
   H=-D(n)*(1-(-D(a)*ln2/n).exp()).ln()/ln2
   N=max(0,int(H.to_integral_value(rounding=ROUND_CEILING))-a-1)
   answer.append(N)
 return answer

def main():
 rows=[]
 for batch,ns in enumerate(((1,2,3,4,8),(16,32,64,128),(256,512,1024,4096))):
  for n in ns:
   B=n*ceil_log2_integer(2*n);L=min(n,ceil_log2_integer(4*B))
   ns50=counts(n,50);ns80=counts(n,80);assert ns50==ns80
   assert all(0<=N<=B for N in ns80)
   if n==2:assert ns80==[2]
   exact_upper=F(0) if L==n else 4*B*F(2)**(-L)
   assert exact_upper<=1
   # Counts are supported by 50/80-digit agreement; this is not interval arithmetic.
   scalar_tail=2*sum((F(N)*F(2)**(-a) for a,N in enumerate(ns80,1) if a>=L),F(0))
   assert scalar_tail<=exact_upper
   rows.append(dict(batch=batch,n=n,B=B,L=L,dimension_free_tail_constant_upper=str(exact_upper),count_based_tail_approx=float(scalar_tail),count_precision_digits=[50,80],count_max=max(ns80,default=0)))
 print(json.dumps(dict(status='passed',batch_counts=[5,4,4],cases=rows,scope='Integer constants and dyadic upper bounds exact. Transcendental admissible counts checked at two precisions, not certified intervals. All-n conclusion uses analytic concavity and geometric-series proof.',main_theorem_proved=False),indent=2))
if __name__=='__main__':main()
