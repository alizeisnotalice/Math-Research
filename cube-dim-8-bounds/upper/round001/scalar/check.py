#!/usr/bin/env python3
"""Recompute all Arb rows and independently bound the winning F using Fractions only."""
from fractions import Fraction as Q
from math import comb,isqrt
from pathlib import Path
from flint import arb,fmpq,ctx
import json, time
D=Path(__file__).resolve().parent; start=time.perf_counter(); ctx.prec=512
r=json.loads((D/'result.json').read_text());trace=json.loads((D/'trace.json').read_text())
assert len(trace)==5+r['optimizer_nfev'] and r['optimizer_nfev']<=80
assert [float.fromhex(v['t_hex']) for v in trace[:5]]==[6,7,8,9,10]
best=min(trace,key=lambda v:float.fromhex(v['F_hex'])); assert best==r['best_observed']
k=round(Q.from_float(float.fromhex(best['t_hex']))*2**20); assert k==r['grid_center_integer']
def coeff(t):
 q=1+Q(1,2)/t;h=(1+1/t)/2
 return 1-h**8-4*sum(Q(comb(7,j))*q**(7-j)*Q(-1,2)**j*(1-1/t**j)/j for j in range(1,8)),8-4*q**7
def a(x):return arb(fmpq(x.numerator,x.denominator))
def endpoint(x):
 m,e=x.man_exp();return Q(int(m))*Q(2)**int(e)
for i,row in enumerate(r['rows']):
 assert row['numerator_on_grid']==k+i-1 and Q(row['t'])==Q(k+i-1,2**20)
 t=Q(row['t']);c,d=coeff(t); B=a(c)+a(d)*a(t).log();assert B>0
 F=(a((1+1/t)**4)+B.sqrt())**2
 for key,v in [('B',B),('F',F)]:
  assert Q(row[key]['lower'])<=endpoint(v.lower())<=endpoint(v.upper())<=Q(row[key]['upper'])
w=min(r['rows'],key=lambda row:Q(row['F']['upper'])); assert w==r['winner']
assert all(Q(w['F']['upper'])<Q(row['F']['lower']) for row in r['rows'] if row!=w)
u=Q(r['strict_rational_upper']); assert Q(w['F']['upper'])<u and u==Q(Q(w['F']['upper']).numerator*10**10//Q(w['F']['upper']).denominator+1,10**10)
# For x>=1, log x = 2 sum z^(2j+1)/(2j+1), z=(x-1)/(x+1).
# Positive tail <= 2*z^(2N+1)/((2N+1)*(1-z*z)).
def loginterval(x,N=180):
 z=(x-1)/(x+1); assert 0<=z<1
 term=z; s=Q(0)
 for j in range(N):s+=term/(2*j+1);term*=z*z
 return 2*s,2*s+2*term/((2*N+1)*(1-z*z))
t=Q(w['t']); power=0; reduced=t
while reduced>=2:reduced/=2;power+=1
l2,u2=loginterval(Q(2));l,u0=loginterval(reduced);l+=power*l2;u0+=power*u2
c,d=coeff(t);blo,bhi=sorted((c+d*l,c+d*u0));assert blo>0
def sqrt_floor(x,bits=600):
 n=isqrt((x.numerator << (2*bits))//x.denominator)
 return Q(n,2**bits),Q(n+1,2**bits)
sl,_=sqrt_floor(blo);_,su=sqrt_floor(bhi);base=(1+1/t)**4
flo=(base+sl)**2;fhi=(base+su)**2
assert flo<=Q(w['F']['upper']) and fhi>=Q(w['F']['lower']) and fhi<u
print(json.dumps({'status':'PASS','arb_rows':3,'strict_winner_separation':True,'fraction_log_terms':180,'integer_sqrt_bits':600,'fraction_winner_overlap_arb':True,'fraction_winner_strict_upper':str(u),'wall_seconds':time.perf_counter()-start},indent=2))
