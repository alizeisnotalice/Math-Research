#!/usr/bin/env python3
"""Independent exact verifier for finite rational p=3 vectors.

Input contract: r_i>0 rational, sum r_i=1; f_i>0 rational, sum r_i f_i=1.
It verifies T=E f^3-1 >= ((11+5 sqrt(5))/8)d^2 via rational A/square
certificate and separately checks the binary coarse-graining lower bound.
"""
from fractions import Fraction as Q
from pathlib import Path
import json

ROOT=Path(__file__).resolve().parent;DIMS=(3,4,8,16,64)

def normalize(w,x):
 m=sum((a*b for a,b in zip(w,x)),Q(0));return w,[z/m for z in x]
def random_vector(n,seed):
 rw=[Q((11*seed+7*i)%97+1,1) for i in range(n)];w=[z/sum(rw,Q(0)) for z in rw]
 x=[Q((17*seed+13*i)%89+1,(seed+3*i)%23+2) for i in range(n)];return normalize(w,x)
def critical_vector(n,q,a):
 lo=max(1,n//2);hi=n-lo;w=[q/Q(lo) for _ in range(lo)]+[(1-q)/Q(hi) for _ in range(hi)]
 b=(1-q*a)/(1-q);x=[a*(Q(9,8) if i%2==0 else Q(7,8)) for i in range(lo)];x[-1]=lo*a-sum(x[:-1],Q(0))
 y=[b*(Q(10,9) if i%2==0 else Q(8,9)) for i in range(hi)];y[-1]=hi*b-sum(y[:-1],Q(0));return w,x+y
def embedded_binary_vector(n,q,a):
 lo=max(1,n//2);hi=n-lo;w=[q/Q(lo) for _ in range(lo)]+[(1-q)/Q(hi) for _ in range(hi)]
 b=(1-q*a)/(1-q);return w,[a for _ in range(lo)]+[b for _ in range(hi)]
def equal_atom_vector(n):
 w=[Q(1,n) for _ in range(n)];x=[Q(1) for _ in range(n)];x[0]=Q(1,2);x[1]=Q(3,2);return w,x
def constant_vector(n): return [Q(1,n) for _ in range(n)],[Q(1) for _ in range(n)]
def law_signature(w,x):
 merged={}
 for r,z in zip(w,x):merged[z]=merged.get(z,Q(0))+r
 return '|'.join('%s:%s' % (z,merged[z]) for z in sorted(merged))
def sqrt5_certificate(T,d):
 # 8T-11d^2 >= 5sqrt(5)d^2 iff A>=0 and A^2>=125d^4.
 A=8*T-11*d*d; square=A*A-125*d**4
 return {'A':str(A),'square_gap':str(square),'A_nonnegative':A>=0,'square_nonnegative':square>=0,
         'pass':A>=0 and square>=0}
def coarsegrain(w,x):
 low=[i for i,z in enumerate(x) if z<1];high=[i for i,z in enumerate(x) if z>=1]
 assert low and high
 q=sum((w[i] for i in low),Q(0));a=sum((w[i]*x[i] for i in low),Q(0))/q;b=sum((w[i]*x[i] for i in high),Q(0))/(1-q)
 Tbin=q*a**3+(1-q)*b**3-1;d=sum((r*abs(z-1) for r,z in zip(w,x)),Q(0));dbin=q*(1-a)+(1-q)*(b-1)
 return q,a,b,Tbin,d,dbin
def verify(label,w,x):
 assert all(r>0 for r in w) and all(z>0 for z in x) and sum(w,Q(0))==1 and sum((r*z for r,z in zip(w,x)),Q(0))==1
 T=sum((r*z**3 for r,z in zip(w,x)),Q(0))-1;d=sum((r*abs(z-1) for r,z in zip(w,x)),Q(0))
 if d==0:
  assert all(z==1 for z in x) and T==0
  return {'label':label,'n':len(w),'law_signature':law_signature(w,x),'T':'0','d':'0','zero_distance_branch':True,'pass':True}
 q,a,b,Tbin,d2,dbin=coarsegrain(w,x);full=sqrt5_certificate(T,d);binary=sqrt5_certificate(Tbin,dbin)
 return {'label':label,'n':len(w),'law_signature':law_signature(w,x),'T':str(T),'d':str(d),'certificate':full,'binary_certificate':binary,
         'q':str(q),'a':str(a),'b':str(b),'jensen_gap':str(T-Tbin),'d_equals_binary_d':d==dbin and d==d2,
         'pass':full['pass'] and binary['pass'] and T>=Tbin and d==dbin}
def main():
 rows=[];failed=[]
 for n in DIMS:
  for k in range(12): rows.append(verify('random',*random_vector(n,1000+n*31+k)))
  for q in (Q(381964,10**6),Q(381966,10**6),Q(381968,10**6)):
   for k in (8,16,32):
    rows.append(verify('critical_embedded',*embedded_binary_vector(n,q,Q(1,2**k))))
    rows.append(verify('critical_spread',*critical_vector(n,q,Q(1,2**k))))
  rows.append(verify('equal_atom',*equal_atom_vector(n)))
  rows.append(verify('constant',*constant_vector(n)))
 for r in rows:
  if not r['pass']:failed.append(r)
 laws={r['law_signature'] for r in rows};nonzero_laws={r['law_signature'] for r in rows if not r.get('zero_distance_branch',False)}
 out={'input_contract':'positive rational r,f; exact sums E1=Ef=1','law_definition':'merge equal f then sort(value,totalmass)','constant':'(11+5sqrt(5))/8','cases':len(rows),'unique_probability_laws':len(laws),'unique_nonconstant_probability_laws':len(nonzero_laws),'failed':failed,'records':rows}
 (ROOT/'exact_failed_inputs.json').write_text(json.dumps(failed,indent=2));(ROOT/'exact_summary.json').write_text(json.dumps(out,indent=2));assert not failed
 print('exact_multi_cases=%d failed=0 dimensions=%s' % (len(rows),DIMS))
if __name__=='__main__':main()
