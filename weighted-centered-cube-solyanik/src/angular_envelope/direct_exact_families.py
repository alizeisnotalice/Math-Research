#!/usr/bin/env python3
"""Exact direct parameter-family, strict-positive outer approximation, and n=1 audits."""
from fractions import Fraction as Q
from pathlib import Path
import json
ROOT=Path(__file__).resolve().parent
def r(q):return 1-3*q+3*q*q
def dq(q):return 3*q*(1-q)*(1-2*q)/r(q)
def G(d,q):return d*d*(3*q*(1-q)+d*(2*q-1)/2)/(4*q*q*(1-q)*(1-q))
def hparam(q):return 27*q*(1-q)*(1-2*q)**2*(1-2*q+2*q*q)/(8*r(q)**3)
def binary(q,d):
 delta=d/2;a=1-delta/q;b=1+delta/(1-q);D=q*a**3+(1-q)*b**3-1;L=q*abs(a-1)+(1-q)*abs(b-1)
 assert q*a+(1-q)*b==1 and L==d
 return a,b,D
def main():
 rows=[];fail=[]
 for q in (Q(17,48),Q(3,8),Q(2,5),Q(5,12),Q(7,16),Q(11,24),Q(23,48)):
  d=dq(q);a,b,D=binary(q,d);ok=(0<a<1<b and D==G(d,q)==hparam(q))
  rows.append({'type':'inner_exact_parameter','q':str(q),'d':str(d),'a':str(a),'b':str(b),'D':str(D),'H_param':str(hparam(q)),'pass':ok})
  if not ok:fail.append(rows[-1])
 for d in (Q(2,3),Q(3,4),Q(1),Q(3,2),Q(2)-Q(1,16)):
  H=d*(4-d)/(2-d)**2
  for k in (8,16,32):
   q=d/2+Q(1,2**k)
   if q>=1:continue
   a,b,D=binary(q,d);ok=(a>0 and D>H)
   rows.append({'type':'outer_strict_positive_approach','d':str(d),'k':k,'q':str(q),'a':str(a),'D_minus_H':str(D-H),'pass':ok})
   if not ok:fail.append(rows[-1])
 for d in (Q(0),Q(1,8),Q(1,3),Q(2,3),Q(7,8)):
  vals=[1-d,1+d];D=sum(z**3 for z in vals)/2-1;L=sum(abs(z-1) for z in vals)/2;ok=(D==3*d*d and L==d and (d==0 or d<1))
  rows.append({'type':'n1_equal_atoms','d':str(d),'D':str(D),'pass':ok})
  if not ok:fail.append(rows[-1])
 # Literal unit atom: complement is {f>=1}, so its mean includes f=1.
 w=[Q(1,3)]*3;x=[Q(1),Q(1,2),Q(3,2)];D=sum(a*b**3 for a,b in zip(w,x))-1;L=sum(a*abs(b-1) for a,b in zip(w,x))
 q=Q(1,3);a=Q(1,2);b=(Q(1)*w[0]+Q(3,2)*w[2])/(w[0]+w[2]);Tbin=q*a**3+(1-q)*b**3-1;dbin=q*(1-a)+(1-q)*(b-1)
 ok=(sum(a*b for a,b in zip(w,x))==1 and D==Q(1,2) and L==Q(1,3) and Tbin==Q(11,32) and L==dbin and D>=Tbin)
 rows.append({'type':'unit_atom','d':str(L),'D':str(D),'Tbin':str(Tbin),'distance_preserved':L==dbin,'pass':ok})
 if not ok:fail.append(rows[-1])
 # Nonnegative boundary f=(0,2): it is intentionally outside strict positivity.
 w=[Q(1,2),Q(1,2)];x=[Q(0),Q(2)];D=sum(a*b**3 for a,b in zip(w,x))-1;L=sum(a*abs(b-1) for a,b in zip(w,x));H=Q(1)*(4-Q(1))/(2-Q(1))**2
 ok=(sum(a*b for a,b in zip(w,x))==1 and L==1 and D==3 and H==3 and D==H)
 rows.append({'type':'nonnegative_boundary_excluded_from_strict_positive','q':'1/2','d':str(L),'values':['0','2'],'D':str(D),'H':str(H),'pass':ok})
 if not ok:fail.append(rows[-1])
 (ROOT/'direct_failed.json').write_text(json.dumps(fail,indent=2));(ROOT/'direct_summary.json').write_text(json.dumps({'rows':rows,'failed':fail},indent=2));assert not fail
 print('direct_exact_rows=%d failed=0' % len(rows))
if __name__=='__main__':main()
