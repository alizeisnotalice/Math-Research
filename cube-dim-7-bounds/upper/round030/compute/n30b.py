#!/usr/bin/env python3
"""Exact rational/log N30B certification; fixed-t search and certification separated."""
import os,sys,time,json,hashlib,resource
from fractions import Fraction as Q
from math import comb,isqrt
from pathlib import Path
sys.dont_write_bytecode=True
root=Path(__file__).resolve().parents[3];sys.path.insert(0,str(root/'revision_v3/numeric'))
from verify import log_bounds, floorq,ceilq
D10=10**180

def sqrtbounds(x):
 assert x>=0
 d=10**190;m=isqrt(x.numerator*d*d//x.denominator)
 return Q(m,d),Q(m+1,d)
def signed_log(c,lo,hi):return (c*lo,c*hi) if c>=0 else (c*hi,c*lo)
def Bbounds(n,t):
 q=1+1/(2*t);h=(1+1/t)/2;c=n-Q(n,2)*q**(n-1)
 terms=sum((Q(comb(n-1,j),j)*q**(n-1-j)*Q(-1,2)**j*(1-t**(-j)) for j in range(1,n)),Q())
 rat=1-h**n-Q(n,2)*terms;lo,hi=signed_log(c,*log_bounds(t))
 return floorq(rat+lo,D10),ceilq(rat+hi,D10)
def Pbounds(n,t,s):
 assert Q(1,t)<s<1
 q=1+1/(2*t);R=t*s;c=n-Q(n,2)*q**(n-1)
 terms=sum((Q(comb(n-1,j),j)*q**(n-1-j)*Q(-1,2)**j*(R**(-j)-t**(-j)) for j in range(1,n)),Q())
 rat=-Q(n,2)*terms;lo,hi=signed_log(c,*log_bounds(1/s))
 return floorq(rat+lo,D10),ceilq(rat+hi,D10)
def f_bounds(A,Bl,Bu):
 l=sqrtbounds(A*Bl)[0];u=sqrtbounds(A*Bu)[1]
 return floorq(A+Bl+2*l,D10),ceilq(A+Bu+2*u,D10)
def certify(n,t,s):
 B,Bu=Bbounds(n,t); A=(1+1/t)**n; D=s**n-((s+1/t)/2)**n; P,Pu=Pbounds(n,t,s)
 assert D>0 and B>D
 tau=sqrtbounds(A*(Bu-D))[1];oldlo,oldhi=f_bounds(A,B,Bu);newlo,newhi=f_bounds(A,B-D,Bu-D)
 return {'B_lower':str(B),'B_upper':str(Bu),'P_lower':str(P),'P_upper':str(Pu),'D':str(D),'tau_upper':str(tau),'P_minus_D_gt_tau_upper':P-D>tau,'new_F_lower':str(newlo),'new_F_upper':str(newhi),'old_F_lower':str(oldlo),'old_F_upper':str(oldhi),'new_upper_lt_old_lower':newhi<oldlo,'gain_lower':str(oldlo-newhi)}
def frozen():
 n=7;t=Q(3391749,524288);s=Q(5,3)/t;R=t*s
 # The frozen proof uses sigma=5/(3t) and a separate guard R=7/4.
 # Here only the limiting R=t*sigma cases 7/4 and 5/3 are compared; see n30b_frozen_guard.py for the original guard.
 results=[]
 for R in [Q(7,4),Q(5,3)]:
  s=R/t;r=certify(n,t,s);B=Q(r['B_lower']);D=Q(r['D']); loglo,loghi=log_bounds(R);coarse=B-1-n*loghi;tau=Q(r['tau_upper'])
  if R==Q(5,3): assert D==Q(61741,2187)/t**7
  r.update({'R':str(R),'sigma':str(s),'coarse_P_lower_minus_D_gt_tau_upper':coarse-D>tau,'coarse_P_lower_decimal':str(float(coarse))})
  results.append(r)
 return results
if __name__=='__main__':
 start=time.monotonic(); cpu=time.process_time();res=frozen();out={'batch':'N30B frozen','case_scope':'Limiting R=t*sigma cases only; frozen sigma=5/(3t) with separate R=7/4 guard is valid and checked in n30b_frozen_guard.py','results':res,'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'pid':os.getpid(),'wall_seconds':time.monotonic()-start,'cpu_seconds':time.process_time()-cpu,'peak_rss_bytes':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
 Path(__file__).with_name('n30b_frozen.json').write_text(json.dumps(out,indent=2)+'\n'); print(json.dumps([{'R':r['R'],'PminusD_gt_tau':r['P_minus_D_gt_tau_upper'],'coarse':r['coarse_P_lower_minus_D_gt_tau_upper'],'new_lt_old':r['new_upper_lt_old_lower'],'gain':float(Q(r['gain_lower']))} for r in res],indent=2))
