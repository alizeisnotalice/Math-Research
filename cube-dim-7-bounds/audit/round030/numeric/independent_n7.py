#!/usr/bin/env python3
"""Independent fixed-point rational audit; no producer imports or search."""
from fractions import Fraction as Q
from math import comb,isqrt
from pathlib import Path
import json,hashlib,time,resource,os
start=time.monotonic(); cpu=time.process_time()
project=Path(__file__).resolve().parents[3]
src=project/'upper/round030/compute'
old=project/'upper/round028/cap_candidate/result.json'
files=[src/'n30b_certificates.json',src/'n30b.py',src/'N30B.md',old,Path(__file__)]
sha={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in files}
n=7;t=Q(3391749,524288);s=Q(2081927,6783498);U=Q(12277313907857,500000000000)
# log x = 2 sum z^(2k+1)/(2k+1); tail <= 2 z^(2N+1)/((2N+1)(1-z^2)).
def log_interval(x,N=500):
 assert x>=1
 z=(x-1)/(x+1); p=z; total=Q(0)
 for k in range(N):
  total+=p/(2*k+1); p*=z*z
 return 2*total,2*total+2*p/((2*N+1)*(1-z*z))
def scale(a,iv): return (a*iv[0],a*iv[1]) if a>=0 else (a*iv[1],a*iv[0])
# ell(rho)=q-1/(2rho). Integrate each monomial explicitly.
q=1+1/(2*t);h=(1+1/t)/2
c=n-Q(n,2)*q**(n-1)
def rat_term(lo,hi):
 return sum((Q(comb(n-1,j),j)*q**(n-1-j)*(-Q(1,2))**j*(lo**(-j)-hi**(-j)) for j in range(1,n)),Q())
def enclosure(rat,coef,li):
 l,u=scale(coef,li);return rat+l,rat+u
B=enclosure(1-h**n-Q(n,2)*rat_term(Q(1),t),c,log_interval(t))
P=enclosure(-Q(n,2)*rat_term(t*s,t),c,log_interval(1/s))
A=(1+1/t)**n;D=s**n-((s+1/t)/2)**n
# Outward decimal rational enclosures; comparisons are exact fractions.
S=10**40
def out(iv):return (Q(iv[0].numerator*S//iv[0].denominator,S),Q(-(-iv[1].numerator*S//iv[1].denominator),S))
B=out(B);P=out(P)
def sqrt_bound(x):
 assert x>=0
 m=isqrt(x.numerator*S*S//x.denominator)
 return Q(m,S),Q(m+1,S)
def cost(bl,bu):
 lo=sqrt_bound(A*bl)[0];hi=sqrt_bound(A*bu)[1]
 return out((A+bl+2*lo,A+bu+2*hi))
new=cost(B[0]-D,B[1]-D); oldF=cost(*B); tau=sqrt_bound(A*(B[1]-D))[1]
raw=json.loads(old.read_text())['candidate_interval']['lower']; mant=int(raw[0]); exp=raw[1]; arb=Q(mant*2**exp) if exp>=0 else Q(mant,2**(-exp))
checks={'one_over_t_lt_sigma_lt_one':1/t<s<1,'D_positive':D>0,'B_lower_gt_D':B[0]>D,'P_lower_minus_D_gt_tau_upper':P[0]-D>tau,'new_upper_lt_U':new[1]<U,'U_lt_round028_arb_lower':U<arb,'new_upper_lt_same_t_old_lower':new[1]<oldF[0]}
assert all(checks.values()),checks
def fmt(iv):return [str(x) for x in iv]
def approx(x):return format(float(x),'.16g')
result={'n':n,'t':str(t),'sigma':str(s),'U':str(U),'A':str(A),'D':str(D),'B':fmt(B),'P':fmt(P),'tau_upper':str(tau),'F_new':fmt(new),'F_old_same_t':fmt(oldF),'round028_arb_lower_exact':str(arb),'round028_encoding':'mantissa times 2**exponent','checks':checks,'approx':{k:approx(v) for k,v in {'B_lower':B[0],'P_lower':P[0],'D':D,'tau_upper':tau,'F_new_upper':new[1],'U':U,'F_old_lower':oldF[0],'P_minus_D_minus_tau':P[0]-D-tau,'same_t_gain':oldF[0]-new[1]}.items()},'log_series_terms':500,'output_denominator':str(S),'sha256':sha,'pid':os.getpid(),'wall_seconds':time.monotonic()-start,'cpu_seconds':time.process_time()-cpu,'peak_rss_bytes':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
Path(__file__).with_name('independent_n7.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({'checks':checks,'approx':result['approx'],'wall_seconds':result['wall_seconds'],'cpu_seconds':result['cpu_seconds'],'peak_rss_bytes':result['peak_rss_bytes'],'pid':result['pid']},indent=2))
