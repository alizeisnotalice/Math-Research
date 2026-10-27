#!/usr/bin/env python3
"""L30A corrected: original integrand and exact rational interval certificates."""
import os,sys,json,time,hashlib,resource
from pathlib import Path
from fractions import Fraction as Q
from math import comb
sys.dont_write_bytecode=True
import mpmath as mp
import sympy as sp
root=Path(__file__).resolve().parents[3];sys.path.insert(0,str(root/'revision_v3/numeric'))
from verify import certify,sqrt_bounds,log_bounds
start=time.monotonic(); cpu=time.process_time();mp.mp.dps=80
n=7;t=Q(3391749,524288);beta=1+1/t;q=(1+beta)/2;Rstar=2/beta
b,s=sp.symbols('beta s');ell=(1+b-s)/2
assert sp.expand((b-s)-(2*ell-1))==0
assert sp.expand((sp.Rational(1,2)+6*ell-3*(b-s)*ell**5)-(sp.Rational(1,2)+6*ell+3*ell**5-6*ell**6))==0
coeff=[(Q(1,2) if j==0 else Q(0))+(6*q if j==0 else Q(-3) if j==1 else Q(0))+(3*comb(5,j)*q**(5-j)*Q(-1,2)**j if j<=5 else Q(0))-6*comb(6,j)*q**(6-j)*Q(-1,2)**j for j in range(7)]
poly=sp.Poly(sp.expand((sp.Rational(1,2)+6*ell-3*(b-s)*ell**5).subs(b,sp.Rational(beta.numerator,beta.denominator))),s)
assert all(poly.nth(j)==sp.Rational(c.numerator,c.denominator) for j,c in enumerate(coeff))
A,Bl,Bu,_,_=certify(n,t);tau_lo=sqrt_bounds(A*Bl)[0];tau_hi=sqrt_bounds(A*Bu)[1]
mt=mp.mpf(t.numerator)/t.denominator; mb=mp.mpf(beta.numerator)/beta.denominator;mq=mp.mpf(q.numerator)/q.denominator
results=[]
for R in [Rstar,Q(7,4),Q(9,5)]:
 assert R>=Rstar and R<t
 mR=mp.mpf(R.numerator)/R.denominator
 def original(x):
  ell0=1+1/(2*mt)-1/(2*x)
  return (mp.mpf('0.5')+6*ell0-3*(mb-1/x)*ell0**5)/x
 v=mp.quad(original,[mR,mt]);ll,lu=log_bounds(t/R)
 rat=sum((coeff[j]*(R**(-j)-t**(-j))/j for j in range(1,7)),Q())
 lo=rat+coeff[0]*(ll if coeff[0]>=0 else lu);hi=rat+coeff[0]*(lu if coeff[0]>=0 else ll)
 assert lo<=hi
 assert abs(v-mp.mpf(lo.numerator)/lo.denominator)<mp.mpf("1e-65")
 results.append({'R':str(R),'quad_80dps':mp.nstr(v,60),'L_lower':str(lo),'L_upper':str(hi),'L_lower_decimal':mp.nstr(mp.mpf(lo.numerator)/lo.denominator,40),'L_upper_decimal':mp.nstr(mp.mpf(hi.numerator)/hi.denominator,40),'L_lower_gt_tau_upper':lo>tau_hi,'L_upper_lt_tau_lower':hi<tau_lo,'gap_lower_minus_tau_upper':str(lo-tau_hi)})
out={'batch':'L30A corrected','previous_error':'previously substituted beta-1/rho=2ell; correct identity is 2ell-1; old L30A exclusion invalid','symbolic_identity_checks':'PASS sympy exact','t':str(t),'beta':str(beta),'Rstar':str(Rstar),'integrand':'(1/2+6ell+3ell^5-6ell^6)/rho','coefficients':[str(c) for c in coeff],'tau_lower':str(tau_lo),'tau_upper':str(tau_hi),'tau_upper_decimal':mp.nstr(mp.mpf(tau_hi.numerator)/tau_hi.denominator,40),'points':results,'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'pid':os.getpid(),'wall_seconds':time.monotonic()-start,'cpu_seconds':time.process_time()-cpu,'peak_rss_bytes':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
Path(__file__).with_name('l30a_corrected_result.json').write_text(json.dumps(out,indent=2)+'\n'); print(json.dumps({'symbolic_identity_checks':out['symbolic_identity_checks'],'tau_upper_decimal':out['tau_upper_decimal'],'points':[{k:v for k,v in x.items() if k not in ('L_lower','L_upper','gap_lower_minus_tau_upper')} for x in results],'pid':out['pid'],'wall_seconds':out['wall_seconds']},indent=2))
