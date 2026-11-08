import os
os.environ['OPENBLAS_NUM_THREADS']='1';os.environ['OMP_NUM_THREADS']='1'
import math,json,sys,time
from fractions import Fraction
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor
import numpy as np
from scipy.special import gammaln,logsumexp
from scipy.stats import binom
from scipy.optimize import brentq
from scipy.integrate import quad
MU=brentq(lambda x:x-1-math.exp(-x),1,2,xtol=1e-15)
def log_bracket(x,terms):
 if x==1:return Fraction(0),Fraction(0)
 t=(x-1)/(x+1);power=t;total=Fraction(0)
 for k in range(terms):
  total+=power/Fraction(2*k+1);power*=t*t
 total*=2;err=2*abs(power)/(Fraction(2*terms+1)*(1-t*t))
 return total-err,total+err

def ratio_sign(n,s,c,k,want):
 g=s-c*k
 if g<=c:return {'certified':want=='nonpositive','endpoint':True}
 a=Fraction(n-k,k+1);b=(g-c)/g
 # Exact ties occur e.g. at N=2, s=2c; detect them without rounding.
 if n<=4096 and a*b**(n-1)==1:return {'certified':True,'exact_tie':True}
 for terms in [32,64,128,256]:
  al,au=log_bracket(a,terms);bl,bu=log_bracket(b,terms);lo=al+(n-1)*bl;hi=au+(n-1)*bu
  if (want=='nonpositive' and hi<0) or (want=='nonnegative' and lo>0):
   return {'certified':True,'terms':terms,'log_ratio_lower_float':float(lo),'log_ratio_upper_float':float(hi)}
 raise AssertionError('Rational mode bracket did not certify requested sign')

def one(j):
 n=j['N'];reg=j['regime'];fac=j.get('refinement',1);z=brentq(lambda x:.5*math.erfc(x/math.sqrt(2))-math.sqrt(MU/n)*x,0,math.sqrt(math.log(n))+2,xtol=1e-14)
 s={'critical':n*MU-math.sqrt(MU*n)*z,'mean':n*MU,'far2':MU*n*n,'far3':MU*n*n*n}[reg];c=MU
 last=min(n,math.ceil(s/c)-1)
 while s-c*last<=0:last-=1
 def lr(k):
  g=s-c*k
  if g<=c:return -math.inf
  return math.log(n-k)-math.log(k+1)+(n-1)*math.log1p(-c/g)
 lo,hi=0,last
 while lo<hi:
  mid=(lo+hi)//2
  if lr(mid)>0:lo=mid+1
  else:hi=mid
 mode=lo;sr=Fraction.from_float(s);cr=Fraction.from_float(c);right_cert=ratio_sign(n,sr,cr,mode,'nonpositive') if mode<last else {'certified':True,'endpoint':True};left_cert=ratio_sign(n,sr,cr,mode-1,'nonnegative') if mode>0 else {'certified':True,'endpoint':True};mode_right=lr(mode) if mode<last else -math.inf;mode_left=lr(mode-1) if mode>0 else math.inf;assert mode_right<=1e-12 and mode_left>=-1e-12
 rad=math.ceil(math.sqrt(110*(n+2)));k=np.arange(max(0,mode-rad),min(last,mode+rad)+1);p=(mode+.5)/(n+1);g0=s-c*mode
 logw=np.log(binom.pmf(k,n,p))-(k-mode)*math.log(p/(1-p))+(n-1)*np.log1p(-c*(k-mode)/g0)
 w=np.exp(logw-logsumexp(logw));g=s-c*k;B=float((n-1)*(w/g).sum());T=min(s,80.)
 def fn(t):
  ok=g>t
  if not np.any(ok):return 0.
  return math.exp(-t)*float(w[ok]@np.exp((n-1)*np.log1p(-t/g[ok])))
 points=g[(g>0)&(g<T)&(w>1e-16)]
 value,err=quad(fn,0,T,points=np.sort(points),epsabs=2e-12*fac,epsrel=2e-11*fac,limit=400)
 return {'N':n,'regime':reg,'refinement':fac,'s':s,'z':z,'mode_certificate_scope':'exact binary-rational s and c; not interval certification of all subsequent floating calculations','mode_certified':right_cert['certified'] and left_cert['certified'],'mode_certificate_right':right_cert,'mode_certificate_left':left_cert,'mode':mode,'mode_log_ratio_right':None if not math.isfinite(mode_right) else mode_right,'mode_log_ratio_left':None if not math.isfinite(mode_left) else mode_left,'retained':len(k),'posterior_tail_bound':(n+1)*math.exp(-220) if len(k)<last+1 else 0.,'normal_root_residual':.5*math.erfc(z/math.sqrt(2))-math.sqrt(MU/n)*z,'D':value,'B_left_at_nodes':B,'D_minus_half':value-.5,'critical_correction_scaled':(value-.5+z/(4*math.sqrt(MU*n)))*n/math.log(n) if reg=='critical' else None,'D_minus_frozen_rate':value-1/(1+B),'integral_tail_bound':math.exp(-T) if T<s else 0.,'quadrature_error_estimate':err,'normalization_residual':float(w.sum()-1),'numeric':bool(math.isfinite(value) and 0<value<1 and err<1e-8)}
if __name__=='__main__':
 inp=Path(sys.argv[1]);out=Path(sys.argv[2]);out.mkdir(exist_ok=True,parents=True);t=time.monotonic()
 with ProcessPoolExecutor(max_workers=12) as pool:rows=list(pool.map(one,json.loads(inp.read_text())['columns']))
 assert all(x['numeric'] for x in rows)
 (out/'columns71.json').write_text(json.dumps({'mu':MU,'seconds':time.monotonic()-t,'workers':12,'rows':rows},indent=2)+'\n');print('PASS',len(rows),'posterior column calculations')
