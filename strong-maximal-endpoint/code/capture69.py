import os
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
os.environ.setdefault('OMP_NUM_THREADS','1')
import json,time,math,sys
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor
import numpy as np
from scipy.special import gammaln,logsumexp,gammainc,gammaincc
from scipy.stats import binom
from scipy.integrate import quad
from scipy.optimize import brentq
MU=brentq(lambda c:c-1-math.exp(-c),1,2,xtol=1e-15)
def one(job):
 n,K=job['N'],job['K'];action=job.get('action','mu');refinement=job.get('refinement',1);start=time.monotonic();c=MU if action=='mu' else MU-math.sqrt(MU)/2*math.sqrt(math.log(n)/n);p=1/(1+math.exp(c));rad=math.sqrt(110*n);k=np.arange(max(0,math.floor(n*p-rad)),min(n,math.ceil(n*p+rad))+1);w=binom.pmf(k,n,p);lw=np.log(w);shift=c*k;lg=gammaln(n);mean=n*MU;sd=math.sqrt(n*MU);L=math.log(n)
 z=brentq(lambda z:.5*math.erfc(z/math.sqrt(2))-math.sqrt(MU/n)*z,0,math.sqrt(L)+2,xtol=1e-14);q0=mean-sd*z
 hole=None if K is None else (q0-K*L,q0+K*L)
 def logf(s):
  g=s-shift;ok=g>0;x=(g[ok]-n)/n;lp=np.log1p(x);const=-.5*math.log(2*math.pi*n)-1/(12*n)+1/(360*n**3)-1/(1260*n**5)+1/(1680*n**7)-1/(1188*n**9);return float(logsumexp(lw[ok]+const+n*(lp-x)-lp))
 def cdf(s):return float(w@gammainc(n,np.maximum(s-shift,0)))
 def logr(s):return s+logf(s)
 logn=math.log(n);logcap=(n-2)*math.log1p(math.log(2))+math.log(n+math.log(2))
 def inv(t):
  if t<=logn:return 1.
  if t>=logcap:return 2.
  y=min(math.log(2),math.expm1((t-logn)/(n-1)))
  for _ in range(8):
   r=(n-2)*math.log1p(y)+math.log(n+y)-t
   y=max(0.,min(math.log(2),y-r/((n-2)/(1+y)+1/(n+y))))
  return math.exp(y)
 def regions(lo,hi):
  if hi<=lo:return []
  if hole is None:return [(lo,hi)]
  a,b=hole
  return [(x,y) for x,y in [(lo,min(hi,a)),(max(lo,b),hi)] if y>x]
 def intquad(fun,lo,hi):
  return sum(quad(fun,a,b,epsabs=2e-12*refinement,epsrel=2e-10*refinement,limit=120)[0] for a,b in regions(lo,hi))
 upper=mean+14*sd
 def thresholds(a):
  la=logr(a)
  soft=brentq(lambda s:logr(s)-la-logn,a,a+5*L+20,xtol=2e-10)
  return la,soft
 def balance(a):
  la,soft=thresholds(a)
  missing=cdf(a)
  if hole is not None:
   lo=max(a,hole[0]);hi=hole[1]
   if hi>lo:missing+=cdf(hi)-cdf(lo)
  excess=intquad(lambda s:math.exp(logf(s))*(inv(logr(s)-la)-1),soft,upper)
  return excess-missing
 lo=max(1e-6,q0-20*L);hi=q0+20*L
 while balance(lo)<0:lo=max(1e-9,lo-20*L)
 while balance(hi)>0:hi+=20*L
 a=brentq(balance,lo,hi,xtol=4e-9*refinement,rtol=1e-13);la,soft=thresholds(a)
 cap=brentq(lambda s:logr(s)-la-logcap,soft,soft+2*n+100,xtol=2e-9)
 plateau=sum(math.exp(a-x)-math.exp(a-y) for x,y in regions(a,soft))
 def ratio(s):
  h=inv(logr(s)-la);y=math.log(h);return math.exp(logf(s))*h*(1+y)/(n+y)
 softint=intquad(ratio,soft,min(cap,upper));tail_bound=float(w@gammaincc(n,np.maximum(upper-shift,0)))
 capcost=math.exp(a-cap+math.log(2)+(n-1)*math.log1p(math.log(2)))
 if hole is not None and hole[1]>cap:
  x=max(cap,hole[0]);y=hole[1];capcost-=math.exp(a-x+math.log(2)+(n-1)*math.log1p(math.log(2)))*(-math.expm1(x-y))
 scaled=plateau+math.exp(a-la)*softint+capcost
 score=a-math.log(scaled)
 stationarity=max(abs((n-2)*math.log1p(math.log(inv(t)))+math.log(n+math.log(inv(t)))-t) for t in np.linspace(logn,logcap,17))
 assert stationarity<2e-8 and abs(balance(a))<1e-7 and scaled>0
 assert cdf(a)<=cdf(soft)+1e-12<=cdf(upper)+1e-12
 # Unrestricted comparison window shares, in scaled coordinates.
 shares={}
 if hole is None:
  for kk in [.5,1,2,4]:
   x=max(a,q0-kk*L);y=q0+kk*L
   pp=max(0.,math.exp(a-x)-math.exp(a-min(y,soft))) if min(y,soft)>x else 0.
   ss=quad(ratio,max(x,soft),min(y,cap),epsabs=2e-12*refinement,epsrel=2e-10*refinement,limit=120)[0] if min(y,cap)>max(x,soft) else 0.
   shares[str(kk)]=(pp+math.exp(a-la)*ss)/scaled
 return {'N':n,'K':K,'action':action,'refinement':refinement,'c':c,'q0':q0,'cutoff':a,'soft_start':soft,'log_lambda':-la,'log_score':score,'mean_residual':balance(a),'normal_root_residual':.5*math.erfc(z/math.sqrt(2))-math.sqrt(MU/n)*z,'prior_normalization_residual':float(w.sum()-1),'upper_integration_limit':upper,'upper_tail_mass_double_diagnostic':tail_bound,'omitted_soft_scaled_cost_bound':0. if cap<=upper else math.exp(a-la)*2*(1+math.log(2))/n*tail_bound,'prior_omitted_bound':2*math.exp(-220),'prior_retained_count':len(k),'inverse_derivative_residual':stationarity,'relative_ideal_dual_truncation_bound':2*math.exp(score-la)*2*math.exp(-220),'cost_scaled_by_exp_cutoff':scaled,'plateau_scaled_cost':plateau,'soft_scaled_cost':math.exp(a-la)*softint,'window_cost_shares':shares,'seconds':time.monotonic()-start,'status':'numerical candidate; no directed quadrature or KKT dual certificate'}
def main():
 inp=Path(sys.argv[1]);out=Path(sys.argv[2]);out.mkdir(exist_ok=True,parents=True);jobs=json.loads(inp.read_text())['jobs'];t=time.monotonic()
 with ProcessPoolExecutor(max_workers=12) as pool:rows=list(pool.map(one,jobs))
 (out/'capture69.json').write_text(json.dumps({'seconds':time.monotonic()-t,'workers':12,'rows':rows},indent=2)+'\n');print('Completed',len(rows),'capture candidates in',time.monotonic()-t,'seconds')
if __name__=='__main__':main()
