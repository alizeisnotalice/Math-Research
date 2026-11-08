import os
os.environ['OPENBLAS_NUM_THREADS']='1';os.environ['OMP_NUM_THREADS']='1'
import math,json,sys,time
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor
import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.special import gammaln,gammainc,gammaincc,logsumexp
from scipy.stats import binom
from scipy.integrate import quad
from scipy.optimize import brentq
MU=brentq(lambda x:x-1-math.exp(-x),1,2,xtol=1e-15)

def law(n,c):
 p=1/(1+math.exp(c));rad=math.ceil(math.sqrt(220*n));k=np.arange(max(0,math.floor(n*p)-rad),min(n,math.ceil(n*p)+rad)+1);lw=binom.logpmf(k,n,p);mass=float(np.exp(lw).sum());lw-=logsumexp(lw);w=np.exp(lw);shift=c*k
 def logf(s):
  g=s-shift;ok=g>0;x=(g[ok]-n)/n;lp=np.log1p(x);const=-.5*math.log(2*math.pi*n)-1/(12*n)+1/(360*n**3)-1/(1260*n**5)+1/(1680*n**7)-1/(1188*n**9)
  return float(logsumexp(lw[ok]+const+n*(lp-x)-lp))
 def cdf(s):return float(w@gammainc(n,np.maximum(s-shift,0)))
 def sf(s):return float(w@gammaincc(n,np.maximum(s-shift,0)))
 return logf,cdf,sf,mass

def post(n,c,s):
 last=min(n,math.ceil(s/c)-1)
 def lr(k):
  g=s-c*k
  return -math.inf if g<=c else math.log(n-k)-math.log(k+1)+(n-1)*math.log1p(-c/g)
 lo,hi=0,last
 while lo<hi:
  k=(lo+hi)//2
  if lr(k)>0:lo=k+1
  else:hi=k
 mode=lo;rad=math.ceil(math.sqrt(110*(n+2)));k=np.arange(max(0,mode-rad),min(last,mode+rad)+1);p=(mode+.5)/(n+1);lw=binom.logpmf(k,n,p)-(k-mode)*math.log(p/(1-p))+(n-1)*np.log1p(-c*(k-mode)/(s-c*mode));w=np.exp(lw-logsumexp(lw));return k,w

def run(j):
 start=time.monotonic();n=j['N'];refine=j['refinement'];L=math.log(n);rho=math.sqrt(L/n);D=1.5 if j['D_type']=='constant' else L**(-.25);cr=MU-math.sqrt(MU)/2*rho;ca=cr+D*rho;cap=1.5;delta=1/(n*L**2);epsilon=1e-8
 fa,Fa,Sa,massa=law(n,ca);fr,Fr,Sr,massr=law(n,cr);logn=math.log(n);logcap=(n-2)*math.log1p(math.log(cap))+math.log(n+math.log(cap));mean=n*MU;sd=math.sqrt(n*MU);z=brentq(lambda z:.5*math.erfc(z/math.sqrt(2))-math.sqrt(MU/n)*z,0,math.sqrt(L)+2);q0=mean-sd*z
 def logr(s):return s+fa(s)
 def inv(t):
  if t<=logn:return 1.
  if t>=logcap:return cap
  y=min(math.log(cap),math.expm1((t-logn)/(n-1)))
  for _ in range(8):y=max(0.,min(math.log(cap),y-((n-2)*math.log1p(y)+math.log(n+y)-t)/((n-2)/(1+y)+1/(n+y))))
  return math.exp(y)
 def threshold(a):
  la=logr(a);soft=brentq(lambda s:logr(s)-la-logn,a,a+5*L+20,xtol=2e-10);return la,soft
 upper=mean+14*sd
 def balance(a):
  la,soft=threshold(a);excess=quad(lambda s:math.exp(fa(s))*(inv(logr(s)-la)-1),soft,upper,epsabs=2e-13*refine,epsrel=2e-11*refine,limit=160)[0];return excess-Fa(a)
 a=brentq(balance,max(1e-6,q0-20*L),q0+20*L,xtol=2e-9*refine,rtol=2e-14);la,soft=threshold(a)
 def h(s):return 0. if s<a else 1. if s<=soft else inv(logr(s)-la)
 def costratio(s):
  hv=h(s);y=math.log(hv);return math.exp(fa(s))*hv*(1+y)/(n+y)
 ss=quad(costratio,soft,upper,epsabs=2e-13*refine,epsrel=2e-11*refine,limit=160)[0];scaled=-math.expm1(a-soft)+math.exp(a-la)*ss;logJ=-a+math.log(scaled)
 # Stable difference of the radial baseline means; subtract the unit profile.
 drift=-(Fr(a)-Fa(a))+quad(lambda s:(math.exp(fr(s))-math.exp(fa(s)))*(h(s)-1),soft,upper,epsabs=2e-14*refine,epsrel=2e-10*refine,limit=160)[0]
 nodes,weights=leggauss(32 if refine==1 else 64);deletion=0.;actualcheck=0.;rows=[]
 for node,weight in zip(nodes,weights):
  s=mean+sd*float(node);ka,wa=post(n,ca,s);kr,wr=post(n,cr,s);tails=np.cumsum(wa[::-1])[::-1];idx=int(np.flatnonzero(tails<=delta)[0]);cut=int(ka[idx]-1);tailabove=float(tails[idx]);boundarymass=float(wa[idx-1]);frac=(delta-tailabove)/boundarymass;refremoved=float(wr[kr>cut].sum())+frac*float(wr[kr==cut].sum());actremoved=tailabove+frac*boundarymass
  deletion+=float(weight)*sd*math.exp(fr(s))*h(s)*(refremoved-delta)/(1-delta)
  actualcheck=max(actualcheck,abs(actremoved-delta));rows.append(dict(S=s,cut=cut,fraction=frac,actual_removed=actremoved,reference_removed=refremoved,weight=float(weight)))
 refmean=(1+epsilon)*(1+balance(a)+drift-deletion);actualmean=(1+epsilon)*(1+balance(a));logcostupper=logJ+n*(math.log1p(epsilon)-math.log1p(-delta));tailbound=cap*(Sa(upper)+Sr(upper))
 return dict(j,c_actual=ca,c_reference=cr,D=D,delta=delta,strict_multiplier=1+epsilon,cap=cap,cutoff=a,soft_start=soft,log_R_cutoff=la,log_scalar_cost=logJ,log_actual_cost_upper=logcostupper,baseline_reference_mean_drift=drift,deletion_reference_mean_loss=deletion,reference_mean=refmean,reference_deficit_over_L_over_N=(1-refmean)*n/L,actual_mean=actualmean,mean_balance_residual=balance(a),conditional_mean_residual=actualcheck,radial_tail_bound_diagnostic=tailbound,prior_mass_residual_actual=massa-1,prior_mass_residual_reference=massr-1,prior_omission_bound=2*math.exp(-440),posterior_omission_bound=(n+1)*math.exp(-220),upper_integration_limit=upper,compact_S_cutoff=2*MU*n,quadrature_nodes=rows,seconds=time.monotonic()-start,scope='full radial floating candidate with posterior soft deletion; no interval certificate')
if __name__=='__main__':
 out=Path(sys.argv[1]);out.mkdir(exist_ok=True,parents=True);jobs=[dict(N=n,D_type=t,refinement=r) for n in [1024,4096,16384,65536,262144,1048576] for t in ['constant','vanishing'] for r in [1,.25]];start=time.monotonic()
 with ProcessPoolExecutor(max_workers=12) as pool:rows=list(pool.map(run,jobs))
 (out/'profile74.json').write_text(json.dumps(dict(workers=12,seconds=time.monotonic()-start,rows=rows),indent=2)+'\n');print('Completed',len(rows),'full radial profiles')
