import itertools,json,math
from fractions import Fraction as F
from pathlib import Path
import numpy as np
from scipy.special import i0e
out=Path(__file__).resolve().parent
rng=np.random.default_rng(202610295199); errs=[]; caps=[]
for n in (1,2,4,8):
 for _ in range(4):
  probs=[F(int(v),100) for v in rng.integers(5,96,n)]
  dp=[F(1)]+[F(0)]*n
  for p in probs:
   dp=[dp[k]*(1-p)+(dp[k-1]*p if k else 0) for k in range(n+1)]
  brute=[F(0)]*(n+1)
  for labels in itertools.product((0,1),repeat=n):
   prob=F(1)
   for p,e in zip(probs,labels):prob*=p if e else 1-p
   brute[sum(labels)]+=prob
  assert dp==brute and sum(dp)==1
  V=sum(p*(1-p) for p in probs)
  caps.append({'n':n,'max_probability':float(max(dp)),'fourier_cap':float(i0e(float(V)))})
  assert float(max(dp))<=i0e(float(V))+1e-14
# True torus n=2: trigonometric quadrature integrates all relevant finite products exactly.
N=128; x=np.arange(N)/N
for L in (10/9,7/6,11/9):
 y=1+.9*np.sinc(L)*np.cos(2*np.pi*x); p=.5+.05*(y-1)/y
 for k in (0,1,2):
  q= ((1-p[:,None])*(1-p[None,:]) if k==0 else p[:,None]*p[None,:] if k==2 else p[:,None]*(1-p[None,:])+(1-p[:,None])*p[None,:])
  G=y[:,None]*y[None,:]*q/(math.comb(2,k)/4)
  errs.append(abs(float(G.mean())-1))
assert max(errs)<1e-13
tau_checks=0
for m,M in ((F(1,100),F(100)),(F(1,10),F(16)),(F(1,10**14),F(16))):
 a=m/(1-m);b=M/(M-1)
 for pi in [F(i,1000) for i in range(1,1000)]+[F(1,10**20),1-F(1,10**20)]:
  for c in (-min((1-pi)*a,pi*b),min(pi*a,(1-pi)*b)):
   assert abs(c)/(pi*(1-pi))<=a+b
   tau_checks+=1
result={'legal_endpoint_tau_bound_exact_cases':tau_checks,'exact_fraction_enumeration_cases':len(caps),'fourier_checks':caps,'n2_periodic_mass_errors':errs,'max_mass_error':max(errs)}
(out/'CHECKS.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({'exact_cases':len(caps),'max_mass_error':max(errs)}))
