"""Exact conditional-binomial variance checks; no Gaussian proxy."""
from pathlib import Path
from fractions import Fraction as F
import math,json,hashlib
HERE=Path(__file__).resolve().parent
cases=0;maxratio=F(0);minpositive=None
for n in [2,4,8,16,64,128]:
 for pi in [F(1,10),F(3,10),F(1,2),F(9,10)]:
  pmf=[F(math.comb(n,k))*pi**k*(1-pi)**(n-k) for k in range(n+1)]
  assert sum(pmf)==1
  for reverse in [False,True]:
   a=b=c=F(0)
   for k in (range(n,-1,-1) if reverse else range(n+1)):
    a+=pmf[k];b+=k*pmf[k];c+=k*k*pmf[k]
    var=c/a-(b/a)**2
    assert 0<=var<=F(n,4)
    maxratio=max(maxratio,var/n);cases+=1
    if a>0:minpositive=a if minpositive is None else min(minpositive,a)
# Rare endpoint gate under exactlypi1/2 is kept as an exact fraction.
rare=F(1,2**128)
out=dict(scope='Fraction checks of conditional Binomial variance; theorem separately proved.',
         tail_cases=cases,max_conditional_variance_over_n=str(maxratio),
         smallest_tested_nonzero_tail_mass=str(minpositive),n128_endpoint_gate_mass=str(rare),
         n128_endpoint_gate_mass_float=float(rare),script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
(HERE/'EXACT_CHECKS.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps({'tail_cases':cases,'max_var_over_n':str(maxratio),'n128_endpoint_mass':str(rare)}))
