"""Independent interior second-arrival implementation; no boundary conditioning in MC."""
import json
import math
from pathlib import Path
import numpy as np
from pressure import tensor_maximal

rng=np.random.default_rng(660811)
cases=[]
for n in [8,32,128]:
    for j in range(200):
        r=float(rng.uniform(1.01,1.99))
        dx=r*rng.uniform(-.5,.5,size=n)
        dx[j % n]=r/2
        z=rng.integers(1,256,size=n)
        # At R=1 each generic interior coordinate contains one integer;
        # its second integer arrives at 2 max(|dx_i|,1-|dx_i|).
        S=np.sort(2*np.maximum(np.abs(dx),1-np.abs(dx)))
        scores=np.arange(1,n+1)*math.log(2)-n*np.log(S)
        last=np.ones(n,dtype=bool);last[:-1]=S[:-1]!=S[1:]
        scores=np.where(last,scores,-np.inf)
        vals=np.r_[0.,scores,0.]
        radii=np.r_[1.,S,2.]
        a=int(np.argmax(vals))
        v,R,_,_=tensor_maximal(dx,z,257)
        assert abs(v-vals[a])<4e-12 and abs(R-radii[a])<2e-12
        cases.append({'n':n,'max_logvalue_abs_error':float(abs(v-vals[a])),'winner_abs_error':float(abs(R-radii[a]))})
result={'status':'PASS_INTERIOR_SECOND_ARRIVAL_FORMULA','cases':len(cases),'max_logvalue_abs_error':max(c['max_logvalue_abs_error'] for c in cases),'max_winner_abs_error':max(c['winner_abs_error'] for c in cases),'scope':'Only generic interior z_i in [1,255], r in (1.01,1.99). Main MC retains uniform complete lattice including boundary sources.','cdf':'F_r(s)=2(s-1)/r for 1<=s<=r, (s+r-2)/r for r<=s<=2; forced face S=r. Formula CDF is analytical; this check validates same-point arrival/max, not a stochastic CDF test.'}
Path(__file__).with_name('internal_formula_validation.json').write_text(json.dumps(result,indent=2))
print(json.dumps(result))
