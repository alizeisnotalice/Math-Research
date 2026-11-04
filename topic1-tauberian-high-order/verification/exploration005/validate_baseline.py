"""Independent analytic baseline for the radius-event exploration.
Not an error certificate for QMC volume estimates or perturbed inputs.
"""
import importlib.util,itertools,json
from pathlib import Path
import numpy as np
base=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('radius_module',base/'lattice_batch.py')
mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
rng=np.random.default_rng(220005)
maxerr=0.;rows=[]
for u in (.004,.006,.007):
    d=u/(1+2*u)
    boxes=[(np.array(ij,dtype=float),np.array([d,d])) for ij in itertools.product(range(-2,3),repeat=2)]
    pts=[np.array([0.,0.]),np.array([d,0.]),np.array([d,d])]
    pts+=list(rng.uniform(-d,d,size=(24,2)))
    for x in pts:
        expected=(d/(1-d-np.max(np.abs(x))))**2
        got=mod.fit(x,boxes)
        err=abs(got-expected);maxerr=max(maxerr,err)
        assert err<1e-11,(u,x,got,expected)
    rows.append({'u':u,'delta':d,'points':len(pts)})
print(json.dumps({'tests':rows,'pass':True,'max_absolute_error':maxerr,
 'scope':'radius search <=1.5, uniform grid analytic exact formula; floating implementation check only'},indent=2))
