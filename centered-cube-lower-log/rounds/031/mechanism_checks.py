"""Exact finite weak-L1 fixtures and periodic truncation arithmetic; not cube lower bounds."""
from pathlib import Path
from fractions import Fraction as F
import numpy as np
import json

rows=[]
for R in [2,4,8,16,32,64,128,256,1024]:
    x=np.arange(R)
    template=R/(np.arange(R)+1.)
    means=np.zeros(R)
    for i in range(R):
        means+=template[(x+i)%R]/R
    H=sum((F(1,k) for k in range(1,R+1)),F(0))
    assert np.max(abs(means-float(H)))<1e-11
    marginal=np.max(np.sort(template)[::-1]*(np.arange(R)+1)/R)
    assert abs(marginal-1)<1e-14
    rows.append(dict(R=R,marginal_weak_coefficient=float(marginal),
                     mean_weak_coefficient=float(means[0]),harmonic_exact=str(H),
                     ratio_to_logR=float(H)/np.log(R)))
boundary=[]
b=F(11,9)
for n in [1,8,32,64,128]:
    N=1000*n
    interior=(1-b/(2*N))**n
    loss=1-interior
    bound=n*b/(2*N)
    assert 0<=loss<=bound
    cell_core=(1-F(1,N))**n
    assert cell_core>=1-F(n,N)
    boundary.append(dict(n=n,N=N,interior_fraction=float(interior),
                         boundary_loss=float(loss),union_bound=str(bound),
                         integer_cell_core_fraction=float(cell_core)))
out=dict(numpy_version=np.__version__,cyclic_fixture=rows,truncation=boundary)
(Path(__file__).resolve().parent/'mechanism_checks.json').write_text(json.dumps(out,indent=2))
print(json.dumps({'R':[x['R'] for x in rows],
                  'harmonic':[x['mean_weak_coefficient'] for x in rows],
                  'boundary_checks':len(boundary)},indent=2))
