"""Recompute 16 actual Sobol points with independent direct intersection and continuous minimization.
Floating-point numerical witness checks; no certified integration or global volume bound.
"""
from pathlib import Path
import itertools,json
import numpy as np
from scipy.stats import qmc
from scipy.optimize import minimize_scalar
base=Path(__file__).resolve().parent; data=json.loads((base/'rescramble'/'selected.json').read_text());n=2;k=data['k'];u=data['u'];cs=np.array(data['centers']);ws=np.array(data['halfwidths']);site=0
pts=cs[site]+(2*qmc.Sobol(2,scramble=True,seed=data['seed']*100+site).random_base2(14)[:16]-1)*ws[site]
# Deliberately enumerate all nearby period shifts, without the search code's prefilter.
boxes=[(c+k*np.array(j),w) for c,w in zip(cs,ws) for j in itertools.product(range(-2,3),repeat=2)]
def val(x,r):
 volume=0.
 for c,w in boxes:
  h=np.maximum(0.,np.minimum(x+r,c+w)-np.maximum(x-r,c-w))
  volume+=float(h[0]*h[1])
 return volume/(2*r)**2
rows=[]
for x in pts:
 breaks={1e-8,1.5}
 for c,w in boxes:
  for z,d in zip(np.abs(x-c),w):
   breaks.update(r for r in (abs(z-d),z+d) if 1e-8<r<1.5)
 breaks=sorted(breaks);best=(1.,0.)
 for lo,hi in zip(breaks[:-1],breaks[1:]):
  for r in (lo,hi):
   f=val(x,r)
   if f<best[0]:best=(f,r)
  z=minimize_scalar(lambda r:val(x,r),bounds=(lo,hi),method='bounded',options={'xatol':1e-12})
  if z.fun<best[0]:best=(float(z.fun),float(z.x))
 rows.append({'point':x.tolist(),'defect_fraction':best[0],'radius':best[1],'strict_margin':u*u-best[0],'strictly_accepted':bool(best[0]<u*u)})
result={'scope':'floating-point direct sum, all radius-event intervals in [1e-8,1.5]; scipy numerical minimization, not certified', 'u':u,'site':site,'examples':rows,'all_strict':all(r['strictly_accepted'] for r in rows),'minimum_margin':min(r['strict_margin'] for r in rows)}
(base/'rescramble'/'witnesses.json').write_text(json.dumps(result,indent=2));print(json.dumps({'all_strict':result['all_strict'],'minimum_margin':result['minimum_margin'],'count':len(rows)}))
