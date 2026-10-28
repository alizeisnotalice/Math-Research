import numpy as np,json,itertools
from pathlib import Path
# Independent exhaustive quadrature on proposal density's constant cells.
rows=[]
for n in [1,2,3]:
 cuts=np.array([.5]+[c for c in [.1,.25,.5,.75,.9] for s in [-4,-2,2,4]])
 shifts=np.array([0]+[s for c in [.1,.25,.5,.75,.9] for s in [-4,-2,2,4]])
 probs=1/(1+np.exp(-(np.log(cuts/(1-cuts))+shifts/np.sqrt(n*cuts*(1-cuts)))))
 mix=np.array([.25]+[.75/20]*20)
 edges=np.array([0,.1,.25,.5,.75,.9,1]);width=np.diff(edges);mid=(edges[:-1]+edges[1:])/2
 mass=weighted=0.;wmax=0.
 for ids in itertools.product(range(6),repeat=n):
  x=mid[list(ids)];vol=np.prod(width[list(ids)])
  density=0.
  for c,p,m in zip(cuts,probs,mix):
   density+=m*np.prod(np.where(x<c,p/c,(1-p)/(1-c)))
  mass+=vol*density;weighted+=vol*density/density;wmax=max(wmax,1/density)
 rows.append(dict(n=n,proposal_integral=mass,weighted_integral=weighted,max_weight=wmax))
 assert abs(mass-1)<1e-12 and abs(weighted-1)<1e-12 and wmax<=4
Path(__file__).with_name('control.json').write_text(json.dumps(rows,indent=2));print(rows)
