"""Frozen mixture importance sampling; common-radius periodic lattice response."""
from pathlib import Path
import numpy as np,json
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'001'))
from kernels import iso
P=Path(__file__).resolve().parent
rows=[]
# Equal mixture: original uniform law and three Bernoulli-count tilts.
# Cut at U=1/4; conditional distribution within each half stays uniform.
for n in [128,512,2048]:
 rng=np.random.default_rng(np.random.SeedSequence([2026102804,n]))
 probs=.5+np.array([0,1,2,3])/np.sqrt(n)
 vals=[];weights=[]
 for batch in range(256):
  component=rng.integers(4,size=64)
  mask=rng.random((64,n))<probs[component,None]
  x=(rng.random((64,n))+(~mask))*.25
  count=mask.sum(axis=1)
  logq=count[:,None]*np.log(2*probs)+(n-count[:,None])*np.log(2*(1-probs))
  mx=logq.max(axis=1)
  logmix=mx+np.log(np.exp(logq-mx[:,None]).mean(axis=1))
  vals.append(iso(x,[16])[16]);weights.append(np.exp(-logmix))
 v=np.concatenate(vals);w=np.concatenate(weights)
 np.savez(P/f'importance_n{n}.npz',log_response=v,weight=w)
 for h in np.arange(1.,17.):
  z=w*(v>h);est=np.exp(h)*z.mean();se=np.exp(h)*z.std(ddof=1)/np.sqrt(len(z))
  rows.append(dict(n=n,log_threshold=h,hits=int((v>h).sum()),estimate=est,mcse=se,normalized=est/np.sqrt(n),tail_ess=float(z.sum()**2/(z@z)) if z@z else 0))
 print(n,'Ew',w.mean(),'SE',w.std(ddof=1)/np.sqrt(len(w)),flush=True)
(P/'importance_results.json').write_text(json.dumps(rows,indent=2))
for n in [128,512,2048]:
 r=[r for r in rows if r['n']==n and r['tail_ess']>=100]
 print('exploratory best (selection biased)',max(r,key=lambda r:r['estimate']))
