from pathlib import Path
import json,numpy as np
P=Path(__file__).resolve().parent
base=json.loads((P.parent/'005/refinement.json').read_text())[-1]
b=np.array(base['coefficients']);m=np.arange(1,len(b)+1);gamma=27/29
# Midpoint quadrature over phase; enumerate exactly admissible nonnegative integers.
x=(np.arange(200000)+.5)/200000
f=np.empty_like(x)
for ids in np.array_split(np.arange(len(x)),200):f[ids]=np.sin(np.pi*x[ids,None]*m)**2@b
rows=[]
for scale in [1,16,256,4096]:
 F=1+scale*f
 for v in [16,64,256,1024]:
  e=1/v;lo=np.sqrt(F/(gamma*(1+e)));hi=np.sqrt(F/(gamma*(1-e)))
  counts=np.maximum(0,np.floor(hi-x)-np.maximum(0,np.ceil(lo-x))+1)
  delta=1/np.sqrt(1-e)-1/np.sqrt(1+e)
  bound=1+delta*np.sqrt((1+scale*b.sum()/2)/gamma)
  rows.append(dict(scale=scale,v=v,grid_points=len(x),measured_length=float(counts.mean()),max_count=int(counts.max()),bound=float(bound),mass=float(scale*b.sum())))
(P/'results.json').write_text(json.dumps(rows,indent=2));print(json.dumps(rows,indent=2))
