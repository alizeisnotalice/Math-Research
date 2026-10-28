from pathlib import Path
import numpy as np,json
P=Path(__file__).resolve().parent;R=32;p=.5;rows=[]
for N in [16,64,256,1024]:
 for kind in ['single','adjacent']:
  freq={}
  for r,weight in ([(R,1.)] if kind=='single' else [(R,.5),(R+1,.5)]):
   for j in range(1,N):freq[j*r]=freq.get(j*r,0)+2*p*weight*(1-j/N)
  m=np.array(sorted(freq));coef=np.array([freq[k] for k in m])
  L=np.linspace(10/9,11/9,1025)
  def features(l):return coef*np.sin(np.pi*l[:,None]*m)/(np.pi*l[:,None]*m)/np.sqrt(2)
  A=features(L);v=np.sum(A*A,axis=1)*R*R
  increments=[]
  for scale in [1,4,16,64,256]:
   delta=1/(R*scale);diff=features(L+delta)-A
   increments.append(dict(scale=scale,delta=delta,scaled_increment_max=float(np.max(np.sum(diff*diff,axis=1))*R*R)))
  rows.append(dict(N=N,kind=kind,scaled_variance_min=float(v.min()),scaled_variance_max=float(v.max()),increments=increments))
(P/'results.json').write_text(json.dumps(rows,indent=2))
for r in rows:print(r['N'],r['kind'],r['scaled_variance_max'],[(q['scale'],round(q['scaled_increment_max'],7)) for q in r['increments']])
