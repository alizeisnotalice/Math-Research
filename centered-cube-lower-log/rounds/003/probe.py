from pathlib import Path
import sys,json,time
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'001'))
from kernels import iso
P=Path(__file__).resolve().parent
start=time.perf_counter();summary=[]
# Freeze proposals, dimensions and thresholds before any samples.
for n in [128,512,2048]:
 cuts=np.array([.5]+[c for c in [.1,.25,.5,.75,.9] for s in [-4,-2,2,4]])
 shifts=np.array([0]+[s for c in [.1,.25,.5,.75,.9] for s in [-4,-2,2,4]])
 probs=1/(1+np.exp(-(np.log(cuts/(1-cuts))+shifts/np.sqrt(n*cuts*(1-cuts)))))
 assert np.all((probs>0)&(probs<1))
 mixing=np.array([.25]+[.75/20]*20)
 rng=np.random.default_rng(np.random.SeedSequence([2026102805,n]));vs=[];ws=[]
 for batch in range(512):
  comp=rng.choice(len(cuts),64,p=mixing);left=rng.random((64,n))<probs[comp,None];u=rng.random((64,n))
  x=np.where(left,u*cuts[comp,None],cuts[comp,None]+u*(1-cuts[comp,None]))
  counts=np.stack([(x<c).sum(axis=1) for c in cuts],axis=1)
  logq=counts*np.log(probs/cuts)+(n-counts)*np.log((1-probs)/(1-cuts))+np.log(mixing)
  mx=logq.max(axis=1);logmix=mx+np.log(np.exp(logq-mx[:,None]).sum(axis=1))
  vs.append(iso(x*.5,[16])[16]);ws.append(np.exp(-logmix))
 v=np.concatenate(vs);w=np.concatenate(ws);assert w.max()<=4+1e-9
 hgrid=np.arange(1.,17.);mid=len(v)//2
 train=np.array([np.exp(h)*np.mean(w[:mid]*(v[:mid]>h)) for h in hgrid]);h=float(hgrid[train.argmax()])
 rows=[]
 for hh in hgrid:
  z=w[mid:]*(v[mid:]>hh);rows.append(dict(h=float(hh),estimate=float(np.exp(hh)*z.mean()),se=float(np.exp(hh)*z.std(ddof=1)/np.sqrt(mid)),ess=float(z.sum()**2/(z@z)) if z@z else 0,hits=int((v[mid:]>hh).sum())))
 chosen=next(r for r in rows if r['h']==h)
 result=dict(n=n,N=len(v),weight_mean=float(w.mean()),weight_se=float(w.std(ddof=1)/np.sqrt(len(w))),max_weight=float(w.max()),frozen_h=h,chosen=chosen,rows=rows)
 summary.append(result);np.savez(P/f'raw_n{n}.npz',response=v,weight=w)
 print(json.dumps({k:val for k,val in result.items() if k!='rows'}),flush=True)
(P/'results.json').write_text(json.dumps(dict(seconds=time.perf_counter()-start,seed=2026102805,results=summary),indent=2))
