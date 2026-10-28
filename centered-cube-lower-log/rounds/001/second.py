from pathlib import Path
import numpy as np,json,time
from kernels import iso,general
P=Path(__file__).resolve().parent; start=time.perf_counter();rows=[]
def threshold(name):
 v=np.load(P/(name+'.npy'));v=np.sort(v[:len(v)//2])[::-1];r=np.arange(1,len(v)+1);ok=r>=np.ceil(.1*len(v));k=np.where(ok)[0][np.argmax((v+np.log(r/len(v)))[ok])];return float(np.nextafter(v[k],-np.inf))
spec=[]
for n in [32,128,512,2048]:spec.append((n,0,0))
for n in [128,512]:
 for si,sg in [(1,.02),(2,.08)]:spec.append((n,si,sg))
for n,si,sg in spec:
 jj=[2,4,8,16,32,64] if sg==0 else [16];names={j:(f'iso_n{n}_J{j}' if sg==0 else f'block_n{n}_s{si}') for j in jj};hh={j:threshold(names[j]) for j in jj};parts={j:[] for j in jj}
 rng=np.random.default_rng(np.random.SeedSequence([2026102802,n,si]));a=np.exp(sg*np.r_[-np.ones(n//2),np.ones(n//2)])
 for _ in range(128):
  x=rng.random((64,n))*.5;out=iso(x,jj) if sg==0 else {16:general(x,a,16)}
  for j in jj:parts[j].append(out[j])
 for j in jj:
  v=np.concatenate(parts[j]);p=float(np.mean(v>hh[j]));lam=float(np.exp(hh[j]));row=dict(name=names[j],n=n,J=j,sigma=sg,N=len(v),threshold=lam,hits=int(np.sum(v>hh[j])),estimate=lam*p,mcse=lam*np.sqrt(p*(1-p)/len(v)),normalized_sqrt=lam*p/np.sqrt(n));rows.append(row);np.save(P/('second_'+names[j]+'.npy'),v)
 print(n,sg,flush=True)
(P/'second_results.json').write_text(json.dumps(dict(seed=2026102802,seconds=time.perf_counter()-start,rows=rows),indent=2));print('done',time.perf_counter()-start)
