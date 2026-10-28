from pathlib import Path
import json,numpy as np,time
P=Path(__file__).resolve().parent;start=time.perf_counter();rows=[]
G=2048;N=32;p=.9;L=np.linspace(1.02,1.48,257);x=(np.arange(G)+.5)/G
for n in [128,512,2048]:
 R=max(1,round(np.sqrt(n)/4));rng=np.random.default_rng(np.random.SeedSequence([2026102806,n]));C=rng.multinomial(n,np.full(G,1/G),size=4096).astype(float)
 configs=[('single',[R],[1]),('adjacent',[R,R+1],[.5,.5]),('octaves',[R,2*R,4*R],[1/3]*3),('adjacent4',list(range(R,R+4)),[.25]*4)]
 for name,periods,weights in configs:
  freq={}
  for r,w in zip(periods,weights):
   for j in range(1,N):freq[j*r]=freq.get(j*r,0)+2*p*w*(1-j/N)
  m=np.array(sorted(freq));coef=np.array([freq[k] for k in m]);aver=1+(np.sin(np.pi*L[:,None]*m)/(np.pi*L[:,None]*m)*coef)@np.cos(2*np.pi*m[:,None]*x)
  assert aver.min()>0
  scores=C@np.log(aver).T;v=scores.max(axis=1);tr=np.sort(v[:2048])[::-1];k=np.arange(1,2049);ok=k>=103;pos=np.where(ok)[0][np.argmax((tr+np.log(k/2048))[ok])];h=float(np.nextafter(tr[pos],-np.inf));hit=v[2048:]>h;q=hit.mean();lam=np.exp(h)
  row=dict(n=n,R=R,name=name,periods=periods,weights=weights,log_threshold=h,hits=int(hit.sum()),estimate=float(lam*q),se=float(lam*np.sqrt(q*(1-q)/2048)),normalized=float(lam*q/np.sqrt(n)),min_average=float(aver.min()))
  rows.append(row);np.save(P/f'response_{n}_{name}.npy',v);print(json.dumps(row),flush=True)
(P/'results.json').write_text(json.dumps(dict(seconds=time.perf_counter()-start,center_grid=G,scale_grid=len(L),N=N,p=p,seed=2026102806,rows=rows),indent=2))
