from pathlib import Path
import numpy as np,json,itertools,time
P=Path(__file__).resolve().parent
src=(P.parent/'017/probe.py').read_text();exec(src.split('rows=[]')[0],globals());P=Path(__file__).resolve().parent
start=time.perf_counter();records=[];a=.5/K
for G in [8192,32768]:
 x=(np.arange(G)+.5)/G
 for name,L in [('near',1.16+np.array([0,1,4])/(16*K)),('far',np.array([1.13,1.17,1.21]))]:
  B=np.cos(2*np.pi*x[:,None]*m)@(np.sqrt(2*b/3)[:,None]*np.sin(np.pi*m[:,None]*L))
  X=(prim(x[:,None]+L/2)-prim(x[:,None]-L/2)+B)/L;mean=X.mean(axis=0);X-=mean;C=X.T@X/G;Z=np.log1p(a*X)
  for v in [16.,32.,64.]:
   n=int(round(v/(a*a*np.diag(C).mean())));h=v/2-2
   def tilted(th):
    u=Z@th;mx=u.max();f=mx+np.log(np.exp(u-mx).mean());q=np.exp(u-mx);q/=q.sum();mu=q@Z;Y=Z-mu;H=(Y.T*q)@Y
    return f,q,n*mu-h,n*H
   solutions=[]
   for mask in range(8):
    ids=np.array([j for j in range(3) if mask>>j&1],dtype=int);th=np.zeros(3)
    if len(ids):
     th[ids]=np.linalg.solve(C[np.ix_(ids,ids)],np.diag(C)[ids]/2+h/(n*a*a))
     for step in range(12):
      f,q,g,H=tilted(th)
      if max(abs(g[ids]))<1e-8:break
      th[ids]-=np.linalg.solve(H[np.ix_(ids,ids)],g[ids])
    f,q,g,H=tilted(th);inactive=[j for j in range(3) if j not in ids]
    if min(th)>=-1e-8 and (not inactive or min(g[inactive])>=-1e-8):solutions.append((th,f,q,g))
   assert len(solutions)==1
   th,f,q,g=solutions[0];rng=np.random.default_rng(np.random.SeedSequence([2026102902,G,int(v),name=='far']));weights=[];hits=0
   for batch in range(16):
    counts=rng.multinomial(n,q,size=128);S=counts@Z;event=np.all(S>=h,axis=1);hits+=int(event.sum());w=np.zeros(128);w[event]=np.exp(-(S[event]-h)@th);weights.extend(w)
   w=np.array(weights);pref=n*f-h*th.sum();est=np.exp(pref)*w.mean();se=np.exp(pref)*w.std(ddof=1)/np.sqrt(len(w))
   row=dict(G=G,name=name,v=v,n=n,h=h,theta=th.tolist(),inactive_gradient=g.tolist(),source_mean_grid_error=mean.tolist(),N=len(w),hits=hits,probability=float(est),se=float(se),chernoff=float(np.exp(pref)),prefactor=float(w.mean()),relative_se=float(se/est),weighted_joint=float(np.exp(h)*est))
   records.append(row);print(json.dumps(row),flush=True)
(P/'centered_results.json').write_text(json.dumps(dict(seconds=time.perf_counter()-start,rows=records),indent=2))
