from pathlib import Path
import numpy as np,json,time
P=Path(__file__).resolve().parent
src=(P.parent/'017/probe.py').read_text();exec(src.split('rows=[]')[0],globals());P=Path(__file__).resolve().parent
out=[];a=.5/K;n=900000000;h=30.;start=time.perf_counter()
for name,L,theta in [('near',1.16+np.array([0,1,4])/(16*K),np.array([.55,0,.55])),('far',np.array([1.13,1.17,1.21]),np.array([.7,.65,.64]))]:
 edges=np.unique(np.concatenate([np.array([0.,1.])]+[np.mod(atoms+edge+shift,1) for ll in L for edge in [-delta/2,delta/2] for shift in [-ll/2,ll/2]]));mid=(edges[:-1]+edges[1:])/2;half=np.diff(edges)/2
 z,w=np.polynomial.legendre.leggauss(4);x=(mid[:,None]+half[:,None]*z).ravel();weights=(half[:,None]*np.broadcast_to(w,(len(mid),4))).ravel();weights/=weights.sum()
 B=np.cos(2*np.pi*x[:,None]*m)@(np.sqrt(2*b/3)[:,None]*np.sin(np.pi*m[:,None]*L));X=(prim(x[:,None]+L/2)-prim(x[:,None]-L/2)+B)/L;Z=np.log1p(a*X);u=Z@theta;F=float(np.log1p(weights@np.expm1(u)));q=weights*np.exp(u-F);q/=q.sum()
 rng=np.random.default_rng(np.random.SeedSequence([2026102903,name=='far']));ws=[];hit=0
 for batch in range(32):
  C=rng.multinomial(n,q,size=64);S=C@Z;ev=np.all(S>=h,axis=1);hit+=int(ev.sum());ww=np.zeros(64);ww[ev]=np.exp(-(S[ev]-h)@theta);ws.extend(ww)
 ws=np.array(ws);factor=np.exp(n*F-h*theta.sum());p=factor*ws.mean();se=factor*ws.std(ddof=1)/np.sqrt(len(ws));r=dict(name=name,n=n,h=h,theta=theta.tolist(),nodes=len(x),N=len(ws),hits=hit,probability=float(p),se=float(se),relative_se=float(se/p),prefactor=float(ws.mean()),chernoff=float(factor));out.append(r);print(json.dumps(r),flush=True)
(P/'tails_results.json').write_text(json.dumps(dict(seconds=time.perf_counter()-start,rows=out),indent=2))
