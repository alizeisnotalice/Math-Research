import os
os.environ['OPENBLAS_NUM_THREADS']='1';os.environ['OMP_NUM_THREADS']='1';os.environ['VECLIB_MAXIMUM_THREADS']='1'
from pathlib import Path
import numpy as np,json,time
P=Path(__file__).resolve().parent
src=(P.parent/'017/probe.py').read_text();exec(src.split('rows=[]')[0],globals());P=Path(__file__).resolve().parent
import argparse
parser=argparse.ArgumentParser();parser.add_argument('--scales',type=int,choices=[65,129],default=65);parser.add_argument('--replications',type=int,default=256);parser.add_argument('--seed',type=int,default=2026102904);parser.add_argument('--output',default='results.json');args=parser.parse_args()
assert args.replications>1
start=time.perf_counter();a=.5/K;n=900000000;h=30.;r=args.scales;L=np.linspace(10/9,11/9,r)
edges=np.unique(np.concatenate([np.array([0.,1.])]+[np.mod(atoms+edge+shift,1) for ll in L for edge in [-delta/2,delta/2] for shift in [-ll/2,ll/2]]));mid=(edges[:-1]+edges[1:])/2;half=np.diff(edges)/2
# Order2 intentionally tests union mechanics first; no continuum error certification.
z,w=np.polynomial.legendre.leggauss(2);x=(mid[:,None]+half[:,None]*z).ravel();base=(half[:,None]*np.broadcast_to(w,(len(mid),2))).ravel();base/=base.sum()
print('nodes',len(x),'precompute',flush=True)
Z=np.empty((len(x),r));means=[]
for j,ll in enumerate(L):
 B=np.cos(2*np.pi*x[:,None]*m)@(np.sqrt(2*b/3)*np.sin(np.pi*m*ll));X=(prim(x+ll/2)-prim(x-ll/2)+B)/ll
 means.append(float(base@X));Z[:,j]=np.log1p(a*X)
F=np.log1p(a*np.array(means));print('nF range',n*F.min(),n*F.max(),'seconds',time.perf_counter()-start,flush=True)
# Bounded initial budget 256 reps, frozen before any outcome. Winner chosen per vector.
rng=np.random.default_rng(args.seed);values={k:[] for k in ([3,9,17,33,65]+([129] if r==129 else []))};hits={k:0 for k in values}
for rep in range(args.replications):
 winner=int(rng.integers(r));q=base*np.exp(Z[:,winner]-F[winner]);q/=q.sum();C=rng.multinomial(n,q);S=C@Z
 loglik=S-n*F;mx=loglik.max();logD=mx+np.log(np.exp(loglik-mx).mean());weight=np.exp(h-logD)
 for k in values:
  ids=np.linspace(0,r-1,k).astype(int);event=bool(S[ids].max()>h);hits[k]+=event;values[k].append(weight if event else 0.)
 if (rep+1)%64==0:print('rep',rep+1,'elapsed',time.perf_counter()-start,flush=True)
rows=[]
for k,v in values.items():
 v=np.array(v);rows.append(dict(scales=k,N=len(v),hits=hits[k],weak_ratio=float(v.mean()),se=float(v.std(ddof=1)/np.sqrt(len(v))),max_weight=float(v.max())))
(P/args.output).write_text(json.dumps(dict(seconds=time.perf_counter()-start,n=n,h=h,a=a,source_K=K,scales=L.tolist(),quadrature_order=2,nodes=len(x),nF_range=[float(n*F.min()),float(n*F.max())],rows=rows),indent=2));print(json.dumps(rows,indent=2))
