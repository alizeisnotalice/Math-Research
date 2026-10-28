from pathlib import Path
import numpy as np,json,time
P=Path(__file__).resolve().parent
src=(P.parent/'017/probe.py').read_text();exec(src.split('rows=[]')[0],globals());P=Path(__file__).resolve().parent
rows=[];a=.5/K;n=900000000;start=time.perf_counter()
for name,L,theta in [('near',1.16+np.array([0,1,4])/(16*K),np.array([.55,0,.55])),('far',np.array([1.13,1.17,1.21]),np.array([.7,.65,.64]))]:
 breaks=[np.array([0.,1.])]
 for length in L:
  for edge in [-delta/2,delta/2]:
   for shift in [-length/2,length/2]:breaks.append(np.mod(atoms+edge+shift,1))
 edges=np.unique(np.concatenate(breaks));left=edges[:-1];right=edges[1:];mid=(left+right)/2;half=(right-left)/2
 for order in [2,4,8,12]:
  z,w=np.polynomial.legendre.leggauss(order);x=(mid[:,None]+half[:,None]*z).ravel();weights=(half[:,None]*np.broadcast_to(w,(len(mid),order))).ravel();weights/=weights.sum()
  B=np.cos(2*np.pi*x[:,None]*m)@(np.sqrt(2*b/3)[:,None]*np.sin(np.pi*m[:,None]*L))
  X=(prim(x[:,None]+L/2)-prim(x[:,None]-L/2)+B)/L
  mean=weights@X;C=(X.T*weights)@X-np.outer(mean,mean);Z=np.log1p(a*X);u=Z@theta
  F=float(np.log1p(weights@np.expm1(u)));q=weights*np.exp(u-F);q/=q.sum();mu=q@Z
  row=dict(name=name,order=order,segments=len(mid),nodes=len(x),weight_sum=float(weights.sum()),min_weight=float(weights.min()),mean=mean.tolist(),covariance=C.tolist(),F=F,nF=n*F,n_gradient=(n*mu).tolist())
  rows.append(row);print(name,order,len(x),max(abs(mean)),row['nF'],row['n_gradient'],flush=True)
(P/'results.json').write_text(json.dumps(dict(n=n,theta_near=[.55,0,.55],theta_far=[.7,.65,.64],seconds=time.perf_counter()-start,rows=rows),indent=2))
