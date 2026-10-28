from pathlib import Path
import numpy as np,json,itertools
P=Path(__file__).resolve().parent
# Reuse only deterministic source definitions, without running round017 output code.
src=(P.parent/'017/probe.py').read_text();exec(src.split('rows=[]')[0],globals());P=Path(__file__).resolve().parent
records=[]
for G in [32768,131072]:
 x=(np.arange(G)+.5)/G
 for name,L in [('near',1.16+np.array([0,1,4])/(16*K)),('far',np.array([1.13,1.17,1.21]))]:
  B=np.cos(2*np.pi*x[:,None]*m)@(np.sqrt(2*b/3)[:,None]*np.sin(np.pi*m[:,None]*L))
  X=(prim(x[:,None]+L/2)-prim(x[:,None]-L/2)+B)/L;X-=X.mean(axis=0);C=X.T@X/G;target=np.diag(C)/2
  for a in [.0004,.0001]:
   Z=np.log1p(a*X);faces=[]
   def evaluate(theta):
    u=Z@theta;mx=u.max();ww=np.exp(u-mx);ww/=ww.sum();mean=ww@Z
    F=(mx+np.log(np.exp(u-mx).mean()))/a**2-target@theta
    grad=mean/a**2-target;Y=Z-mean;H=(Y.T*ww)@Y/a**2
    return F,grad,H
   for mask in range(8):
    ids=np.array([i for i in range(3) if mask>>i&1],dtype=int);theta=np.zeros(3)
    if len(ids):
     theta[ids]=np.linalg.solve(C[np.ix_(ids,ids)],np.diag(C)[ids])
     for step in range(15):
      F,grad,H=evaluate(theta)
      if abs(grad[ids]).max()<2e-10:break
      change=np.linalg.solve(H[np.ix_(ids,ids)],grad[ids]);theta[ids]-=change
    F,grad,H=evaluate(theta);inactive=[i for i in range(3) if i not in ids]
    if np.all(theta>=-1e-8) and (not inactive or min(grad[inactive])>=-1e-8):faces.append((F,theta,grad,H,mask))
   assert len(faces)==1,len(faces)
   F,theta,grad,H,mask=faces[0]
   row=dict(G=G,name=name,a=a,mask=mask,theta=theta.tolist(),gradient=grad.tolist(),objective=F,tilted_hessian_eigenvalues=np.linalg.eigvalsh(H).tolist())
   records.append(row);print(G,name,a,theta,grad,flush=True)
(P/'nonlinear_results.json').write_text(json.dumps(records,indent=2))
