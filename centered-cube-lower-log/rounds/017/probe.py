from pathlib import Path
import numpy as np,json
P=Path(__file__).resolve().parent;K=1023;z=1;ex=[]
for _ in range(K):
 ex.append(z);z<<=1
 if z&1024:z^=0x409
ex=np.array(ex);log=np.full(1024,-1);log[ex]=np.arange(K)
def tr(a):
 cur=a.copy();out=a.copy()
 for _ in range(9):cur=ex[2*log[cur]%K];out^=cur
 return out
U=(tr(ex)+2*tr(ex[(log[ex]+log[6])%K])+.5)/4;atoms=(np.arange(K)+U)/K;delta=1/(100*K)
def prim(t):
 y=np.mod(t,1);lo=np.searchsorted(atoms,y-delta/2);hi=np.searchsorted(atoms,y+delta/2)
 return lo+np.where(hi>lo,(y-atoms[np.minimum(lo,K-1)]+delta/2)/delta,0)-K*y
m=np.arange(1,17);b=4*(729/29)*((np.sin(np.pi*m*10/27)/(np.pi*m))**4+3*(np.sin(np.pi*m/9)/(np.pi*m))**4)
rows=[]
for G in [32768,131072]:
 x=(np.arange(G)+.5)/G
 for name,L in [('near',1.16+np.array([0,1,4])/(16*K)),('far',np.array([1.13,1.17,1.21]))]:
  B=np.cos(2*np.pi*x[:,None]*m)@(np.sqrt(2*b/3)[:,None]*np.sin(np.pi*m[:,None]*L))
  X=(prim(x[:,None]+L/2)-prim(x[:,None]-L/2)+B)/L;mean=X.mean(axis=0);Y=X-mean;cov=Y.T@Y/G
  # Centering isolates quadrature mean error; it is recorded, not hidden.
  for a in [.0004,.0002,.0001]:
   theta=np.array([1.,1.,1.]);Z=np.log1p(a*Y)@theta
   exact=float(np.log1p(np.expm1(Z).mean()))
   quadratic=float(a*a*(theta@cov@theta-np.dot(theta,np.diag(cov)))/2)
   rows.append(dict(G=G,name=name,a=a,mean=mean.tolist(),covariance=cov.tolist(),eigenvalues=np.linalg.eigvalsh(cov).tolist(),log_mgf=exact,quadratic=quadratic,error=exact-quadratic,error_over_a3=(exact-quadratic)/a**3))
(P/'results.json').write_text(json.dumps(rows,indent=2))
for r in rows:
 if r['G']==131072:print(r['name'],r['a'],r['error_over_a3'],r['eigenvalues'])
