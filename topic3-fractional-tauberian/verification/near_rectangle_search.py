"""Both holes and additions in four corner boxes; actual second-order flat model."""
import numpy as np,json,time,platform
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
class Model:
 def __init__(self,rho,N=12,K=1.5,nd=64):
  self.N=N;self.K=K;h=2*K/N;self.area=h*h;edges=np.linspace(-K,K,N+1);self.baseline=((edges[:-1,None]>=0)&(edges[None,:-1]>=0)).astype(float)
  cuts=sorted(set([0.,.5,(1-rho)/2,(1+rho)/2,1.]));Z,W=np.polynomial.legendre.leggauss(nd//4);x=np.concatenate([(Z+1)*(b-a)/2+a for a,b in zip(cuts,cuts[1:])]);self.w=np.concatenate([W*(b-a)/2 for a,b in zip(cuts,cuts[1:])]);y=1-x;self.x=x
  lo=np.maximum(rho/2-x,-rho/2-y);hi=np.minimum(rho/2+x,-rho/2+y)
  U=np.concatenate([lo[:,None],hi[:,None],rho/2+x[:,None]-edges,rho/2-x[:,None]-edges,-rho/2+y[:,None]-edges,-rho/2-y[:,None]-edges],axis=1);U=np.clip(U,lo[:,None],hi[:,None]);self.nc=U.shape[1]
  self.base=U*U+U+x[:,None]*y[:,None]+rho*(x[:,None]-.5)+rho*rho/4
  self.TX=[np.clip((edges[None,None,1:]-(rho/2+sgn*x[:,None]-U)[:,:,None])/h,0,1).reshape(-1,N) for sgn in [-1,1]]
  self.TY=[np.clip((edges[None,None,1:]-(-rho/2+sgn*y[:,None]-U)[:,:,None])/h,0,1).reshape(-1,N) for sgn in [-1,1]]
 def evaluate(self,f,update=False):
  g=f-self.baseline;total=g.sum()*self.area;B=0;picks=[]
  for c in range(4):
   G=self.base.copy()-.5*total
   for j in range(4):
    tx=self.TX[(c^j)&1];ty=self.TY[((c^j)>>1)&1];G+=np.sum((tx@g[j])*ty,axis=1).reshape(G.shape)*self.area
   pick=G.argmax(axis=1);picks.append(np.arange(len(self.x))*self.nc+pick);B+=self.w@G[np.arange(len(self.x)),pick]
  if not update:return float(B)
  nf=np.empty_like(f)
  for j in range(4):
   coeff=np.full((self.N,self.N),-2.)
   for c in range(4):
    tx=self.TX[(c^j)&1][picks[c]];ty=self.TY[((c^j)>>1)&1][picks[c]];coeff+=(tx.T*self.w)@ty
   nf[j]=np.where(abs(coeff)<1e-12,f[j],coeff>0)
  return float(B),nf

def run(job):
 rho,seed=job;rng=np.random.default_rng(seed);m=Model(rho);best=-100;iterations=0
 for k in range(100):
  if k==0:f=np.tile(m.baseline,(4,1,1))
  elif k%3==0:f=np.tile(m.baseline,(4,1,1))*(rng.random((4,m.N,m.N))>rng.random())
  else:f=(rng.random((4,m.N,m.N))<rng.random()).astype(float)
  for _ in range(100):
   B,nf=m.evaluate(f,True);iterations+=1
   if np.array_equal(f,nf):break
   f=nf
  B=m.evaluate(f)
  if B>best:best=B;bestf=f.copy()
 return {'rho':rho,'seed':seed,'iterations':iterations,'B':best,'refined':Model(rho,nd=1024).evaluate(bestf),'nonbaseline_cells':np.argwhere(abs(bestf-m.baseline)>1e-12).tolist()}
if __name__=='__main__':
 st=time.monotonic();jobs=[(r,20261105+i) for i,r in enumerate([0.,.1,.25,.4,.5,.6,.75,.9,.95,.99,.999,.9999])]
 with ProcessPoolExecutor(max_workers=6) as p:rows=list(p.map(run,jobs))
 out={'rows':rows,'seconds':time.monotonic()-st,'starts':1200,'iterations':sum(r['iterations'] for r in rows),'python':platform.python_version(),'numpy':np.__version__,'parameters':{'N':12,'K':1.5,'quadrature_budget':64,'refined_budget':1024,'starts_per_rho':100,'iteration_cap':100,'processes':6},'scope':'signed corner perturbations, finite contact scales and split Gaussian exploration'};(Path(__file__).resolve().parents[1]/'results/near_rectangle_search.json').write_text(json.dumps(out,indent=2));print({'seconds':out['seconds'],'rows':rows})
