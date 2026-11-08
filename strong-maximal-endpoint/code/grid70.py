import os
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['OMP_NUM_THREADS']='1'
import sys,json,time,math,itertools
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor
import numpy as np
from scipy.integrate import quad
from scipy.optimize import milp,LinearConstraint,Bounds
from scipy.sparse import csr_matrix,hstack,vstack,eye
from scipy.stats import irwinhall
from scipy.special import gammaln
MU=1.2784645427610738

def run(j):
 start=time.monotonic();n=j['N'];m=j['source_side'];r=j['receiver_outer'];step=j['step'];mode=j['mode'];budget=j['budget'];c=MU if mode=='finite' else None
 src=np.array(list(itertools.product(range(m),repeat=n)),dtype=int)
 rec=np.array(list(itertools.product(range(-r,m),repeat=n)),dtype=int)
 ss=src.sum(1);rs=rec.sum(1);sw=np.exp(-step*ss)*(-math.expm1(-step))**n;rw=np.exp(-step*rs)*(-math.expm1(-step))**n
 top=int(n*(m+r));fact=math.factorial(n-1);binoms=np.array([math.comb(n,k) for k in range(n+1)])
 def H(s):
  if c is None:return s**(n-1)/fact
  return float(binoms@np.maximum(s-c*np.arange(n+1),0)**(n-1))/fact
 def baseprob(M):
  # The bin thresholds are aligned with c in finite mode.
  def f(x):
   s=step*(M+x)
   if s==0:return step
   return step**n*math.exp(-s)/H(s)*float(irwinhall.pdf(x,n))
  return sum(quad(f,a,a+1,epsabs=2e-12,epsrel=2e-10,limit=80)[0] for a in range(n))
 probs=np.array([baseprob(M) for M in range(top+1)])
 ii=[];jj=[];vv=[]
 mc=round(c/step) if c is not None else 0
 for a,y in enumerate(rec):
  d=src-y;ok=np.all(d>=0,axis=1);ids=np.flatnonzero(ok);ds=d[ids];v=probs[ds.sum(1)]
  if c is not None:v=v*2.**np.sum(ds>=mc,axis=1)
  vals=budget*v/sw[ids]
  keep=vals>1e-14
  ii.extend([a]*int(keep.sum()));jj.extend(ids[keep]);vv.extend(vals[keep])
 P=csr_matrix((vv,(ii,jj)),shape=(len(rec),len(src)))
 # Source variables are fractions of the unit mass budget; selected receiver
 # cells are binary. This is collocation, not a certified continuum lower bound.
 A=hstack([P,-eye(len(rec),format='csr')],format='csr')
 mass=csr_matrix(np.r_[np.ones(len(src)),np.zeros(len(rec))][None,:])
 cons=[LinearConstraint(A,0,np.inf),LinearConstraint(mass,0,1)]
 obj=np.r_[np.zeros(len(src)),-rw/budget]
 integrality=np.r_[np.zeros(len(src),dtype=int),np.ones(len(rec),dtype=int)]
 result=milp(obj,integrality=integrality,bounds=Bounds(np.zeros(len(obj)),np.ones(len(obj))),constraints=cons,options={'time_limit':8.,'mip_rel_gap':.01,'node_limit':2000})
 ans={'spec':j,'status':int(result.status),'message':result.message,'seconds':time.monotonic()-start,'probabilities':probs.tolist(),'source_cells':len(src),'receiver_cells':len(rec),'scope':'finite collocation search only; no continuum upper/lower certificate'}
 if result.x is not None:
  massx=result.x[:len(src)];chosen=result.x[len(src):]>.5;cost=float(massx.sum());outputs=P@massx
  ans.update({'collocation_ratio':float(rw[chosen].sum()/budget/max(cost,1e-300)),'threshold_residual':float(np.min(outputs[chosen]-1)) if np.any(chosen) else None,'mass_budget':cost,'selected_cells':int(chosen.sum()),'solver_gap':float(getattr(result,'mip_gap',float('nan'))),'source_mass':massx.tolist(),'selected_indices':np.flatnonzero(chosen).tolist()})
 return ans
if __name__=='__main__':
 inp=Path(sys.argv[1]);out=Path(sys.argv[2]);out.mkdir(exist_ok=True,parents=True);jobs=json.loads(inp.read_text())['grid'];start=time.monotonic()
 with ProcessPoolExecutor(max_workers=12) as pool:rows=list(pool.map(run,jobs))
 (out/'grid70.json').write_text(json.dumps({'seconds':time.monotonic()-start,'workers':12,'rows':rows},indent=2)+'\n');print('Completed',len(rows),'grid searches')
