import numpy as np, scipy, json,time,platform,warnings
from scipy.optimize import linprog,differential_evolution,OptimizeWarning
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
warnings.filterwarnings('ignore',message='Unrecognized options detected.*',category=OptimizeWarning)
def capacity(C,ell,a):
 R=np.c_[C[:,0]-ell/2,C[:,0]+ell/2,C[:,1]-ell/2,C[:,1]+ell/2];X=np.unique(R[:,:2]);Y=np.unique(R[:,2:]);xx,yy=np.meshgrid((X[:-1]+X[1:])/2,(Y[:-1]+Y[1:])/2,indexing='ij');ar=np.diff(X)[:,None]*np.diff(Y)[None,:]
 inside=(xx.ravel()[None,:]>R[:,0,None])&(xx.ravel()[None,:]<R[:,1,None])&(yy.ravel()[None,:]>R[:,2,None])&(yy.ravel()[None,:]<R[:,3,None]);mask=(inside*(1<<np.arange(8))[:,None]).sum(0);cap=np.bincount(mask,weights=ar.ravel(),minlength=256);ids=np.where(cap[1:]>1e-18)[0]+1;c=cap[ids];B=((ids[None,:]>>np.arange(8)[:,None])&1).astype(float);n=len(c)
 res=linprog(np.r_[np.zeros(n),-1],A_ub=np.vstack([np.r_[np.ones(n),0],np.c_[-B,np.ones(8)]]),b_ub=np.r_[1,-a*ell],bounds=[(0,v) for v in c]+[(None,None)],method='highs',options={'threads':1,'primal_feasibility_tolerance':1e-10,'dual_feasibility_tolerance':1e-10});assert res.success
 return res.x[-1]
def run(job):
 j,t,sigma=job;a=1-t;A=(1/a-a)/2*(1-sigma*t);B=(a-np.sqrt(2*a*a-1))/2*(1-sigma*t);C=np.array([[A,0],[-A,0],[0,A],[0,-A],[B,B],[-B,B],[B,-B],[-B,-B]])
 def scales(v):
  rho,m,f=v;lo=-1+abs(rho);hi=1-abs(rho);ld=lo+(hi-lo)*f+t*(2*f-1)
  lam=np.array([-rho,-rho,rho+t*m,rho+t*m,ld,ld,ld,ld]);return np.clip(1+t*lam,a+1e-10,1/a-1e-10)
 def obj(v):return -capacity(C,scales(v),a)/t**2
 res=differential_evolution(obj,[(-1,1),(-4,4),(0,1)],seed=20261105+j,popsize=16,maxiter=220,polish=False,tol=1e-6,atol=1e-7)
 return {'t':t,'sigma':sigma,'params':res.x.tolist(),'scales':scales(res.x).tolist(),'centers':C.tolist(),'eta_over_t2':-res.fun,'evaluations':res.nfev,'success':bool(res.success),'seed':20261105+j}
if __name__=='__main__':
 jobs=[(j,t,s) for j,(t,s) in enumerate(( (t,s) for t in [.02,.01,.005,.001] for s in [1/3,.36,.375,.39,.42]))];st=time.monotonic()
 with ProcessPoolExecutor(max_workers=8) as pool:rows=list(pool.map(run,jobs))
 out={'rows':rows,'workers':8,'seconds':time.monotonic()-st,'evaluations':sum(r['evaluations'] for r in rows),'python':platform.python_version(),'numpy':np.__version__,'scipy':scipy.__version__,'scope':'reflection-symmetric scale search only; negative values do not certify infeasibility'};(Path(__file__).resolve().parents[1]/'results/eight_seeded_search.json').write_text(json.dumps(out,indent=2));print([(r['t'],r['sigma'],r['eta_over_t2'],r['params']) for r in rows]);print(out['seconds'],out['evaluations'])
