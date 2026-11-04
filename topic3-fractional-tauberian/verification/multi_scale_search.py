import numpy as np,scipy,json,time,platform,warnings
from scipy.optimize import differential_evolution,linprog,OptimizeWarning
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
warnings.filterwarnings('ignore',message='Unrecognized options detected.*',category=OptimizeWarning)

def capacity(C,ell,alpha):
 R=np.c_[C[:,0]-ell/2,C[:,0]+ell/2,C[:,1]-ell/2,C[:,1]+ell/2];X=np.unique(R[:,:2]);Y=np.unique(R[:,2:]);xx,yy=np.meshgrid((X[:-1]+X[1:])/2,(Y[:-1]+Y[1:])/2,indexing='ij');area=np.diff(X)[:,None]*np.diff(Y)[None,:]
 inside=(xx.ravel()[None,:]>R[:,0,None])&(xx.ravel()[None,:]<R[:,1,None])&(yy.ravel()[None,:]>R[:,2,None])&(yy.ravel()[None,:]<R[:,3,None]);masks=(inside*(1<<np.arange(8))[:,None]).sum(axis=0);cap=np.bincount(masks,weights=area.ravel(),minlength=256);ss=np.where(cap[1:]>1e-15)[0]+1;c=cap[ss];A=((ss[None,:]>>np.arange(8)[:,None])&1).astype(float);m=len(c)
 res=linprog(np.r_[np.zeros(m),-1],A_ub=np.vstack([np.r_[np.ones(m),0],np.c_[-A,np.ones(8)]]),b_ub=np.r_[1,-alpha*ell],bounds=[(0,z) for z in c]+[(None,None)],method='highs',options={'threads':1,'primal_feasibility_tolerance':1e-9,'dual_feasibility_tolerance':1e-9});assert res.success
 return float(res.x[-1])
def run(job):
 i,al,sigma=job;t=1-al;c=2*al*al-1;A=(1/al-al)/2;B=(al-np.sqrt(c))/2;shrink=1-sigma*t;C=shrink*np.array([[A,0],[-A,0],[0,A],[0,-A],[B,B],[-B,B],[B,-B],[-B,-B]]);L=al;U=1/al;count=0
 def obj(f):return -capacity(C,L+(U-L)*f,al)/t**2
 res=differential_evolution(obj,[(1e-6,1-1e-6)]*8,seed=20261105+i,popsize=8,maxiter=140,tol=2e-5,atol=1e-8,polish=False,workers=1)
 return {'alpha':al,'sigma':sigma,'centers':C.tolist(),'scale_fractions':res.x.tolist(),'scales':(L+(U-L)*res.x).tolist(),'eta':-res.fun*t*t,'eta_over_t2':-res.fun,'evaluations':int(res.nfev),'optimizer_success':bool(res.success),'seed':20261105+i}
if __name__=='__main__':
 jobs=[(i,al,s) for i,(al,s) in enumerate(( (al,s) for al in [.9,.95,.98] for s in [0,.2,1/3,.375,.5]))];start=time.monotonic()
 with ProcessPoolExecutor(max_workers=8) as pool:rows=list(pool.map(run,jobs))
 out={'rows':rows,'environment':{'python':platform.python_version(),'numpy':np.__version__,'scipy':scipy.__version__},'workers':8,'seconds':time.monotonic()-start,'evaluations':sum(r['evaluations'] for r in rows),'scope':'finite nonsymmetric independent-scale optimization; numerical candidates only, not a global upper bound'};Path(__file__).with_suffix('.json').write_text(json.dumps(out,indent=2));print([(r['alpha'],r['sigma'],r['eta_over_t2']) for r in rows],out['seconds'])
