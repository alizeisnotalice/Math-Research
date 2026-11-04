from eight_seeded_search import capacity
import numpy as np,scipy,json,time,platform
from scipy.optimize import differential_evolution
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
seeded=json.loads((Path(__file__).resolve().parents[1]/'results/eight_seeded_search.json').read_text())['rows']
def run(job):
 j,t,sigma=job;a=1-t;A=(1/a-a)/2*(1-sigma*t);B=(a-np.sqrt(2*a*a-1))/2*(1-sigma*t);C=np.array([[A,0],[-A,0],[0,A],[0,-A],[B,B],[-B,B],[B,-B],[-B,-B]])
 prev=min(seeded,key=lambda r:abs(r['t']-t)+abs(r['sigma']-sigma));rr,mm,ff=prev['params'];lo=-1+abs(rr);hi=1-abs(rr);ld=lo+(hi-lo)*ff+t*(2*ff-1);fnew=(ld-lo+4*t)/(hi-lo+8*t);start=np.array([rr,mm,0,0,fnew,fnew,fnew,fnew]);rng=np.random.default_rng(20261105+j)
 bounds=np.array([[-1,1],[-4,4],[-4,4],[-4,4]]+[[0,1]]*4,float);init=rng.uniform(bounds[:,0],bounds[:,1],size=(96,8));init[:64]=np.clip(start+rng.normal(size=(64,8))*np.array([.02,.4,.2,.2,.02,.02,.02,.02]),bounds[:,0],bounds[:,1]);init[0]=start
 def scales(v):
  rho,m,nx,ny=v[:4];fs=v[4:];lo=-1+abs(rho);hi=1-abs(rho);ds=lo+(hi-lo)*fs+4*t*(2*fs-1);lam=np.r_[-rho+t*nx,-rho-t*nx,rho+t*(m+ny),rho+t*(m-ny),ds];return np.clip(1+t*lam,a+1e-10,1/a-1e-10)
 def obj(v):return -capacity(C,scales(v),a)/t**2
 res=differential_evolution(obj,bounds,init=init,seed=20261105+j,maxiter=320,polish=False,tol=1e-6,atol=1e-7)
 return {'t':t,'sigma':sigma,'params':res.x.tolist(),'scales':scales(res.x).tolist(),'centers':C.tolist(),'eta_over_t2':-res.fun,'evaluations':res.nfev,'success':bool(res.success),'seed':20261105+j}
if __name__=='__main__':
 jobs=[(j,t,s) for j,(t,s) in enumerate(( (t,s) for t in [.01,.005,.001] for s in [.36,.375,.39]))];st=time.monotonic()
 with ProcessPoolExecutor(max_workers=8) as pool:rows=list(pool.map(run,jobs))
 out={'rows':rows,'workers':8,'seconds':time.monotonic()-st,'evaluations':sum(r['evaluations'] for r in rows),'python':platform.python_version(),'numpy':np.__version__,'scipy':scipy.__version__,'scope':'independent-scale finite-box numerical search; not a global infeasibility proof'};(Path(__file__).resolve().parents[1]/'results/eight_full_search.json').write_text(json.dumps(out,indent=2));print([(r['t'],r['sigma'],r['eta_over_t2'],r['params']) for r in rows]);print(out['seconds'],out['evaluations'])
