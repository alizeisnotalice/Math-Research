from eight_seeded_search import capacity
import numpy as np,scipy,json,time,platform
from scipy.optimize import differential_evolution
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

def run(job):
 j,t,sigma=job;al=1-t;aa=.75;bb=.25;xx=1-al*al;r=xx/(al+np.sqrt(al*al-4*aa*bb*xx))*(1-sigma*t);C=r*np.array([[aa,bb],[-aa,bb],[aa,-bb],[-aa,-bb],[bb,aa],[-bb,aa],[bb,-aa],[-bb,-aa]])
 def scales(v):return np.clip(1+t*np.repeat(v,4),al+1e-10,1/al-1e-10)
 def obj(v):return -capacity(C,scales(v),al)/t**2
 res=differential_evolution(obj,[(-1,1/(1-t))]*2,seed=20261105+j,popsize=16,maxiter=180,polish=False,tol=1e-7,atol=1e-7)
 return {'t':t,'sigma':sigma,'lambdas':res.x.tolist(),'scales':scales(res.x).tolist(),'centers':C.tolist(),'eta_over_t2':-res.fun,'evaluations':res.nfev,'seed':20261105+j}
if __name__=='__main__':
 jobs=[(j,t,s) for j,(t,s) in enumerate(( (t,s) for t in [.02,.01,.005,.001] for s in [.1,.2,.25,.3,.375]))];st=time.monotonic()
 with ProcessPoolExecutor(max_workers=8) as pool:rows=list(pool.map(run,jobs))
 out={'rows':rows,'workers':8,'seconds':time.monotonic()-st,'evaluations':sum(r['evaluations'] for r in rows),'python':platform.python_version(),'numpy':np.__version__,'scipy':scipy.__version__,'scope':'two-group scale search; negative values not a global upper certificate'};(Path(__file__).resolve().parents[1]/'results/eight_smooth_search.json').write_text(json.dumps(out,indent=2));print([(r['t'],r['sigma'],r['eta_over_t2'],r['lambdas']) for r in rows]);print(out['seconds'],out['evaluations'])
