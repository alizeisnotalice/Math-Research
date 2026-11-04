"""Final direction scan; each LP has one HiGHS thread. Exploratory only."""
import sys,json,time,platform,scipy
from pathlib import Path
from eight_seeded_search import capacity
import numpy as np
from scipy.optimize import differential_evolution
from concurrent.futures import ProcessPoolExecutor

def run(job):
 j,aa,t,s=job;b=1-aa;al=1-t;x=1-al*al;r=x/(al+np.sqrt(al*al-4*aa*b*x))*(1-s*t);C=r*np.array([[aa,b],[-aa,b],[aa,-b],[-aa,-b],[b,aa],[-b,aa],[b,-aa],[-b,-aa]])
 def scales(v):return np.r_[np.repeat(1+v[0]*t,4),np.repeat(1+v[1]*t,4)]
 def obj(v):return -capacity(C,scales(v),al)/t**2
 res=differential_evolution(obj,[(-1+1e-8,1/(1-t)-1e-8)]*2,seed=20261105+j,popsize=22,maxiter=250,tol=1e-8,atol=1e-8,polish=False)
 return dict(seed=20261105+j,a=aa,t=t,sigma=s,params=res.x.tolist(),slack=-res.fun,n=res.nfev)
if __name__=='__main__':
 jobs=[]
 for a in [.55,.6,.65,.7,.72,.73,.74,.75,.8,.85,.9,.95]:
  b=1-a;q=(2-np.sqrt(4-12*a*b))/3;sig=min((1-4*b*b)/2,2*q-2*q*q)
  for t in [.002,.001]:jobs.append((len(jobs),a,t,sig));st=time.monotonic()
 with ProcessPoolExecutor(max_workers=8) as p:rows=list(p.map(run,jobs))
 out=dict(rows=rows,workers=8,python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__,seconds=time.monotonic()-st,evaluations=sum(x['n'] for x in rows),scope='bounded two-group floating search; negative slack is not an upper certificate');(Path(__file__).resolve().parents[1]/'results/two_group_direction_search.json').write_text(json.dumps(out,indent=2));print(out)
