import numpy as np,scipy,mpmath as m,json,time,platform
from scipy.optimize import differential_evolution
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

def run(job):
 i,n,beta,al=job;p=n-beta;L=al**(1/beta);U=al**(-1/p)
 def F(l,s):return (l+s)/2+s**(1-n)*(1-al*(l**p+s**p))
 def obj(z):
  s=L+(U-L)*z[0];l=s+(U-s)*z[1];return -F(l,s)
 res=differential_evolution(obj,[(0,1),(0,1)],seed=20261105+i,popsize=10,maxiter=140,tol=1e-11,polish=True,workers=1)
 WL=L**(1-n)-L;WU=U-U**(1-n);end=max(WL,WU);val=max(-res.fun,F(L,L),F(U,U));return {'n':n,'beta':beta,'alpha':al,'endpoint_max':end,'numeric_max':val,'excess':val-end,'seed':20261105+i,'evaluations':int(res.nfev)}
if __name__=='__main__':
 jobs=[]
 for n in range(1,7):
  bs=sorted({round(x,8) for x in [n*.1,n*.25,n*.5,n*.75,n*.9,.49,.5,.51,n-.49,n-.5,n-.51,n/2-.001,n/2+.001] if 0<x<n})
  for b in bs:
   for al in [.95,.99,.999]:
    if al>2**(-b/n):jobs.append((len(jobs),n,b,al))
 start=time.monotonic()
 with ProcessPoolExecutor(max_workers=8) as pool:rows=list(pool.map(run,jobs))
 m.mp.dps=70;exact=[]
 for n,b in [(1,m.mpf('.5')),(2,m.mpf(1)),(3,m.mpf('1.5')),(4,m.mpf(1)),(4,m.mpf(3))]:
  for t in ['.01','.001','.0001']:
   t=m.mpf(t);al=1-t;p=n-b;WL=al**(-(n-1)/b)-al**(1/b);WU=al**(-1/p)-al**((n-1)/p);W=max(WL,WU);exact.append({'n':n,'beta':str(b),'t':str(t),'W':str(W),'leading_ratio':str(W/t),'branch':'L' if WL>=WU else 'U'})
 out={'environment':{'python':platform.python_version(),'numpy':np.__version__,'scipy':scipy.__version__,'mpmath':m.__version__},'cases':len(rows),'seed_base':20261105,'workers':8,'optimizer':{'popsize_parameter':10,'max_iterations':140,'tolerance':1e-11},'evaluations':sum(r['evaluations'] for r in rows),'maximum_excess':max(r['excess'] for r in rows),'largest_excess_cases':sorted(rows,key=lambda r:r['excess'],reverse=True)[:10],'failures':[r for r in rows if r['excess']>1e-9],'high_precision_formula_cases':exact,'precision':70,'seconds':time.monotonic()-start,'scope':'finite numerical search and formula evaluation, not a global optimization certificate'}
 (Path(__file__).resolve().parents[1]/'results'/'sharp_width_python.json').write_text(json.dumps(out,indent=2));print({k:v for k,v in out.items() if k not in ['largest_excess_cases','high_precision_formula_cases']})
