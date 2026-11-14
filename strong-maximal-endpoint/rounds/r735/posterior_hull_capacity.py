#!/usr/bin/env python3
import os
for key in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS']:os.environ[key]='1'
import json,time,sys,importlib.util,hashlib
from pathlib import Path
from fractions import Fraction as F
from concurrent.futures import ProcessPoolExecutor
import numpy as np
from scipy.optimize import minimize
HERE=Path(__file__).resolve().parent
BASE=Path(sys.argv[1]) if len(sys.argv)>1 else HERE.parent/'r732'
spec=importlib.util.spec_from_file_location('capacity_math',BASE/'capacity.py');cap=importlib.util.module_from_spec(spec);spec.loader.exec_module(cap)

def run(r):
 mass=list(map(F,r['atom_masses']));vol=list(map(F,r['row_volumes']));masks=r['atom_masks'];pi=list(map(F,r['pi']));active=r['active_rows'];pir=[pi[i] for i in active];om=sum(mass);q=[m/om for m in mass];k=len(active);a=len(mass)
 amat=[[mass[j]/vol[i] if masks[j]>>i&1 else F(0) for j in range(a)] for i in active]
 pm=[sum(pir[i]*amat[i][j] for i in range(k)) for j in range(a)]
 alpha=[[pir[i]*amat[i][j]/pm[j] for i in range(k)] for j in range(a)]
 ar=np.array([[float(x) for x in row] for row in amat]);al=np.array([[float(x) for x in row] for row in alpha]);pr=np.array([float(x) for x in pir]);qr=np.array([float(x) for x in q]);z0=np.array([float(x) for x in pm])
 def fg(z):
  nu=z@al;pp=nu@ar;g=ar@np.log(pp/qr)-np.log(nu/pr);value=float(nu@g);return value,al@g
 starts=[z0,np.full(a,1/a)]+[.98*np.eye(a)[j]+.02*z0 for j in range(a)]
 out=[]
 for x0 in starts:
  rr=minimize(lambda z:-fg(z)[0],x0,jac=lambda z:-fg(z)[1],bounds=[(1e-14,1)]*a,constraints=[dict(type='eq',fun=lambda z:sum(z)-1,jac=lambda z:np.ones(a))],method='SLSQP',options=dict(ftol=1e-14,maxiter=3000))
  if np.isfinite(rr.fun):out.append(rr)
 rr=min(out,key=lambda z:z.fun);assert -rr.fun>=fg(z0)[0]-1e-8
 counts=[max(1,int(round(x*10**12))) for x in rr.x];z=[F(c,sum(counts)) for c in counts];nu=[sum(z[j]*alpha[j][i] for j in range(a)) for i in range(k)];pn=[sum(nu[i]*amat[i][j] for i in range(k)) for j in range(a)]
 flo,fhi=cap.linear_interval([(pn[j],pn[j]/q[j]) for j in range(a)]+[(-nu[i],nu[i]/pir[i]) for i in range(k)])
 gis=[cap.linear_interval([(amat[i][j],pn[j]/q[j]) for j in range(a)]+[(F(-1),nu[i]/pir[i])]) for i in range(k)]
 his=[sum(alpha[j][i]*gis[i][1] for i in range(k)) for j in range(a)];upper=max(his)
 full_lo=F(r['capacity_lower']);full_hi=F(r['capacity_upper']);gap=full_lo-upper
 return dict(case=r['case'],target=r['target'],space_weights=[str(x) for x in z],induced_nu=[str(x) for x in nu],hull_capacity_lower=cap.decimal_bound(flo),hull_capacity_upper=cap.decimal_bound(upper,True),hull_interval_width_upper=cap.decimal_bound(upper-flo,True),full_capacity_lower=r['capacity_lower'],full_capacity_upper=r['capacity_upper'],full_minus_hull_lower=cap.decimal_bound(gap),certified_strict_loss=gap>0,solver_success=bool(rr.success),solver_status=str(rr.message))
if __name__=='__main__':
 start=time.monotonic();raw=(BASE/'capacity.json').read_bytes();rows=json.loads(raw)['results']
 with ProcessPoolExecutor(max_workers=12) as pool:res=list(pool.map(run,rows))
 data=dict(status='complete',scope='Fixed LFpi finite catalogue capacities restricted to arbitrary convex mixtures of Bayes atom posteriors; no central asymptotic claim',source_sha256=hashlib.sha256(raw).hexdigest(),workers=12,log_terms=cap.TERMS,certificate='F lower; max_j alpha(j) dot gradient upper, exact rational log series',wall_seconds=time.monotonic()-start,results=res)
 (HERE/'posterior_hull_capacity.json').write_text(json.dumps(data,indent=2)+'\n');print(json.dumps(dict(cases=len(res),wall_seconds=data['wall_seconds'],strict_losses=[(r['case'],r['target'],r['full_minus_hull_lower']) for r in res if r['certified_strict_loss']],max_interval_width=max(float(r['hull_interval_width_upper']) for r in res))))
