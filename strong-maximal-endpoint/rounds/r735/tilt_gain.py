#!/usr/bin/env python3
import os
for k in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS']:os.environ[k]='1'
import json,time,sys,importlib.util,hashlib
from pathlib import Path
from fractions import Fraction as F
from concurrent.futures import ProcessPoolExecutor
HERE=Path(__file__).resolve().parent
BASE=Path(sys.argv[1]) if len(sys.argv)>1 else HERE.parent/'r732'
spec=importlib.util.spec_from_file_location('capacity_math',BASE/'capacity.py');cap=importlib.util.module_from_spec(spec);spec.loader.exec_module(cap)
def run(r):
 mass=list(map(F,r['atom_masses']));vol=list(map(F,r['row_volumes']));masks=r['atom_masks'];pi=list(map(F,r['pi']));active=r['active_rows'];pir=[pi[i] for i in active];om=sum(mass);q=[m/om for m in mass];k=len(active);a=len(mass)
 ker=[[mass[j]/vol[i] if masks[j]>>i&1 else F(0) for j in range(a)] for i in active]
 pm=[sum(pir[i]*ker[i][j] for i in range(k)) for j in range(a)];wp=[pm[j]/q[j] for j in range(a)];alpha=[[pir[i]*ker[i][j]/pm[j] for i in range(k)] for j in range(a)]
 gs=[cap.linear_interval([(ker[i][j],wp[j]) for j in range(a)]) for i in range(k)]
 rounded=[]
 for j in range(a):
  midpoint=sum(alpha[j][i]*(gs[i][0]+gs[i][1])/2 for i in range(k));scaled=midpoint*10**12+F(1,2);rounded.append(F(scaled.numerator//scaled.denominator,10**12))
 center=sum(pm[j]*rounded[j] for j in range(a));phi=[x-center for x in rounded];mx=max(map(abs,phi));t=min(F(1,2),1/(2*mx)) if mx else F(1,2)
 z=[pm[j]*(1+t*phi[j]) for j in range(a)];assert min(z)>0 and sum(z)==1 and t*mx<=F(1,2)
 nu=[sum(z[j]*alpha[j][i] for j in range(a)) for i in range(k)];pn=[sum(nu[i]*ker[i][j] for i in range(k)) for j in range(a)]
 flo,fhi=cap.linear_interval([(pn[j],pn[j]/q[j]) for j in range(a)]+[(-nu[i],nu[i]/pir[i]) for i in range(k)])
 gi=[cap.linear_interval([(ker[i][j],pn[j]/q[j]) for j in range(a)]+[(F(-1),nu[i]/pir[i])]) for i in range(k)]
 upper=max(sum(alpha[j][i]*gi[i][1] for i in range(k)) for j in range(a))
 linear_coeff=[sum((nu[i]-pir[i])*ker[i][j] for i in range(k)) for j in range(a)];dl,du=cap.linear_interval(list(zip(linear_coeff,wp)));penalty=t*t*sum(pm[j]*phi[j]**2 for j in range(a));analytic=dl-penalty
 gain=flo-F(r['KL_pi_upper']);gap=F(r['capacity_lower'])-upper
 return dict(case=r['case'],target=r['target'],space_weights=[str(x) for x in z],induced_nu=[str(x) for x in nu],rounded_B=[str(x) for x in rounded],phi=[str(x) for x in phi],step_t=str(t),hull_capacity_lower=cap.decimal_bound(flo),hull_capacity_upper=cap.decimal_bound(upper,True),hull_interval_width_upper=cap.decimal_bound(upper-flo,True),full_capacity_lower=r['capacity_lower'],full_capacity_upper=r['capacity_upper'],full_minus_hull_lower=cap.decimal_bound(gap),certified_strict_loss=gap>0,direct_gain_lower=cap.decimal_bound(gain),analytic_direction_gain_lower=cap.decimal_bound(analytic),certified_strict_gain=gain>0)
if __name__=='__main__':
 start=time.monotonic();raw=(BASE/'capacity.json').read_bytes();rows=json.loads(raw)['results']
 with ProcessPoolExecutor(max_workers=12) as pool:res=list(pool.map(run,rows))
 data=dict(status='complete',scope='Explicit rational Bayes-space tilts at fixed LFpi; certified gain only, not optimality or central asymptotics',workers=12,source_sha256=hashlib.sha256(raw).hexdigest(),rounding_scale=10**12,log_terms=48,wall_seconds=time.monotonic()-start,results=res)
 (HERE/'tilt_gain.json').write_text(json.dumps(data,indent=2)+'\n');print(json.dumps(dict(status='complete',cases=len(res),strict_gains=sum(x['certified_strict_gain'] for x in res),max_gain_lower=max(float(x['direct_gain_lower']) for x in res),analytic_positive=sum(F(x['analytic_direction_gain_lower'])>0 for x in res),wall_seconds=data['wall_seconds'])))
