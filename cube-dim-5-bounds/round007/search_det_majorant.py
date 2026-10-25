from pathlib import Path
from fractions import Fraction as F
import numpy as np
from scipy.optimize import linprog,brentq
import json,time,hashlib
p=Path(__file__).resolve().parent;start=time.monotonic();qU=F(2457631,10**6);m1=F(648,1463);m2=F(126404864,474619145)
nodes=list(np.exp(-np.linspace(0,80,161)/140));nodes.append(0.)
rows=[]
for step in range(4):
 z=np.array(nodes);mat=np.column_stack([np.ones_like(z),z**30,z**170,z**310]);res=linprog([float(qU),1,float(m1),float(m2)],A_ub=-mat,b_ub=-z**9,bounds=[(0,10)]+[(-100,100)]*3,method='highs',options={'time_limit':30.,'dual_feasibility_tolerance':1e-9,'primal_feasibility_tolerance':1e-9})
 if not res.success:raise RuntimeError(res.message)
 alpha,A,B,C=res.x
 def h(u):return 310*C*u**43+170*B*u**23+30*A*u**3-9
 mesh=np.linspace(0,1,4001);roots=[]
 for u,v in zip(mesh[:-1],mesh[1:]):
  if h(u)*h(v)<0:roots.append(brentq(h,u,v,xtol=1e-15))
 trial=[0.,1.]+[u**(1/7) for u in roots];values=[alpha+A*x**30+B*x**170+C*x**310-x**9 for x in trial];worst=int(np.argmin(values));rows.append({'step':step,'sampled_constraints':len(nodes),'coefficients':res.x.tolist(),'objective':res.fun,'numeric_candidate_critical_points':trial,'numeric_R_values':values})
 nodes.append(trial[worst])
coeff=[str(F(round(float(v)*10**8),10**8)) for v in res.x]
out={'status':'numerical_candidate_not_global_certificate','candidate_rationals_alpha_A_B_C':coeff,'exchange':rows,'elapsed_seconds':time.monotonic()-start,'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
(p/'det_majorant_candidate.json').write_text(json.dumps(out,indent=2)+'\n');print('candidate',coeff,'objective',res.fun,'elapsed',out['elapsed_seconds'])
