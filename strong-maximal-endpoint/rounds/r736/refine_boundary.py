import os
for k in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS']:os.environ[k]='1'
import json,time,sys,importlib.util
from pathlib import Path
from fractions import Fraction as F
import numpy as np
HERE=Path(__file__).resolve().parent;BASE=Path(sys.argv[1]) if len(sys.argv)>1 else HERE.parent/'r732';spec=importlib.util.spec_from_file_location('capmath',BASE/'capacity.py');cap=importlib.util.module_from_spec(spec);spec.loader.exec_module(cap)
d=json.loads((HERE/'lf_capacities.json').read_text());changed=[]
for r in d['vertex_certificates']:
 if F(r['certified_capacity_gap_upper'])<F(1,1000000):continue
 started=time.monotonic();mass=list(map(F,r['atom_masses']));vol=list(map(F,r['row_volumes']));masks=r['atom_masks'];pi=list(map(F,r['pi']));active=r['active_rows'];pir=[pi[i] for i in active];q=[m/sum(mass) for m in mass];ker=[[m/vol[i] if s>>i&1 else F(0) for m,s in zip(mass,masks)] for i in active]
 ar=np.array([[float(x) for x in v] for v in ker]);pr=np.array([float(x) for x in pir]);qr=np.array([float(x) for x in q]);nu=np.array([float(F(x)) for x in r['nu']]);trace=[]
 for iteration in range(100001):
  p=nu@ar;g=ar@np.log(p/qr)-np.log(nu/pr);value=float(nu@g);gap=float(max(g)-value)
  if iteration in [0,10,100,1000,10000,100000] or gap<1e-10:trace.append(dict(iteration=iteration,value=value,gradient_gap=gap,nu=nu.tolist()))
  if gap<1e-10:break
  logits=np.log(pr)+ar@np.log(p/qr);logits-=max(logits);nu=np.maximum(np.exp(logits),1e-280);nu/=sum(nu)
 # Preserve the converged ratios in small components instead of flooring all labels equally.
 rn=[F(str(x)) for x in nu];rn=[x/sum(rn) for x in rn];pn=[sum(rn[i]*ker[i][j] for i in range(len(rn))) for j in range(len(q))]
 lo,hi=cap.linear_interval([(pn[j],pn[j]/q[j]) for j in range(len(q))]+[(-rn[i],rn[i]/pir[i]) for i in range(len(rn))]);gg=[cap.linear_interval([(ker[i][j],pn[j]/q[j]) for j in range(len(q))]+[(F(-1),rn[i]/pir[i])]) for i in range(len(rn))];upper=max(t[1] for t in gg)
 r.update(nu=[str(x) for x in rn],capacity_lower=cap.decimal_bound(lo),capacity_upper=cap.decimal_bound(upper,True),certified_capacity_gap_upper=cap.decimal_bound(upper-lo,True),solver_status='positive multiplicative capacity iteration with gradient-gap stopping; preserves small-component ratios',solver_success=gap<1e-10)
 changed.append(dict(source_case=r['source_case'],vertex=r['vertex'],iterations=iteration,trace=trace,final_interval_width=r['certified_capacity_gap_upper'],elapsed_seconds=time.monotonic()-started));print(json.dumps(changed[-1]),flush=True)
for s in d['LF_capacity_summaries']:
 rows=[r for r in d['vertex_certificates'] if r['source_case']==s['case']];lo=max(F(r['capacity_lower']) for r in rows);hi=max(F(r['capacity_upper']) for r in rows);s.update(LF_capacity_lower=cap.decimal_bound(lo),LF_capacity_upper=cap.decimal_bound(hi,True),LF_capacity_interval_width_upper=cap.decimal_bound(hi-lo,True))
d['boundary_refinement']=d.get('boundary_refinement',[])+changed;(HERE/'lf_capacities.json').write_text(json.dumps(d,indent=2)+'\n')
