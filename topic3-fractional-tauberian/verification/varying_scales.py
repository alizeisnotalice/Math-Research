from fractions import Fraction as F
from witness_capacity import atoms,solve
from pathlib import Path
import itertools,json,time,platform
import numpy as np,scipy

def run():
 start=time.monotonic();rows=[]
 for al in [F(9,10),F(19,20),F(99,100)]:
  t=1-al
  for rr in [F(19,20),F(1),F(41,40),F(21,20),F(11,10)]:
   a=t*rr;C=[(a,0),(-a,0),(0,a),(0,-a),(a/2,a/2),(-a/2,a/2),(a/2,-a/2),(-a/2,-a/2)]
   for fs in itertools.product([F(1,20),F(1,2),F(19,20)],repeat=3):
    ex,ey,ed=[al+f*(1/al-al) for f in fs];ell=[ex,ex,ey,ey,ed,ed,ed,ed];R=[(x-l/2,x+l/2,y-l/2,y+l/2) for (x,y),l in zip(C,ell)]
    res=solve(atoms(R),[al*l for l in ell])
    if res['status']=='unresolved':
     res=solve(atoms(R),[al*l for l in ell],denominator_limit=10**12,solver_options={'primal_feasibility_tolerance':1e-10,'dual_feasibility_tolerance':1e-10});res['refined']=True
    rows.append({'alpha':str(al),'radius':str(a),'scales':[str(ex),str(ey),str(ed)],**res})
 summary={'cases':len(rows),'counts':{s:sum(r['status']==s for r in rows) for s in ['exact-feasible','exact-infeasible','unresolved']},'seconds':time.monotonic()-start,'environment':{'python':platform.python_version(),'numpy':np.__version__,'scipy':scipy.__version__},'parameters':{'alpha':['9/10','19/20','99/100'],'radius_over_t':['19/20','1','41/40','21/20','11/10'],'scale_fraction':['1/20','1/2','19/20']},'scope':'fixed eight-center families with independently chosen x-axis/y-axis/diagonal scales; not a continuum upper bound'}
 best=[]
 for al in ['9/10','19/20','99/100']:
  good=[r for r in rows if r['alpha']==al and r['status']=='exact-feasible']
  if good:best.append(max(good,key=lambda r:(F(r['radius']),F(r['primal_strict_slack']))))
 summary['selected_certificates']=best
 summary['unresolved_cases']=[r for r in rows if r['status']=='unresolved']
 summary['refined_cases']=[r for r in rows if r.get('refined')]
 summary['refinement']={'rational_denominator_limit':10**12,'primal_dual_feasibility_tolerance':1e-10}
 (Path(__file__).resolve().parents[1]/'results'/'varying_scales.json').write_text(json.dumps(summary,indent=2));print(summary)
if __name__=='__main__':run()
