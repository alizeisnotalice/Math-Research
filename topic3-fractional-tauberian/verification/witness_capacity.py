from fractions import Fraction as F
import itertools,json,time,platform,scipy
from pathlib import Path
import numpy as np
from scipy.optimize import linprog

def atoms(boxes):
 xs=sorted({x for b in boxes for x in b[:2]});ys=sorted({y for b in boxes for y in b[2:]});caps={}
 for a,b in zip(xs,xs[1:]):
  for c,d in zip(ys,ys[1:]):
   x=(a+b)/2;y=(c+d)/2;mask=sum(1<<i for i,(l,r,lo,hi) in enumerate(boxes) if l<x<r and lo<y<hi)
   if mask:caps[mask]=caps.get(mask,F(0))+(b-a)*(d-c)
 return caps

def solve(caps,thresholds,denominator_limit=10**6,solver_options=None):
 n=len(thresholds);ss=list(caps);c=np.array([float(caps[s]) for s in ss]);A=np.array([[int(s>>i&1) for s in ss] for i in range(n)],float);b=np.array(list(map(float,thresholds)))
 obj=np.r_[np.zeros(len(ss)),-1];pr=linprog(obj,A_ub=np.vstack([np.r_[np.ones(len(ss)),0],np.c_[-A,np.ones(n)]]),b_ub=np.r_[1,-b],bounds=[(0,z) for z in c]+[(None,None)],method='highs',options=solver_options);assert pr.success
 # dual variables lambda(n), z, h(atom count)
 ob=np.r_[-b,1,c];DU=np.c_[A.T,-np.ones(len(ss)),-np.eye(len(ss))];eq=np.r_[np.ones(n),0,np.zeros(len(ss))][None,:]
 du=linprog(ob,A_ub=DU,b_ub=np.zeros(len(ss)),A_eq=eq,b_eq=[1],bounds=(0,None),method='highs',options=solver_options);assert du.success
 lam=[F(float(z)).limit_denominator(denominator_limit) for z in du.x[:n]];tot=sum(lam);lam=[x/tot for x in lam];z=max(F(0),F(float(du.x[n])).limit_denominator(denominator_limit))
 dual=z+sum(caps[s]*max(F(0),sum(lam[i] for i in range(n) if s>>i&1)-z) for s in ss)-sum(x*y for x,y in zip(lam,thresholds))
 m=[min(caps[s],max(F(0),F(float(x)).limit_denominator(denominator_limit))) for s,x in zip(ss,pr.x[:-1])];scale=max(F(1),sum(m));m=[x/scale for x in m]
 slack=min(sum(v for s,v in zip(ss,m) if s>>i&1)-thresholds[i] for i in range(n))
 return {'eta_float':float(pr.x[-1]),'dual_float':float(du.fun),'duality_gap':abs(float(pr.x[-1])-float(du.fun)),'dual_certificate':str(dual),'primal_strict_slack':str(slack),'status':'exact-infeasible' if dual<=0 else ('exact-feasible' if slack>0 else 'unresolved'),'primal_masses':{str(s):str(v) for s,v in zip(ss,m) if v},'lambda':[str(x) for x in lam],'z':str(z)}
if __name__=='__main__':
 start=time.monotonic();rows=[]
 for alpha in [F(9,10),F(19,20),F(99,100)]:
  t=1-alpha
  for ff in [F(1,100),F(1,10),F(1,2),F(9,10),F(99,100)]:
   ell=alpha+ff*(1/alpha-alpha)
   for rr in [F(4,5),F(19,20),F(1),F(41,40),F(21,20),F(11,10),F(6,5)]:
    a=t*rr;centers=[(a,0),(-a,0),(0,a),(0,-a),(a/2,a/2),(-a/2,a/2),(a/2,-a/2),(-a/2,-a/2)]
    boxes=[(x-ell/2,x+ell/2,y-ell/2,y+ell/2) for x,y in centers];cp=atoms(boxes);res=solve(cp,[alpha*ell]*8)
    if res['status']=='exact-feasible':
     mm={int(s):F(v) for s,v in res['primal_masses'].items()}
     radii=[(sum(v for s,v in mm.items() if s>>i&1)-alpha*ell)/(2*ell) for i in range(8)]
     lower_boxes=[(x-r,x+r,y-r,y+r) for (x,y),r in zip(centers,radii)]
     res['certified_halo_lower_bound']=str(sum(atoms(lower_boxes).values()))
    rows.append({'alpha':str(alpha),'ell':str(ell),'axis_radius':str(a),'atoms':len(cp),**res})
 example=solve({3:F(3,10),5:F(3,10),6:F(3,10),1:F(1,10),2:F(1,10),4:F(1,10)},[F(16,25)]*3)
 out={'environment':{'python':platform.python_version(),'numpy':np.__version__,'scipy':scipy.__version__},'solver':'HiGHS primal/dual LP','rational_denominator_limit':10**6,'cases':len(rows),'rows':rows,'example':example,'elapsed_seconds':time.monotonic()-start,'counts':{s:sum(r['status']==s for r in rows) for s in ['exact-feasible','exact-infeasible','unresolved']},'max_duality_gap':max(r['duality_gap'] for r in rows),'scope':'finite fixed witness families; does not certify continuous-scale halo upper bound'}
 (Path(__file__).resolve().parents[1]/'results'/'witness_capacity.json').write_text(json.dumps(out,indent=2));print({k:v for k,v in out.items() if k!='rows'})
