#!/usr/bin/env python3
"""Finite numerical exploration only: 81 fixed-R LPs with full angular bounds."""
from __future__ import annotations
import json
import argparse
from fractions import Fraction as Q
from pathlib import Path
from time import monotonic
import numpy as np
import scipy
from scipy.optimize import linprog, minimize_scalar
from scipy.special import xlogy

ROOT = Path(__file__).resolve().parent
BASELINE = Q(21773155252261,5615797500000)
MAX_FIXED_R_LPS = 81
# 0.00015 mesh, with the angular transition inserted exactly.
RADII = np.unique(np.r_[np.linspace(.1,1.,6001),.75])

def coefficients(r, R):
    r=np.asarray(r,dtype=float)
    s=1-r
    b2=8*(2-r-xlogy(s,s)+(2-r)*np.log(r/2))
    b3=12*(3*(1-r/2)**2+(s*s+1)*np.log(r/2)-s*xlogy(s,s))
    tail=2*(2*R/r**2*(-np.log1p(-r/R)-r/R)+1/R)/(R-1)
    a2=64*s*s/(2-r)**3
    c=3-4*np.log(2)+4*np.log(r)+1/(2*r*r)
    h=1/(R**3*(R-1))
    return b2,b3,tail,a2,c,h

def margins(r,R,e2,e3):
    b2,b3,tail,a2,c,h=coefficients(r,R)
    diag=tail-e2*b2-e3*b3
    angular=diag+(e2*a2+h-3*e3*c)*np.asarray(r)**2
    return diag,angular

def make_lp(R):
    b2,b3,tail,a2,c,h=coefficients(RADII,R)
    angular=RADII>=.75
    A=np.vstack([np.column_stack([b2,b3]),
                 np.column_stack([(b2-a2*RADII**2)[angular],
                                  (b3+3*c*RADII**2)[angular]]),
                 [R*(R-1)*(2-2*np.log(2))]*2,
                 [(R-1)/4]*2])
    b=np.r_[tail,(tail+h*RADII**2)[angular],1.,R]
    labels=([{'family':'diag','r':float(r)} for r in RADII]
            +[{'family':'angular','r':float(r)} for r in RADII[angular]]
            +[{'family':'outer'},{'family':'dominance'}])
    return A,b,labels

def numerical_minima(R,e2,e3):
    out={}
    # Local scalar minimization is diagnostic only, never a certificate.
    for kind,lo,hi in [('diag',.1,1.),('angular',.75,1.)]:
        j=0 if kind=='diag' else 1
        grid=np.linspace(lo,hi,2001)
        vals=margins(grid,R,e2,e3)[j]
        idx=int(np.argmin(vals))
        left=grid[max(0,idx-1)]; right=grid[min(len(grid)-1,idx+1)]
        opt=minimize_scalar(lambda r:float(margins(r,R,e2,e3)[j]),
                            bounds=(left,right),method='bounded',
                            options={'xatol':1e-14,'maxiter':100})
        candidates=[(float(vals[idx]),float(grid[idx])),(float(opt.fun),float(opt.x)),
                    (float(margins(lo,R,e2,e3)[j]),lo),(float(margins(hi,R,e2,e3)[j]),hi)]
        value,r=min(candidates)
        out[kind]={'numerical_minimum':value,'r':r,'not_a_proof':True}
    return out

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--output',type=Path,default=ROOT/'search_result.json')
    path=ap.parse_args().output
    if path.exists():
        raise SystemExit('Search already recorded; refusing to spend another LP budget.')
    start=monotonic(); records=[]
    for i in range(MAX_FIXED_R_LPS):
        rq=Q(17,10)+Q(i,200)
        R=float(rq);A,b,labels=make_lp(R)
        result=linprog([-1/6,-1/10],A_ub=A,b_ub=b,bounds=[(0,None),(0,None)],
                       method='highs',options={'dual_feasibility_tolerance':1e-9,
                                               'primal_feasibility_tolerance':1e-9})
        row={'index':i,'R':str(rq),'status':int(result.status),'success':bool(result.success),
             'message':result.message,'iterations':int(result.nit)}
        if result.success:
            e2,e3=map(float,result.x)
            residual=b-A@result.x
            dual=np.asarray(result.ineqlin.marginals)
            row.update({'e2':e2,'e3':e3,'cost':R*R/(R-1)-e2/6-e3/10,
                        'minimum_grid_slack':float(np.min(residual)),
                        'active_constraints':[{**labels[k],'slack':float(residual[k]),
                                               'dual_marginal':float(dual[k])}
                                              for k in np.flatnonzero(residual<1e-7)],
                        'numeric_dual_objective':float(b@dual),
                        'numeric_primal_objective':float(result.fun),
                        'numeric_dual_gap':float(result.fun-b@dual),
                        'numeric_dual_feasibility_residual':float(np.max(A.T@dual-[-1/6,-1/10]))})
        records.append(row)
    feasible=[x for x in records if x['success']]
    best=min(feasible,key=lambda x:x['cost'])
    # A fixed inward contraction gives room for exact interval dependency widths.
    factor=Q(999,1000)
    e2q=Q(int(np.floor(best['e2']*float(factor)*10**6)),10**6)
    e3q=Q(int(np.floor(best['e3']*float(factor)*10**6)),10**6)
    rq=Q(best['R']);cost=rq*rq/(rq-1)-e2q/6-e3q/10
    candidate={'R':str(rq),'e2':str(e2q),'e3':str(e3q),
               'selection':'best of 81 fixed-R numerical LPs; multiply e2,e3 by 999/1000 then floor at denominator 1000000',
               'cost':str(cost),'cost_display':float(cost),
               'improvement_over_round008':str(BASELINE-cost),
               'improvement_over_round008_display':float(BASELINE-cost),
               'numerical_minima':numerical_minima(float(rq),float(e2q),float(e3q))}
    data={'status':'numerical_exploration_only','solver':'scipy.optimize.linprog / HiGHS',
          'scipy_version':scipy.__version__,'numpy_version':np.__version__,
          'LP_count':len(records),'LP_budget':MAX_FIXED_R_LPS,
          'R_range':['17/10','21/10'],'R_step':'1/200','radial_grid_count':len(RADII),
          'radial_grid':'linspace(0.1,1,6001) union {0.75}',
          'constraints':'e2,e3>=0; Gdiag>=0 for all grid r; Gdiag+k*r^2>=0 for grid r>=.75; outer; kernel dominance',
          'baseline_round008':str(BASELINE),'best_numerical':best,
          'rational_candidate':candidate,'records':records,'elapsed_seconds':monotonic()-start}
    path.write_text(json.dumps(data,indent=2)+'\n')
    print(json.dumps({k:v for k,v in data.items() if k!='records'},indent=2))
if __name__=='__main__':main()
