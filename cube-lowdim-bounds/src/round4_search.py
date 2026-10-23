"""Round 4: finite discovery over an additional scaled quadratic bump.
Grid feasibility is discovery ONLY; candidates need continuous certification.
Budget: 72 coarse pairs + bounded differential evolution, fixed seed.
"""
import json,time
from pathlib import Path
import numpy as np
from scipy.optimize import linprog,differential_evolution,minimize_scalar
ROOT=Path(__file__).resolve().parents[1]
b=2-2*np.log(2)
def J(s):return .5-2/s-4/s**2*np.log1p(-s/2)
def T(r):return 2*(2-r)*(1-np.log(2)+np.log(r))-2*(1-r)*np.log1p(-r)
def F(r):
    out=np.full_like(r,b);inside=r<1;out[inside]=T(r[inside]);return out
def solve(R,a,n=500):
    r=np.unique(np.r_[np.linspace(.3,.9999,n),a*np.linspace(.5,.99999,n//2),a,1.])
    A=R*(R-1)*np.column_stack((F(r),F(r/a)/a))
    A=np.vstack((A,R*(R-1)*b*np.array([1.,1/a])))
    rhs=np.r_[J(2*r/R),1.]
    # The dominance constraint is uniform and generous.
    A=np.vstack((A,np.ones(2)));rhs=np.r_[rhs,4*R/(R-1)]
    res=linprog(-np.array([1.,a*a]),A_ub=A,b_ub=rhs,bounds=[(0,None)]*2,method='highs')
    if not res.success:return 1e3,None
    return R*R/(R-1)+res.fun/6,res.x
def run():
    start=time.time();best=[10,None];calls=0
    def objective(p):
        nonlocal calls,best
        calls+=1;val,e=solve(*p)
        if val<best[0]:
            best=[val,{'R':float(p[0]),'a':float(p[1]),'e':e.tolist()}]
            print(calls,best,flush=True)
        return val
    for R in np.linspace(1.8,2,8):
        for a in np.linspace(.65,.95,9):objective((R,a))
    res=differential_evolution(objective,[(1.7,2.1),(.5,.98)],seed=20261024,popsize=8,maxiter=35,tol=1e-10,polish=False)
    val,e=solve(*res.x,n=5000)
    out={'best_grid_500':best,'refined_grid_5000':{'R':float(res.x[0]),'a':float(res.x[1]),'e':e.tolist(),'cost':float(val)},'calls':calls,'elapsed':time.time()-start,'seed':20261024,'status':'numerical_candidate_only'}
    (ROOT/'results/round4_search.json').write_text(json.dumps(out,indent=2));print(json.dumps(out,indent=2))
if __name__=='__main__':run()
