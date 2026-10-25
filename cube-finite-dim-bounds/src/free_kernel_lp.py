"""Finite-difference discovery only. No continuous maximal-operator bound.

Hyperoctahedral averaging is exact for this finite convex LP. It does not
certify that finite differences imply a distributional Laplacian inequality.
"""
from pathlib import Path
import itertools, math, json, time
from fractions import Fraction as F
import numpy as np
from scipy.sparse import coo_matrix
from scipy.optimize import linprog
from scipy.io import savemat

ROOT=Path(__file__).resolve().parents[1]

def model(q, closed=True):
    # h=1/q, R=3, all interior coordinates have absolute integer index <3q.
    m=3*q
    points=list(itertools.combinations_with_replacement(range(m),3))
    lookup={p:i for i,p in enumerate(points)}
    weights=[]
    rr=[];cc=[];vv=[];rows=[]
    for i,p in enumerate(points):
        multiplicity=len(set(itertools.permutations(p)))*2**sum(t!=0 for t in p)
        weights.append(multiplicity/q**3)
        if p==(0,0,0):continue
        row=len(rows);rows.append(p);rr.append(row);cc.append(i);vv.append(6)
        for axis in range(3):
            for d in (-1,1):
                x=list(p);x[axis]+=d;x=tuple(sorted(map(abs,x)))
                if x[-1]<m:rr.append(row);cc.append(lookup[x]);vv.append(-1)
    A=coo_matrix((vv,(rr,cc)),shape=(len(rows),len(points))).tocsr()
    low=np.array([.125 if (p[-1]<=q if closed else p[-1]<q) else 0 for p in points])
    return A,np.array(weights),low,points

def solve(q,closed):
    A,c,low,points=model(q,closed);start=time.monotonic()
    r=linprog(c,A_ub=A,b_ub=np.zeros(A.shape[0]),bounds=np.column_stack((low,np.full(len(low),np.inf))),method='highs',options={'time_limit':60.,'dual_feasibility_tolerance':1e-9,'primal_feasibility_tolerance':1e-9})
    out={'q':q,'h':1/q,'R':3,'closed_obstacle':closed,'variables':len(c),'constraints':A.shape[0],'status':int(r.status),'message':r.message,'seconds':time.monotonic()-start}
    if r.success:
        out.update(cost=float(c@r.x),max_laplacian_violation=float(max(0,np.max(A@r.x))),min_obstacle_margin=float(np.min(r.x-low)),dual_value=float(low@r.lower.marginals),dual_stationarity_residual=float(np.max(np.abs(c-A.T@r.ineqlin.marginals-r.lower.marginals))),origin_value=float(r.x[0]))
        if q==4 and closed:
            # External scratch only: numerical matrix witness for local MATLAB.
            savemat('/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/round005_lp.mat',{'A':A,'c':c,'lower':low,'x':r.x,'dual_ineq':r.ineqlin.marginals,'dual_lower':r.lower.marginals})
            # Exact certificate only for this finite LP, not for continuous K.
            y=[F(str(float(v)))*F(999999,1000000) for v in r.ineqlin.marginals]
            assert all(v<=0 for v in y)
            z=[F(float(v)) for v in c]
            for row in range(A.shape[0]):
                for k in range(A.indptr[row],A.indptr[row+1]):
                    z[int(A.indices[k])]-=int(A.data[k])*y[row]
            assert all(v>=0 for v in z)
            exact_dual=sum((F(float(v))*zz for v,zz in zip(low,z)),F(0))
            assert exact_dual>8
            (ROOT/'results/round005_lp_dual.json').write_text(json.dumps({'scope':'finite grid LP only; not a continuous obstruction','q':q,'closed_obstacle':closed,'y':[str(v) for v in y],'dual_bound':str(exact_dual),'strictly_above_8':True},indent=2)+'\n')
            out['exact_finite_dual_bound']=str(exact_dual)
    return out

if __name__=='__main__':
    data={'scope':'numerical finite-difference LP only; no continuous upper/lower bound','model':'min h^3 sum K; K>=cube obstacle; 6K-sum neighbors<=0 off origin; exterior zero at R=3','symmetry':'signs and coordinate permutations; convex group averaging is exact for this finite LP','runs':[]}
    for q in (4,8):
        for closed in (False,True):
            result=solve(q,closed);data['runs'].append(result);print(json.dumps(result),flush=True)
    (ROOT/'results/round005_free_kernel_lp.json').write_text(json.dumps(data,indent=2)+'\n')
