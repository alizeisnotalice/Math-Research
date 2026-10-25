"""Exact integer Bernstein check; no solver and no imports of agent code."""
from pathlib import Path
from itertools import product
from fractions import Fraction as F
from math import comb
import json, hashlib
import numpy as np
import sympy as s
u=s.symbols('u'); D=18

def Bpoly(cell,shift=0):
    z=u+cell-shift
    j=cell-shift
    if j==-2:return (z+2)**3/6
    if j==-1:return (4-6*z*z-3*z**3)/6
    if j==0:return (4-6*z*z+3*z**3)/6
    if j==1:return (2-z)**3/6
    return s.Integer(0)

def bern(p):
    a=s.Poly(s.expand(p),u)
    vals=[sum((F(a.nth(j))*F(comb(i,j),comb(3,j)) for j in range(i+1)),F(0)) for i in range(4)]
    assert all((v*D).denominator==1 for v in vals)
    return np.array([int(v*D) for v in vals],dtype=np.int64)

def tensor(rows):
    out=np.array(1,dtype=np.int64)
    for row in rows:out=np.multiply.outer(out,row)
    return out

def run():
    rows=[]
    for k in range(-3,3):
        b=Bpoly(k); phi=s.Rational(2,3)*(Bpoly(k,1)+Bpoly(k,-1))-b/3
        rows.append([bern(b),bern(s.diff(b,u,2)),bern(phi),bern(s.diff(phi,u,2))])
    assert all(min(r[0])>=0 and min(r[2])>=0 for r in rows)
    assert s.diff(s.Rational(2,3)*(Bpoly(0,1)+Bpoly(0,-1))-Bpoly(0)/3,u,2).subs(u,0)==2
    result=[]
    for n in (3,4):
        minpsi=0;maxlap=0;minlap=0;integral=0;cells=0
        # A crude integer bound for intermediate sums is far below signed int64.
        maxentry=max(abs(int(v)) for row in rows for arr in row for v in arr)
        assert n*n*maxentry**n<2**63
        for c in product(range(6),repeat=n):
            rr=[rows[j] for j in c];p=np.zeros((4,)*n,dtype=np.int64);lap=p.copy()
            for j in range(n):
                p+=tensor([rr[k][2 if k==j else 0] for k in range(n)])
                for ell in range(n):
                    lap+=tensor([rr[k][3 if k==j==ell else 2 if k==j else 1 if k==ell else 0] for k in range(n)])
            minpsi=min(minpsi,int(p.min()));maxlap=max(maxlap,int(lap.max()));minlap=min(minlap,int(lap.min()));cells+=1
            if all(j in (2,3) for j in c):integral+=int(lap.sum())
        assert minpsi>=0
        C=F(maxlap,D**n);scale=1/max(F(1),C)
        integ=F(integral,D**n*4**n)
        hess=scale*2*F(2,3)**(n-1)
        assert hess>0
        result.append({'n':n,'full_cells':cells,'coefficients_per_polynomial':cells*4**n,'psi_Bernstein_min':str(F(minpsi,D**n)),'raw_Delta_Bernstein_max':str(C),'raw_Delta_Bernstein_min':str(F(minlap,D**n)),'scale':str(scale),'scaled_origin_Hessian_diagonal':str(hess),'cube_Delta_integral_raw':str(integ),'D_exact':str(1-scale*integ/F(2**n)),'origin_value':'0','scope':'continuous feasible signed-coefficient example, not a strong dual or C_n bound'})
    return {'status':'PASS','basis_denominator':D,'results':result,'libraries':{'numpy':np.__version__,'sympy':s.__version__},'checker_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
if __name__=='__main__':
    out=run();p=Path(__file__).resolve().parents[1]/'results/round007_signed_interface.json';p.write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
