"""Exact signed-spline certificate checker; does not import discovery/model.
Uses Python integer objects for tensor arithmetic. Evaluates canonical cells
from full signed lattice, then explicitly maps every support cell by symmetry.
"""
from pathlib import Path
from fractions import Fraction as F
from itertools import product,combinations_with_replacement,permutations
from math import comb,factorial
import json,time,hashlib
import numpy as np
OUT=Path(__file__).resolve().parent
POWER={-2:[0,0,0,1],-1:[1,3,3,-3],0:[4,0,-6,3],1:[1,-3,3,-1]}

def derivative(p,order,t):return sum((F(p[i],6)*F(factorial(i),factorial(i-order))*t**(i-order) for i in range(order,4)),F(0))
def tables():
    B={};DD={}
    for m,p in POWER.items():
        b=[sum((F(p[j],6)*F(comb(i,j),comb(3,j)) for j in range(i+1)),F(0)) for i in range(4)]
        d=[F(p[2],3)+F(i,3)*p[3] for i in range(4)]
        B[m]=np.array([int(v*6) for v in b],dtype=object);DD[m]=np.array([int(v*3) for v in d],dtype=object)
        assert min(b)>=0
        for order in range(3):
            if m<1:assert derivative(p,order,F(1))==derivative(POWER[m+1],order,F(0))
    for order in range(3):assert derivative(POWER[-2],order,F(0))==0 and derivative(POWER[1],order,F(1))==0
    return B,DD

def outer(v):
    a=v[0]
    for b in v[1:]:a=np.multiply.outer(a,b)
    assert a.dtype==object
    return a.reshape(-1)

def check(n,q):
    start=time.monotonic();path=OUT/f'cb8_certificate_n{n}_q{q}.json';r=json.loads(path.read_text());K=r['K'];den=int(r['common_coefficient_denominator']);nums=list(map(int,r['common_coefficient_numerators']))
    points=list(combinations_with_replacement(range(K+1),n));assert r['orbits']==list(map(list,points));assert len(nums)==len(points)==len(r['coefficients'])
    assert all(F(v,den)==F(a) for v,a in zip(nums,r['coefficients']))
    lookup={p:i for i,p in enumerate(points)};weights=[0]*len(points);full={}
    for k in product(range(-K,K+1),repeat=n):
        idx=lookup[tuple(sorted(map(abs,k)))];weights[idx]+=1;full[k]=nums[idx]
    assert weights==r['orbit_multiplicities']
    origin=[comb(n,sum(v==1 for v in p))*2**(n-sum(v==1 for v in p)) if max(p)<=1 else 0 for p in points]
    assert origin==r['origin_integer_row'] and sum(v*a for v,a in zip(origin,nums))==0
    B,DD=tables();pd=6**n;md=3*6**(n-1);basis={}
    for mm in product(range(-2,2),repeat=n):
        bb=[B[m] for m in mm];dd=[DD[m] for m in mm];pp=outer(bb);ll=np.zeros(4**n,dtype=object)
        for j in range(n):ll+=outer([dd[i] if i==j else bb[i] for i in range(n)])
        basis[mm]=(pp,ll)
    cache={}
    for cell in combinations_with_replacement(range(K+2),n):
        pp=np.zeros(4**n,dtype=object);ll=np.zeros(4**n,dtype=object)
        for k in product(*(range(max(-K,v-1),min(K,v+2)+1) for v in cell)):
            a=full[k]
            if not a:continue
            bp,bl=basis[tuple(v-w for v,w in zip(cell,k))];pp+=a*bp;ll+=a*bl
        cache[cell]=(min(pp),max(ll),sum(ll),pp,ll)
    # Reflection/permutation of an invariant spline permutes its Bernstein list.
    # Each full cell is explicitly represented. Coefficient min/max/sum unchanged.
    minpsi=None;maxlap=None;cube_sum=0;cells=0
    for cell in product(range(-K-2,K+2),repeat=n):
        canonical=tuple(sorted(v if v>=0 else -v-1 for v in cell));lo,hi,sm,_,_=cache[canonical]
        assert lo>=0 and q*q*hi<=md*den
        minpsi=lo if minpsi is None else min(minpsi,lo);maxlap=hi if maxlap is None else max(maxlap,hi);cells+=1
        if all(-q<=v<q for v in cell):cube_sum+=sm
    assert cells==(2*K+4)**n
    integral_cells=F(q*q*cube_sum,md*den*4**n*q**n)
    def integ(k):return sum((F(p[j],6*(j+1)) for ell in range(-q,q) if (p:=POWER.get(ell-k)) is not None for j in range(4)),F(0))
    def prime(t):return derivative(POWER[t],1,F(0)) if t in POWER else F(0)
    coeff=[]
    for p,w in zip(points,weights):
        ii=[integ(v) for v in p];jj=[prime(q-v)-prime(-q-v) for v in p]
        v=F(0)
        for j in range(n):
            term=jj[j]
            for k in range(n):
                if k!=j:term*=ii[k]
            v+=term
        coeff.append(-F(w,2**n*q**(n-2))*v)
    assert list(map(str,coeff))==r['objective_coefficients']
    value=1+sum((a*F(v,den) for a,v in zip(coeff,nums)),F(0));assert value==1-integral_cells/2**n==F(r['D_exact'])
    lap0=F(q*q*int(cache[(0,)*n][4][0]),md*den);hdiag=lap0/n
    jet=F(q*q,3**n)*sum((F(c*a,den)*(F(9*sum(v==1 for v in p),n)-3) for p,c,a in zip(points,origin,nums)),F(0))
    assert hdiag==jet==F(r['Hessian_diagonal_at_origin_exact']) and 0<=n*hdiag<=1
    assert F(minpsi,pd*den)==F(r['minimum_psi_Bernstein_exact']) and F(q*q*maxlap,md*den)==F(r['maximum_Delta_Bernstein_exact'])
    alpha=F(r['positive_scale_alpha']);assert value==1+alpha*(F(r['raw_D_exact'])-1)
    return {'n':n,'q':q,'status':'PASS exact complete-cell symmetry-map checker','source_sha256':hashlib.sha256(path.read_bytes()).hexdigest(),'complete_signed_centers':len(full),'canonical_cells_evaluated_from_full_lattice':len(cache),'all_cells_explicitly_mapped':cells,'each_polynomial_coefficients_covered':cells*4**n,'minimum_psi_Bernstein_exact':str(F(minpsi,pd*den)),'maximum_Delta_Bernstein_exact':str(F(q*q*maxlap,md*den)),'origin_value_exact':'0','Hessian_diagonal_exact':str(hdiag),'cube_Delta_integral_exact':str(integral_cells),'D_exact':str(value),'exceeds_target':value>2**n,'seconds':time.monotonic()-start,'shared_dependencies':'original B formula, exact symmetry lemma and Python int/Fraction shared; does not import LP/rebuild or use float signs'}
if __name__=='__main__':
    rows=[]
    for n,q in ((3,1),(3,2),(4,1)):
        row=check(n,q);rows.append(row);print(json.dumps(row),flush=True)
    (OUT/'cb8_check_results.json').write_text(json.dumps({'runs':rows,'scope':'exact feasible points; no exact LP optimality or C_n claim'},indent=2)+'\n')
