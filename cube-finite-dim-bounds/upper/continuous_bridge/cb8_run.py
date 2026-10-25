"""Three preregistered signed continuous-spline LPs with exact facet rebuilds.
No imported earlier discovery/checker; shared cardinal spline definitions.
"""
from pathlib import Path
from fractions import Fraction as F
from itertools import combinations_with_replacement, permutations, product
from math import comb, gcd, lcm
import json,time,hashlib
import numpy as np
import sympy as sp
from scipy.sparse import csr_matrix,vstack
from scipy.optimize import linprog
OUT=Path(__file__).resolve().parent
BN=np.array([[0,0,0,1],[1,2,4,4],[4,4,2,1],[1,0,0,0]],dtype=np.int64)
DN=np.array([[0,1,2,3],[3,0,-3,-6],[-6,-3,0,3],[3,2,1,0]],dtype=np.int64)

def tensor(vecs):
    a=vecs[0]
    for b in vecs[1:]:a=np.multiply.outer(a,b)
    return a

def model(n,q,K):
    points=list(combinations_with_replacement(range(K+1),n));lookup={p:i for i,p in enumerate(points)}
    cells=list(combinations_with_replacement(range(K+2),n));pc=[];mc=[]
    for cell in cells:
        P=np.zeros((4**n,len(points)),dtype=np.int64);M=P.copy()
        choices=[range(max(-K,l-1),min(K,l+2)+1) for l in cell]
        for k in product(*choices):
            col=lookup[tuple(sorted(map(abs,k)))];mm=[l-v+2 for l,v in zip(cell,k)]
            bs=[BN[m] for m in mm];ds=[DN[m] for m in mm]
            P[:,col]+=tensor(bs).reshape(4**n)
            val=np.zeros((4,)*n,dtype=np.int64)
            for j in range(n):val+=tensor([ds[i] if i==j else bs[i] for i in range(n)])
            M[:,col]+=val.reshape(4**n)
        pc.append(P);mc.append(M)
    P=np.vstack(pc);M=np.vstack(mc);pd=6**n;md=3*6**(n-1)
    boundp=4**(2*n);boundm=n*6*4**(2*n-1)
    assert np.max(np.abs(P))<=boundp and np.max(np.abs(M))<=boundm
    # Exact psi(0) row times 3^n, including every central orbit multiplicity.
    origin=[comb(n,sum(v==1 for v in p))*2**(n-sum(v==1 for v in p)) if max(p)<=1 else 0 for p in points]
    elim=lookup[(1,)*n];assert origin[elim]==1
    free=[i for i in range(len(points)) if i!=elim]
    Pred=P[:,free]-P[:,elim,None]*np.array(origin,dtype=np.int64)[free]
    Mred=M[:,free]-M[:,elim,None]*np.array(origin,dtype=np.int64)[free]
    assert max(boundp,boundm)*(1+max(origin))*4**n <2**63
    def integ(k):return sum((F(int(sum(BN[l-k+2])),24) for l in range(-q,q) if -2<=l-k<=1),F(0))
    def prime(t):return F(1,2) if t==-1 else -F(1,2) if t==1 else F(0)
    c=[];weights=[]
    for p in points:
        weight=len(set(permutations(p)))*2**sum(v!=0 for v in p);weights.append(weight)
        ii=[integ(v) for v in p];jj=[prime(q-v)-prime(-q-v) for v in p]
        integral=sum((jj[j]*__import__('functools').reduce(lambda a,b:a*b,[ii[l] for l in range(n) if l!=j],F(1)) for j in range(n)),F(0))
        c.append(-F(weight,2**n*q**(n-2))*integral)
    alt=[F(0) for _ in points]
    for cell,block in zip(cells,mc):
        if cell[-1]>=q:continue
        mult=2**n*len(set(permutations(cell)))
        for j in range(len(points)):alt[j]-=F(mult*int(block[:,j].sum()),2**n*q**(n-2)*md*4**n)
    assert c==alt and sum(weights)==(2*K+1)**n
    cred=[c[i]-origin[i]*c[elim] for i in free]
    return points,cells,P,M,Pred,Mred,origin,elim,free,c,cred,weights,pd,md

def dot(row,nums):return sum(int(v)*nums[j] for j,v in enumerate(row) if v)

def modular_basis(rows,nvar):
    prime=1000003;basis={};actual=[]
    for row in rows:
        v=np.array(row,dtype=np.int64)%prime
        for piv,bb in sorted(basis.items()):
            if v[piv]:v=(v-int(v[piv])*bb)%prime
        nz=np.flatnonzero(v)
        if len(nz):
            piv=int(nz[0]);v=(v*pow(int(v[piv]),-1,prime))%prime;basis[piv]=v;actual.append(list(map(int,row)))
            if len(actual)==nvar:break
    return actual

def exact_face_rebuild(active,z0,roundden):
    nvar=len(z0);basis=modular_basis(active,nvar);closure=0
    while True:
        if basis:rref,piv=sp.Matrix(basis).rref()
        else:rref=sp.zeros(0,nvar);piv=()
        free=[j for j in range(nvar) if j not in piv]
        b=[F(0) for _ in range(nvar)]
        for j in free:b[j]=F(int(round(float(z0[j])*roundden)),roundden)
        for rr,j in enumerate(piv):b[j]=-sum((F(rref[rr,k])*b[k] for k in free),F(0))
        den=lcm(*(a.denominator for a in b));nums=[int(a*den) for a in b]
        missed=next((row for row in active if dot(row,nums)!=0),None)
        if missed is None:return nums,den,{'active_rank':len(piv),'active_unique_rows':len(active),'free_parameters':len(free),'rank_closure_additions':closure}
        basis.append(missed);closure+=1
        assert closure<=nvar

def solve(n,q,K):
    start=time.monotonic();data=model(n,q,K)
    points,cells,P,M,Pr,Mr,origin,elim,free,c,cr,weights,pd,md=data
    # Deduplication is exact; full rows remain for the final checks.
    Pu=np.unique(Pr,axis=0);Pu=Pu[np.any(Pu!=0,axis=1)]
    Mu=np.unique(Mr,axis=0);Mu=Mu[np.any(Mu!=0,axis=1)]
    # New matrix regression only: the already-audited q1 signed example.
    regression=None
    if q==1:
        a=[F(0) for _ in points];a[points.index((0,)*n)]=F(-3,8) if n==3 else F(-9,16);a[points.index((0,)*(n-1)+(1,))]=F(1,4) if n==3 else F(9,32)
        rd=lcm(*(v.denominator for v in a));an=[int(v*rd) for v in a]
        assert sum(v*w for v,w in zip(origin,an))==0
        assert min(dot(row,an) for row in P)>=0 and max(q*q*dot(row,an) for row in M)<=md*rd
        value=1+sum((v*aa for v,aa in zip(c,a)),F(0));expected=F(1079,1024) if n==3 else F(13135,12288);assert value==expected
        regression={'status':'PASS old019 example in new matrix only','D_exact':str(value),'not_new_progress':True}
    A=vstack([csr_matrix(-Pu.astype(float)/pd),csr_matrix(Mu.astype(float)*(q*q/md))],format='csr')
    rhs=np.r_[np.zeros(len(Pu)),np.ones(len(Mu))]
    res=linprog(-np.array([float(v) for v in cr]),A_ub=A,b_ub=rhs,bounds=[(None,None)]*len(free),method='highs',options={'primal_feasibility_tolerance':1e-9,'dual_feasibility_tolerance':1e-9})
    tag=f'n{n}_q{q}';base={'n':n,'q':q,'h':str(F(1,q)),'K':K,'orbit_variables':len(points),'free_variables':len(free),'complete_centers':sum(weights),'cell_orbits':len(cells),'all_support_cells':(2*K+4)**n,'full_psi_constraints':len(P),'full_Delta_constraints':len(M),'unique_psi_rows_solved':len(Pu),'unique_Delta_rows_solved':len(Mu),'origin_eliminated_orbit':list(points[elim]),'origin_integer_row':origin,'solver_status':int(res.status),'solver_message':res.message,'old_example_regression':regression,'objective_flux_vs_cell_integral':'PASS exact each orbit','seconds':time.monotonic()-start}
    if not res.success:
        (OUT/f'cb8_result_{tag}.json').write_text(json.dumps(base,indent=2)+'\n');print(json.dumps(base),flush=True);return base
    base['numerical_D']=float(1+np.array([float(v) for v in cr])@res.x);base['numerical_min_psi_Bernstein']=float(np.min(Pr@res.x)/pd);base['numerical_max_Delta_Bernstein']=float(np.max(Mr@res.x)*q*q/md)
    attempts=[];add_indices=set();candidate=None;activities=Pu@res.x/pd
    for attempt,(tol,roundden) in enumerate(((1e-8,10**12),(1e-6,10**14)),1):
        indices=sorted(set(int(i) for i in np.flatnonzero(activities<=tol))|add_indices)
        nums,den,info=exact_face_rebuild([Pu[i].tolist() for i in indices],res.x,roundden)
        pvals=[dot(row,nums) for row in Pu];neg=[i for i,v in enumerate(pvals) if v<0]
        info.update(attempt=attempt,threshold=tol,rounding_denominator=str(roundden),negative_psi_rows=len(neg),minimum_psi_Bernstein_exact=str(F(min([0]+pvals),pd*den)))
        attempts.append(info)
        if neg:add_indices.update(neg);continue
        # Only now repair the Delta constraint by positive scaling.
        rawmax=max([0]+[q*q*dot(row,nums) for row in Mr]);alpha=min(F(1),F(md*den,rawmax)) if rawmax>0 else F(1)
        finalnums=[md*v for v in nums] if alpha<1 else nums;finalden=rawmax if alpha<1 else den
        fullnums=[0]*len(points)
        for j,i in enumerate(free):fullnums[i]=finalnums[j]
        fullnums[elim]=-sum(origin[i]*fullnums[i] for i in free)
        assert sum(v*w for v,w in zip(origin,fullnums))==0
        minpsi=min(dot(row,fullnums) for row in P);maxdelta=max(q*q*dot(row,fullnums) for row in M)
        assert minpsi>=0 and maxdelta<=md*finalden
        aa=[F(v,finalden) for v in fullnums];value=1+sum((v*a for v,a in zip(c,aa)),F(0));rawvalue=1+sum((v*F(num,den) for v,num in zip(cr,nums)),F(0))
        assert value==1+alpha*(rawvalue-1)
        lap0=F(q*q*dot(M[0],fullnums),md*finalden)
        cert={'scope':'continuous signed C2 spline dual feasible point only; not C_n bound','n':n,'q':q,'h':str(F(1,q)),'K':K,'orbits':list(map(list,points)),'orbit_multiplicities':weights,'origin_integer_row':origin,'origin_eliminated_index':elim,'common_coefficient_numerators':list(map(str,fullnums)),'common_coefficient_denominator':str(finalden),'coefficients':list(map(str,aa)),'objective_coefficients':list(map(str,c)),'D_exact':str(value),'exceeds_target':value>2**n,'minimum_psi_Bernstein_exact':str(F(minpsi,pd*finalden)),'maximum_Delta_Bernstein_exact':str(F(maxdelta,md*finalden)),'origin_value_exact':'0','Hessian_diagonal_at_origin_exact':str(lap0/n),'raw_max_Delta_exact':str(F(rawmax,md*den)),'raw_D_exact':str(rawvalue),'positive_scale_alpha':str(alpha),'reconstruction_attempts':attempts,'source_code_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
        (OUT/f'cb8_certificate_{tag}.json').write_text(json.dumps(cert,indent=2)+'\n');candidate=cert;break
    base['reconstruction_attempts']=attempts;base['exact_reconstruction_succeeded']=candidate is not None
    if candidate:base.update(D_exact=candidate['D_exact'],D_float_display=float(F(candidate['D_exact'])),exceeds_target=candidate['exceeds_target'],minimum_psi_Bernstein_exact=candidate['minimum_psi_Bernstein_exact'],maximum_Delta_Bernstein_exact=candidate['maximum_Delta_Bernstein_exact'],Hessian_diagonal_at_origin_exact=candidate['Hessian_diagonal_at_origin_exact'])
    base['seconds']=time.monotonic()-start
    (OUT/f'cb8_result_{tag}.json').write_text(json.dumps(base,indent=2)+'\n');print(json.dumps(base),flush=True);return base
if __name__=='__main__':
    rows=[]
    for n,q,K in ((3,1,3),(3,2,6),(4,1,3)):
        row=solve(n,q,K);rows.append(row)
        if row.get('exceeds_target'):break # prioritize exact above-target certificate before any remaining solve
    (OUT/'cb8_results.json').write_text(json.dumps({'runs':rows,'scope':'three preregistered LPs at most; no method barrier without D>target','source_code_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()},indent=2)+'\n')
