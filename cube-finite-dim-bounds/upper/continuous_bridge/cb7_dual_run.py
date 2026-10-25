"""Two preregistered continuous C2 B-spline LPs; rational repair once per solve.

Finite nonnegative coefficient class only. A large exact D would obstruct
low-cost Laplacian finite-positive-source kernels, not lower-bound C_n.
"""
from pathlib import Path
from fractions import Fraction as F
from itertools import combinations_with_replacement, permutations, product
import json,time,hashlib
import numpy as np
from scipy.sparse import csr_matrix
from scipy.optimize import linprog
OUT=Path(__file__).resolve().parent
BN=np.array([[0,0,0,1],[1,2,4,4],[4,4,2,1],[1,0,0,0]],dtype=np.int64) # /6
DN=np.array([[0,1,2,3],[3,0,-3,-6],[-6,-3,0,3],[3,2,1,0]],dtype=np.int64) # /3
DEN=10**10

def orbit(p):
    return sorted({tuple(a*b for a,b in zip(v,sgn)) for v in set(permutations(p)) for sgn in product((-1,1),repeat=3)})

def model(q,K):
    points=[p for p in combinations_with_replacement(range(K+1),3) if p[-1]>=2]
    lookup={p:i for i,p in enumerate(points)}
    cells=list(combinations_with_replacement(range(K+2),3))
    blocks=[]
    for cell in cells:
        block=np.zeros((64,len(points)),dtype=np.int64)
        choices=[range(max(-K,l-1),min(K,l+2)+1) for l in cell]
        for k in product(*choices):
            p=tuple(sorted(map(abs,k)))
            if p[-1]<2:continue
            mm=[l-v+2 for l,v in zip(cell,k)]
            bx,by,bz=(BN[m] for m in mm);dx,dy,dz=(DN[m] for m in mm)
            val=np.einsum('i,j,k->ijk',dx,by,bz)+np.einsum('i,j,k->ijk',bx,dy,bz)+np.einsum('i,j,k->ijk',bx,by,dz)
            block[:,lookup[p]]+=val.reshape(64)
        blocks.append(block)
    M=np.vstack(blocks)
    def integral(k):
        return sum((F(int(sum(BN[l-k+2])),24) for l in range(-q,q) if -2<=l-k<=1),F(0))
    def prime(t):return F(1,2) if t==-1 else -F(1,2) if t==1 else F(0)
    def flux(k):return prime(q-k)-prime(-q-k)
    c=[];weights=[]
    for p in points:
        oo=orbit(p);weights.append(len(oo));total=F(0)
        for k in oo:
            ints=[integral(v) for v in k];jj=[flux(v) for v in k]
            total+=sum((jj[j]*ints[(j+1)%3]*ints[(j+2)%3] for j in range(3)),F(0))/q
        c.append(-total/8)
    # Independent exact integration route: sum all cube cells' Bernstein averages.
    cc=[F(0) for _ in points]
    for cell,block in zip(cells,blocks):
        if cell[-1]>=q:continue
        mult=8*len(set(permutations(cell)))
        for j in range(len(points)):
            cc[j]-=F(mult*int(block[:,j].sum()),8*q*108*64)
    assert cc==c
    assert sum(weights)==(2*K+1)**3-27
    return points,cells,M,c,weights

def integer_rows(M):
    return [[(int(j),int(row[j])) for j in np.flatnonzero(row)] for row in M]

def solve(q,K):
    start=time.monotonic();points,cells,M,c,weights=model(q,K)
    A=csr_matrix(M.astype(float)*(q*q/108))
    result=linprog(-np.array([float(v) for v in c]),A_ub=A,b_ub=np.ones(len(M)),bounds=(0,None),method='highs',options={'primal_feasibility_tolerance':1e-9,'dual_feasibility_tolerance':1e-9})
    base={'q':q,'h':str(F(1,q)),'K':K,'variables':len(points),'center_count':sum(weights),'cell_orbits':len(cells),'bernstein_constraints':len(M),'all_support_cell_count':(2*K+4)**3,'solver_status':int(result.status),'solver_message':result.message,'seconds':time.monotonic()-start,'objective_crosscheck':'PASS exact boundary flux = all cube-cell Bernstein integrals'}
    if not result.success:
        (OUT/f'cb7_dual_result_q{q}.json').write_text(json.dumps(base,indent=2)+'\n');print(json.dumps(base),flush=True);return base
    rows=integer_rows(M)
    nums=[max(0,int(round(float(v)*DEN))) for v in result.x]
    weighted=[q*q*sum(a*nums[j] for j,a in row) for row in rows]
    raw_max=max(weighted)
    alpha=min(F(1),F(108*DEN,raw_max)) if raw_max>0 else F(1)
    if alpha<1:finalnums=[108*v for v in nums];finalden=raw_max
    else:finalnums=nums;finalden=DEN
    finalmax=max(q*q*sum(a*finalnums[j] for j,a in row) for row in rows)
    assert F(finalmax,108*finalden)<=1 and min(finalnums)>=0
    a=[F(v,finalden) for v in finalnums]
    value=1+sum((v*aa for v,aa in zip(c,a)),F(0))
    rawvalue=1+sum((v*F(nn,DEN) for v,nn in zip(c,nums)),F(0))
    cert={'scope':'continuous C2 B-spline psi certificate; lower bound for admissible Laplacian-kernel cost only, not for C_n','dimension':3,'q':q,'h':str(F(1,q)),'K':K,'orbits':list(map(list,points)),'orbit_multiplicities':weights,'coefficients':list(map(str,a)),'common_coefficient_numerators':list(map(str,finalnums)),'common_coefficient_denominator':str(finalden),'exact_objective_coefficients':list(map(str,c)),'D_exact':str(value),'strictly_above_8':value>8,'origin_value_exact':'0','maximum_Bernstein_coefficient_exact':str(F(finalmax,108*finalden)),'complete_support_cells':'ell_j=-K-2,...,K+1; all 64 coefficients per reflection/permutation representative','repair':'rounded nonnegative to 1e-10 once; globally scaled at most once','scale_alpha':str(alpha),'raw_rounded_D_exact':str(rawvalue),'source_code_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    (OUT/f'cb7_dual_certificate_q{q}.json').write_text(json.dumps(cert,indent=2)+'\n')
    base.update(numerical_D=float(1+np.array([float(v) for v in c])@result.x),numerical_max_Bernstein=float(np.max(A@result.x)),rational_D_exact=str(value),rational_D_float_display=float(value),strictly_above_8=value>8,raw_max_Bernstein_exact=str(F(raw_max,108*DEN)),scale_alpha=str(alpha),repaired_max_Bernstein_exact=str(F(finalmax,108*finalden)),positive_orbit_coefficients=sum(aa>0 for aa in a),seconds=time.monotonic()-start)
    (OUT/f'cb7_dual_result_q{q}.json').write_text(json.dumps(base,indent=2)+'\n');print(json.dumps(base),flush=True)
    return base
if __name__=='__main__':
    results=[solve(1,3),solve(2,6)]
    (OUT/'cb7_dual_results.json').write_text(json.dumps({'runs':results,'scope':'two finite continuous psi classes only; no new C_n bound','source_code_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()},indent=2)+'\n')
