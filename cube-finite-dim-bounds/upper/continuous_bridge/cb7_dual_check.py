"""Stdlib exact all-cell verification of saved continuous spline certificates.
Rebuilds full signed lattice, without NumPy, SciPy or reduced LP matrix.
Shares the explicit B-spline formula and integer/Fraction arithmetic.
"""
from pathlib import Path
from fractions import Fraction as F
from itertools import permutations,product
import json,time
OUT=Path(__file__).resolve().parent
# Exact power polynomials for B(t+m), 0<=t<=1, numerators /6.
POWER={-2:[0,0,0,1],-1:[1,3,3,-3],0:[4,0,-6,3],1:[1,-3,3,-1]}
# Basis Bernstein coefficients /6, second derivatives elevated to degree3 /3.
B={-2:[0,0,0,1],-1:[1,2,4,4],0:[4,4,2,1],1:[1,0,0,0]}
DD={-2:[0,1,2,3],-1:[3,0,-3,-6],0:[-6,-3,0,3],1:[3,2,1,0]}

def derivative_value(p,order,t):
    return sum((F(p[i],6)*F(__import__('math').factorial(i),__import__('math').factorial(i-order))*t**(i-order) for i in range(order,len(p))),F(0))

def basis_checks():
    for order in range(3):
        assert derivative_value(POWER[-2],order,F(0))==0
        assert derivative_value(POWER[1],order,F(1))==0
        for m in (-2,-1,0):
            assert derivative_value(POWER[m],order,F(1))==derivative_value(POWER[m+1],order,F(0))
    # Re-derive degree3 Bernstein tables from power coefficients.
    from math import comb
    for m,p in POWER.items():
        bern=[sum((F(p[j],6)*F(comb(i,j),comb(3,j)) for j in range(i+1)),F(0)) for i in range(4)]
        assert bern==[F(v,6) for v in B[m]] and min(bern)>=0
        sec=[F((j+2)*(j+1)*p[j+2],6) for j in range(2)]
        elevated=[sec[0]+F(i,3)*sec[1] for i in range(4)]
        assert elevated==[F(v,3) for v in DD[m]]

def verify(q):
    start=time.monotonic();rec=json.loads((OUT/f'cb7_dual_certificate_q{q}.json').read_text())
    assert rec['q']==q and rec['dimension']==3
    K=rec['K'];den=int(rec['common_coefficient_denominator']);full={}
    for p,nn,a,mult in zip(rec['orbits'],rec['common_coefficient_numerators'],rec['coefficients'],rec['orbit_multiplicities']):
        nn=int(nn);assert nn>=0 and F(nn,den)==F(a)
        assert max(p)>=2
        centers={tuple(vv*ss for vv,ss in zip(v,sgn)) for v in set(permutations(p)) for sgn in product((-1,1),repeat=3)}
        assert len(centers)==mult
        for k in centers:
            assert k not in full;full[k]=nn
    assert len(full)==(2*K+1)**3-27
    assert sum((F(nn,den)*__import__('functools').reduce(lambda a,b:a*b,[F(2,3) if v==0 else F(1,6) if abs(v)==1 else F(0) for v in k],F(1)) for k,nn in full.items()),F(0))==0
    # Full support box, without using symmetry of cells or reduced matrix.
    maximum=0;maximum_cell=None;cube_integral_num=0;count=0
    for cell in product(range(-K-2,K+2),repeat=3):
        values=[0]*64
        choices=[range(max(-K,l-1),min(K,l+2)+1) for l in cell]
        for k in product(*choices):
            nn=full.get(k,0)
            if not nn:continue
            m=[l-v for l,v in zip(cell,k)];bb=[B[v] for v in m];dd=[DD[v] for v in m]
            for i,j,l in product(range(4),repeat=3):
                val=dd[0][i]*bb[1][j]*bb[2][l]+bb[0][i]*dd[1][j]*bb[2][l]+bb[0][i]*bb[1][j]*dd[2][l]
                values[16*i+4*j+l]+=q*q*nn*val
        mm=max(values);count+=1
        if mm>maximum:maximum=mm;maximum_cell=cell
        assert mm<=108*den
        if all(-q<=v<q for v in cell):cube_integral_num+=sum(values)
    assert count==(2*K+4)**3
    by_cells=F(cube_integral_num,108*den*64*q**3)
    # Independent objective route: exact integrals of power polynomials and boundary B'.
    def integral(k):
        return sum((F(p[j],6*(j+1)) for ell in range(-q,q) if (p:=POWER.get(ell-k)) is not None for j in range(4)),F(0))
    def prime(t):
        if -2<=t<=1:return derivative_value(POWER[t],1,F(0))
        return F(0)
    by_flux=F(0)
    for k,nn in full.items():
        if not nn:continue
        ints=[integral(v) for v in k];jump=[prime(q-v)-prime(-q-v) for v in k]
        by_flux+=F(nn,den*q)*sum((jump[j]*ints[(j+1)%3]*ints[(j+2)%3] for j in range(3)),F(0))
    assert by_cells==by_flux
    value=1-by_flux/8
    assert value==F(rec['D_exact'])
    assert F(maximum,108*den)==F(rec['maximum_Bernstein_coefficient_exact'])
    return {'q':q,'K':K,'status':'PASS full signed-lattice exact checker','all_cells_checked':count,'Bernstein_coefficients_checked':64*count,'maximum_Bernstein_exact':str(F(maximum,108*den)),'maximum_cell_example':maximum_cell,'origin_value_exact':'0','cube_Delta_integral_exact':str(by_flux),'D_exact':str(value),'strictly_above_8':value>8,'seconds':time.monotonic()-start,'independence':'does not import discovery or reduced matrix; shares B-spline formula and Fraction arithmetic, not a wholly independent audit'}
if __name__=='__main__':
    basis_checks();rows=[verify(1),verify(2)]
    out={'C2_nonnegative_basis_check':'PASS power-polynomial knot derivatives and Bernstein tables','runs':rows};(OUT/'cb7_dual_check_results.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
