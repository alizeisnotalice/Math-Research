"""Saved certificates only: independent orbit and power-first all-cell audit."""
from pathlib import Path
from fractions import Fraction as F
from itertools import product
from collections import Counter,defaultdict
from functools import lru_cache
from math import comb
import hashlib,json,argparse,time
import sympy as s
HERE=Path(__file__).resolve().parent
EXPECTED={1:'b9dad20dd2c7408add7e95e6078a9e322db17127624e64fba206c072aed86c21',2:'dd4e8b98bab49e6c15a10ad2063ef5119a3557d7edb7dbe1ca5968e9e543875e'}
AUTHOR_SHA='89cbeab8bcd5733fcf640c6e0e979d5bf7d6b26d27bf1fa9cc1c4b9b5a161e75'
IDX=list(product(range(4),repeat=3))
flat=lambda i,j,k:16*i+4*j+k

def basis_power():
    # Derive each cell from the centered B definition; no imported/hardcoded tables.
    t,u=s.symbols('t u');polys=[(2+t)**3/6,(4-6*t*t-3*t**3)/6,(4-6*t*t+3*t**3)/6,(2-t)**3/6];table={}
    for a,p in zip(range(-2,2),polys):
        poly=s.Poly(s.expand(6*p.subs(t,a+u)),u);table[a]=[int(poly.nth(j)) for j in range(4)]
    for order in range(3):
        for j in range(3):assert s.diff(polys[j],t,order).subs(t,j-1)==s.diff(polys[j+1],t,order).subs(t,j-1)
        assert s.diff(polys[0],t,order).subs(t,-2)==s.diff(polys[-1],t,order).subs(t,2)==0
    return table
P=basis_power()
SECOND={a:[2*p[2],6*p[3],0,0] for a,p in P.items()}
# 6 times the power->cubic-Bernstein conversion matrix, all integer.
M=[[int(6*F(comb(i,j),comb(3,j))) if j<=i else 0 for j in range(4)] for i in range(4)]
def tensor_bernstein_numerator(v):
    for axis in range(3):
        out=[0]*64
        for ijk in IDX:
            r=ijk[axis];acc=0
            for j in range(r+1):
                at=list(ijk);at[axis]=j;acc+=M[r][j]*v[flat(*at)]
            out[flat(*ijk)]=acc
        v=out
    return v
@lru_cache(None)
def local_power(m):
    out=[]
    for ijk in IDX:
        value=sum(SECOND[m[j]][ijk[j]]*P[m[(j+1)%3]][ijk[(j+1)%3]]*P[m[(j+2)%3]][ijk[(j+2)%3]] for j in range(3))
        if value:out.append((flat(*ijk),value))
    return out
@lru_cache(None)
def basis_integral(k,q):
    return sum((F(p[j],6*(j+1)) for a in range(-q-k,q-k) if (p:=P.get(a)) is not None for j in range(4)),F(0))
def bprime(t):
    return F(P[t][1],6) if t in P else F(0)
def flux_coefficient(k,q):
    vals=[basis_integral(z,q) for z in k];jumps=[bprime(q-z)-bprime(-q-z) for z in k]
    return F(1,q)*sum((jumps[j]*vals[(j+1)%3]*vals[(j+2)%3] for j in range(3)),F(0))
def verify(path,q):
    started=time.monotonic();sha=hashlib.sha256(path.read_bytes()).hexdigest();assert sha==EXPECTED[q];r=json.loads(path.read_text());K=r['K'];den=int(r['common_coefficient_denominator'])
    assert r['dimension']==3 and r['q']==q and F(r['h'])==F(1,q) and den>0
    fields=['orbits','orbit_multiplicities','coefficients','common_coefficient_numerators','exact_objective_coefficients'];assert len({len(r[k]) for k in fields})==1
    orbit={}
    for i,k in enumerate(r['orbits']):
        assert all(type(v) is int for v in k) and tuple(k)==tuple(sorted(k)) and 0<=min(k)<=max(k)<=K and max(k)>=2
        assert tuple(k) not in orbit
        val=F(r['coefficients'][i]);assert val>=0 and val==F(int(r['common_coefficient_numerators'][i]),den)
        orbit[tuple(k)]=i
    full={};multiplicity=Counter();objective=[F(0)]*len(orbit);origin=F(0)
    for k in product(range(-K,K+1),repeat=3):
        canonical=tuple(sorted(map(abs,k)))
        if max(canonical)<2:assert canonical not in orbit;continue
        i=orbit[canonical];multiplicity[i]+=1;nn=int(r['common_coefficient_numerators'][i]);full[k]=nn
        objective[i]-=flux_coefficient(k,q)/8
        origin+=F(nn,den)*__import__('math').prod(F(2,3) if v==0 else F(1,6) if abs(v)==1 else F(0) for v in k)
    assert len(full)==(2*K+1)**3-27 and origin==0
    assert [multiplicity[i] for i in range(len(orbit))]==r['orbit_multiplicities']
    assert list(map(str,objective))==r['exact_objective_coefficients']
    cells=defaultdict(lambda:[0]*64)
    for k,nn in full.items():
        if not nn:continue
        for m in product(range(-2,2),repeat=3):
            cell=tuple(k[j]+m[j] for j in range(3));a=cells[cell]
            for idx,v in local_power(m):a[idx]+=q*q*nn*v
    # We explicitly loop the complete rectangular support, including empty cells.
    maximum=0;minval=0;cube_power_integral_num=0;nonzero_cells=0
    for cell in product(range(-K-2,K+2),repeat=3):
        a=cells.get(cell)
        if a is None:continue
        nonzero_cells+=1;b=tensor_bernstein_numerator(a);maximum=max(maximum,max(b));minval=min(minval,min(b));assert max(b)<=46656*den
        if all(-q<=ell<q for ell in cell):
            cube_power_integral_num+=sum(a[flat(i,j,k)]*(12//(i+1))*(12//(j+1))*(12//(k+1)) for i,j,k in IDX)
    cube_integral=F(cube_power_integral_num,216*den*12**3*q**3)
    flux=sum((F(nn,den)*flux_coefficient(k,q) for k,nn in full.items()),F(0));assert cube_integral==flux
    D=1-cube_integral/8;assert str(D)==r['D_exact'];maxbern=F(maximum,46656*den);assert str(maxbern)==r['maximum_Bernstein_coefficient_exact']
    alpha=F(r['scale_alpha']);assert 0<alpha<=1
    raw=[F(z)/alpha for z in r['coefficients']];assert all((z*10**10).denominator==1 for z in raw)
    Draw=1+sum((a*b for a,b in zip(raw,objective)),F(0));assert str(Draw)==r['raw_rounded_D_exact']
    C=maxbern/alpha;assert alpha==(1/C if C>1 else F(1));assert D==1+alpha*(Draw-1)
    if C>1:assert D==1+(Draw-1)/C
    assert D<8
    return {'q':q,'certificate_sha256':sha,'status':'PASS independent power-first all-cell audit','orbit_count':len(orbit),'full_lattice_centers':len(full),'nonzero_lattice_centers':sum(nn>0 for nn in full.values()),'all_support_cells_checked':(2*K+4)**3,'tensor_coefficients_checked':64*(2*K+4)**3,'cells_accumulated_from_nonzero_terms':nonzero_cells,'origin_exact':str(origin),'maximum_Bernstein_exact':str(maxbern),'minimum_Bernstein_exact':str(F(minval,46656*den)),'cube_integral_by_power_and_flux':str(cube_integral),'D_exact':str(D),'D_display_only':float(D),'raw_rounded_D_exact':str(Draw),'raw_Bernstein_max_exact':str(C),'scale_alpha_exact':str(alpha),'repair_formula_verified':True,'raw_coefficients_on_1e_10_grid':True,'all_orbit_objective_coefficients_verified':True,'seconds':time.monotonic()-started}
def main(root):
    assert hashlib.sha256((root/'cb7_dual_check.py').read_bytes()).hexdigest()==AUTHOR_SHA
    return {'status':'PASS','runs':[verify(root/f'cb7_dual_certificate_q{q}.json',q) for q in (1,2)],'author_checker_sha256_read_only':AUTHOR_SHA,'implementation':'independent abs-sort orbit enumeration; generated power polynomials; power-first cell sum then 3-axis Bernstein conversion; direct power integral plus boundary flux; no author imports/LP solve','shared':'same mathematical B-spline definition and Python int/Fraction; SymPy '+s.__version__}
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--source-dir',type=Path,default=HERE.parent/'upper/continuous_bridge');ap.add_argument('--save',action='store_true');args=ap.parse_args();out=main(args.source_dir);out['checker_sha256']=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    if args.save:HERE.joinpath('018_results.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out,indent=2))
