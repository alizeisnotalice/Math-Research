"""Independent full signed tensor certificate check, exact object-integer transforms."""
from pathlib import Path
from fractions import Fraction as F
from itertools import product
from collections import Counter
from math import comb,prod
from functools import lru_cache
import numpy as np
import sympy as s
import json,hashlib,argparse,time,gc
HERE=Path(__file__).resolve().parent
EXPECTED={(3,1):'a36b9f52081528685c866d3a2714d0ca6c2be3580752e1394d281bbd9287ebaf',(3,2):'67f2c0b6644ef252846fbd4d4f761619d7afc21ae3e6bba76bf0f2cb181e42a8',(4,1):'46ee76cfd5f0a9cdf7979e699e77bb7c0eb88c4343f4d1029628af72129df2eb'}
u=s.symbols('u')
def cardinal(cell):
    return s.expand(sum(((-1)**i*comb(4,i)*(u+cell+2-i)**3 for i in range(5) if cell+2-i>=0),s.Integer(0))/6)
def table():
    b={};dd={};powers={}
    for m in range(-2,2):
      p=cardinal(m);P=s.Poly(p,u);powers[m]=[F(P.nth(j)) for j in range(4)]
      for out,poly in [(b,p),(dd,s.diff(p,u,2))]:
        P=s.Poly(poly,u);v=[6*sum((F(P.nth(j))*F(comb(i,j),comb(3,j)) for j in range(i+1)),F(0)) for i in range(4)]
        assert all(x.denominator==1 for x in v);out[m]=list(map(int,v))
    return b,dd,powers
B,DD,P=table()
def operators(K,tab):
    rows=[]
    for ell in range(-K-2,K+2):
      for i in range(4):rows.append([(k+K,tab[ell-k][i]) for k in range(-K,K+1) if ell-k in tab and tab[ell-k][i]])
    return rows
def apply_axis(arr,axis,rows):
    assert arr.dtype==object
    moved=np.moveaxis(arr,axis,-1);out=np.zeros(moved.shape[:-1]+(len(rows),),dtype=object)
    for i,row in enumerate(rows):
      for j,w in row:out[...,i]+=moved[...,j]*w
    return np.moveaxis(out,-1,axis)
@lru_cache(None)
def integ(k,q):return sum((p[j]/(j+1) for ell in range(-q,q) if (p:=P.get(ell-k)) is not None for j in range(4)),F(0))
def prime(t):return P[t][1] if t in P else F(0)
def flux(k,q):
    n=len(k);v=[integ(z,q) for z in k];jump=[prime(q-z)-prime(-q-z) for z in k]
    return F(1,q**(n-2))*sum((jump[j]*prod(v[l] for l in range(n) if l!=j) for j in range(n)),F(0))
def origin_basis(k,which=None):
    out=F(1)
    for j,z in enumerate(k):
      value=(F(-2) if z==0 else F(1) if abs(z)==1 else F(0)) if j==which else (F(2,3) if z==0 else F(1,6) if abs(z)==1 else F(0))
      out*=value
    return out
def verify(path,n,q):
    started=time.monotonic();sha=hashlib.sha256(path.read_bytes()).hexdigest();assert sha==EXPECTED[n,q];r=json.loads(path.read_text());K=r['K'];den=int(r['common_coefficient_denominator']);assert r['n']==n and r['q']==q and F(r['h'])==F(1,q)
    lengths=[len(r[k]) for k in ('orbits','orbit_multiplicities','common_coefficient_numerators','coefficients','objective_coefficients','origin_integer_row')];assert len(set(lengths))==1
    look={};coeff=[F(x) for x in r['coefficients']];nums=list(map(int,r['common_coefficient_numerators']))
    for i,k in enumerate(r['orbits']):
      assert len(k)==n and all(type(x)is int for x in k) and list(sorted(k))==k and min(k)>=0 and max(k)<=K and tuple(k) not in look
      look[tuple(k)]=i;assert coeff[i]==F(nums[i],den)
    assert len(look)==comb(K+n,n)
    A=np.zeros((2*K+1,)*n,dtype=object);counts=Counter();objective=[F(0)]*len(look);originrow=[F(0)]*len(look);hess=[F(0)]*n;origin=F(0);nonzero=0
    for k in product(range(-K,K+1),repeat=n):
      i=look[tuple(sorted(map(abs,k)))];nn=nums[i];A[tuple(z+K for z in k)]=nn;counts[i]+=1;nonzero+=nn!=0
      objective[i]-=flux(k,q)/F(2**n);ov=origin_basis(k);originrow[i]+=ov;origin+=coeff[i]*ov
      for j in range(n):hess[j]+=q*q*coeff[i]*origin_basis(k,j)
    assert [counts[i] for i in range(len(look))]==r['orbit_multiplicities']
    assert list(map(str,objective))==r['objective_coefficients'];assert [int(v*3**n) for v in originrow]==r['origin_integer_row']
    elim=r['origin_eliminated_index'];assert r['orbits'][elim]==[1]*n and r['origin_integer_row'][elim]==1
    assert origin==0 and all(z==F(r['Hessian_diagonal_at_origin_exact']) for z in hess) and 0<=n*hess[0]<=1
    psi=A;lap=np.zeros(A.shape,dtype=object);opB=operators(K,B);opDD=operators(K,DD)
    for axis in range(n):
      newlap=apply_axis(lap,axis,opB);addition=apply_axis(psi,axis,opDD);newlap+=addition;del addition
      newpsi=apply_axis(psi,axis,opB);psi,lap=newpsi,newlap;gc.collect()
    assert psi.size==lap.size==(4*(2*K+4))**n
    lo=int(psi.min());hi=int(lap.max());lolo=int(lap.min());assert lo>=0 and hi*q*q<=den*6**n
    assert F(lo,den*6**n)==F(r['minimum_psi_Bernstein_exact']);assert F(hi*q*q,den*6**n)==F(r['maximum_Delta_Bernstein_exact'])
    take=tuple(slice(4*(K+2-q),4*(K+2+q)) for _ in range(n));integral=F(q*q*sum(lap[take].flat),den*6**n*4**n*q**n)
    fluxsum=-F(2**n)*sum((a*b for a,b in zip(coeff,objective)),F(0));assert integral==fluxsum
    D=1-integral/F(2**n);assert D==F(r['D_exact']) and D<2**n
    alpha=F(r['positive_scale_alpha']);assert 0<alpha<=1;rawC=F(hi*q*q,den*6**n)/alpha;rawD=1+(D-1)/alpha
    assert rawC==F(r['raw_max_Delta_exact']) and rawD==F(r['raw_D_exact']) and alpha==(1/rawC if rawC>1 else F(1));assert D==1+alpha*(rawD-1)
    del psi,lap,newpsi,newlap,A;gc.collect()
    return {'n':n,'q':q,'status':'PASS exact independent full tensor audit','certificate_sha256':sha,'orbits':len(look),'full_centers':(2*K+1)**n,'nonzero_centers':nonzero,'full_cells':(2*K+4)**n,'coefficients_per_polynomial':(4*(2*K+4))**n,'minimum_psi_Bernstein_exact':str(F(lo,den*6**n)),'maximum_Delta_Bernstein_exact':str(F(hi*q*q,den*6**n)),'minimum_Delta_Bernstein_exact':str(F(lolo*q*q,den*6**n)),'origin_exact':str(origin),'Hessian_diagonal_physical_exact':str(hess[0]),'cube_Delta_integral_power_flux_and_Bernstein':str(integral),'all_orbit_objective_coefficients_verified':True,'origin_row_and_all_ones_elimination_verified':True,'D_exact':str(D),'D_display_only':float(D),'raw_max_Delta_exact':str(rawC),'raw_D_exact':str(rawD),'positive_scale_alpha':str(alpha),'repair_formula_verified':True,'seconds':time.monotonic()-started}
def run(root):
    return {'status':'PASS signed continuous feasible dual certificates','runs':[verify(root/f'cb8_certificate_n{n}_q{q}.json',n,q) for n,q in ((3,1),(3,2),(4,1))],'independence':'independent abs-sort full lattice; truncated-power generated exact one-dimensional matrices; full tensor transforms for psi and Laplacian; no author LP/checker imports or symmetry-reduced cells','arithmetic':'NumPy object arrays contain Python arbitrary-precision integers only; no float/int64 coefficient operations','libraries':{'numpy':np.__version__,'sympy':s.__version__,'Fraction':'Python standard library'},'scope':'three saved continuous feasible psi only, no C_n lower bound, LP optimality or whole-method obstruction'}
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--source-dir',type=Path,default=HERE.parent/'upper/continuous_bridge');ap.add_argument('--save',action='store_true');args=ap.parse_args();out=run(args.source_dir);out['checker_sha256']=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    if args.save:HERE.joinpath('021_results.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out,indent=2))
