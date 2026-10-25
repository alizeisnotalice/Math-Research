"""Independent signed-spline audit via truncated powers and scalar integer DP."""
from pathlib import Path
from itertools import product
from fractions import Fraction as F
from math import comb
import sympy as s
import json,hashlib,argparse,time
HERE=Path(__file__).resolve().parent
HASHES={'proofs/signed_dual_interface.md':'6f1501b2850670e2d4b1a61c2339a470cd6ea52e89681e393ebf96c535c4c2e1','src/check_signed_dual_interface.py':'9683d6b91e0b5ec181f28f34a0c6ad424ad830b35a5f71876d5b675266fa9391','results/round007_signed_interface.json':'6128935e6522d16d78460006a1960b222594cc14f33a8e6b96396d7b49cad03f'}
u=s.symbols('u');DEN=18

def cardinal(cell,shift=0):
    # B(t)=1/6 sum (-1)^i binom(4,i) (t+2-i)_+^3.
    return s.expand(sum(((-1)**i*comb(4,i)*(u+cell-shift+2-i)**3 for i in range(5) if cell-shift+2-i>=0),s.Integer(0))/6)
def bern(poly):
    p=s.Poly(poly,u);v=[sum((F(p.nth(j))*F(comb(i,j),comb(3,j)) for j in range(i+1)),F(0)) for i in range(4)]
    assert all((DEN*x).denominator==1 for x in v)
    return [int(DEN*x) for x in v]
def run(root):
    for name,expected in HASHES.items():assert hashlib.sha256((root/name).read_bytes()).hexdigest()==expected
    author=json.loads((root/'results/round007_signed_interface.json').read_text())['results']
    polynomials={};states=[]
    for cell in range(-3,3):
        b=cardinal(cell);phi=s.expand(s.Rational(2,3)*(cardinal(cell,1)+cardinal(cell,-1))-b/3)
        polynomials[cell]=(b,phi);bs=bern(b);bd=bern(s.diff(b,u,2));fs=bern(phi);fd=bern(s.diff(phi,u,2));assert min(bs)>=0 and min(fs)>=0
        for i in range(4):states.append((cell,bs[i],bd[i],fs[i],fd[i]))
    # C2 joins of both independently generated functions, and zero extension.
    for which in (0,1):
      for order in range(3):
        for cell in range(-3,2):assert s.diff(polynomials[cell][which],u,order).subs(u,1)==s.diff(polynomials[cell+1][which],u,order).subs(u,0)
        assert s.diff(polynomials[-3][which],u,order).subs(u,0)==s.diff(polynomials[2][which],u,order).subs(u,1)==0
    phi0=polynomials[0][1]
    assert s.expand(phi0-(u*u-s.Rational(11,18)*u**3))==0
    assert phi0.subs(u,0)==s.diff(phi0,u).subs(u,0)==0 and s.diff(phi0,u,2).subs(u,0)==2
    # Direct univariate volume integrals and independent boundary derivative jumps.
    Ib=sum((F(s.integrate(polynomials[c][0],(u,0,1))) for c in (-1,0)),F(0));If=sum((F(s.integrate(polynomials[c][1],(u,0,1))) for c in (-1,0)),F(0))
    Ibd=sum((F(s.integrate(s.diff(polynomials[c][0],u,2),(u,0,1))) for c in (-1,0)),F(0));Ifd=sum((F(s.integrate(s.diff(polynomials[c][1],u,2),(u,0,1))) for c in (-1,0)),F(0))
    for which,integ in ((0,Ibd),(1,Ifd)):
      jump=F(s.diff(polynomials[1][which],u).subs(u,0)-s.diff(polynomials[-1][which],u).subs(u,0));assert jump==integ
    assert (Ib,If,Ibd,Ifd)==(F(11,12),F(13,36),F(-1),F(1,3))
    runs=[]
    for n,expectedC,scale,expectedD in [(3,F(8,3),F(3,8),F(1079,1024)),(4,F(64,27),F(27,64),F(13135,12288))]:
      start=time.monotonic();minpsi=0;maxlap=0;minlap=0;sumcube=0;count=0
      # Coefficient of t and t*v in product_j(B_j+t phi_j+v B''_j+t*v phi''_j).
      # Sequential scalar DP is independent of the author's NumPy outer products.
      for items in product(states,repeat=n):
        z00,z10,z01,z11=1,0,0,0
        for cell,b,bd,f,fd in items:
          z00,z10,z01,z11=z00*b,z10*b+z00*f,z01*b+z00*bd,z11*b+z10*bd+z01*f+z00*fd
        minpsi=min(minpsi,z10);maxlap=max(maxlap,z11);minlap=min(minlap,z11);count+=1
        if all(cell in (-1,0) for cell,*_ in items):sumcube+=z11
      assert count==6**n*4**n and minpsi>=0
      C=F(maxlap,DEN**n);assert C==expectedC and scale*C==1
      integral=F(sumcube,DEN**n*4**n)
      direct=n*Ifd*Ib**(n-1)+n*(n-1)*If*Ibd*Ib**(n-2);assert integral==direct
      hess=scale*2*F(2,3)**(n-1);assert hess==F(1,n)
      D=1-scale*direct/F(2**n);assert D==expectedD and D<2**n and D>1
      a=next(x for x in author if x['n']==n)
      for field,value in [('full_cells',6**n),('coefficients_per_polynomial',count),('psi_Bernstein_min',str(F(minpsi,DEN**n))),('raw_Delta_Bernstein_max',str(C)),('raw_Delta_Bernstein_min',str(F(minlap,DEN**n))),('scale',str(scale)),('scaled_origin_Hessian_diagonal',str(hess)),('cube_Delta_integral_raw',str(direct)),('D_exact',str(D))]:assert a[field]==value
      runs.append({'n':n,'status':'PASS','full_cells':6**n,'coefficients_per_polynomial':count,'nonnegative_psi_Bernstein_min':str(F(minpsi,DEN**n)),'raw_Delta_Bernstein_max':str(C),'raw_Delta_Bernstein_min':str(F(minlap,DEN**n)),'scale':str(scale),'scaled_origin_Hessian_diagonal':str(hess),'cube_Delta_integral_raw':str(direct),'D_exact':str(D),'D_display_only':float(D),'nonzero_signed_center_coefficient':str(-scale*F(n,3)),'signed_axis_neighbor_coefficients':str(scale*F(2,3)),'seconds':time.monotonic()-start})
    return {'status':'PASS independent signed-interface audit','source_sha256':HASHES,'univariate_integrals':{'B':str(Ib),'phi':str(If),'B_second':str(Ibd),'phi_second':str(Ifd)},'runs':runs,'shared_libraries':{'Fraction':'standard Python','SymPy':s.__version__},'independence':'truncated-power generation and arbitrary-integer scalar coefficient recurrence; no author imports, NumPy, solver or new LP'}
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--project',type=Path,default=HERE.parent);ap.add_argument('--save',action='store_true');args=ap.parse_args();out=run(args.project);out['checker_sha256']=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    if args.save:HERE.joinpath('019_results.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out,indent=2))
