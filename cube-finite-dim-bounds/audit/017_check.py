"""Independent finite exact audit: no author spline/LP/cost checker imports."""
from fractions import Fraction as F
from itertools import product
from math import comb
from pathlib import Path
import json,hashlib,argparse
import sympy as sp
HERE=Path(__file__).resolve().parent
SOURCE=HERE.parent/'proofs/continuous_laplacian_dual.md'
COST=HERE.parent/'upper/u06_parameter_check.py'
def bern3(poly,u):
    p=sp.Poly(poly,u); a=[F(p.nth(k)) for k in range(4)]
    return [sum((a[k]*F(comb(i,k),comb(3,k)) for k in range(i+1)),F(0)) for i in range(4)]
def mul(xs):
    v=F(1)
    for x in xs:v*=x
    return v
def run(source=SOURCE,cost_source=COST):
    sha=hashlib.sha256(source.read_bytes()).hexdigest();assert sha=='d87178e5f35374eb3fd51ebee976a2aa7e3b901b2c98a29807e493e498a54af8'
    costsha=hashlib.sha256(cost_source.read_bytes()).hexdigest();assert costsha=='5f41dabdc023ea00c5213e1e812e89bf797fabdfc57f8e722f58793fd2fe9707'
    u,t=sp.symbols('u t');parts=[(t+2)**3/6,(4-6*t*t-3*t**3)/6,(4-6*t*t+3*t**3)/6,(2-t)**3/6]
    # Check C2 joins, endpoint vanishing, and unit integral independently.
    for j in range(3):
      at=j-1
      for d in range(3):assert sp.diff(parts[j],t,d).subs(t,at)==sp.diff(parts[j+1],t,d).subs(t,at)
    for p,at in [(parts[0],-2),(parts[-1],2)]:
      for d in range(3):assert sp.diff(p,t,d).subs(t,at)==0
    assert sum(sp.integrate(p,(t,k,k+1)) for p,k in zip(parts,range(-2,2)))==1
    one=[]
    for p,k in zip(parts,range(-2,2)):
      v=bern3(sp.expand(p.subs(t,k+u)),u);w=bern3(sp.expand(sp.diff(p,t,2).subs(t,k+u)),u)
      assert min(v)>=0
      one.append((v,w))
    rows=[];h=F(2,5)
    for n in (3,4):
      # Single tensor coefficient at k=(2,0,...,0); B(-2)=0 exactly.
      c=h*h/F(4*n);maximum=None;minimum=None;count=0
      for cells in product(range(4),repeat=n):
        for bs in product(range(4),repeat=n):
          val=c/h**2*sum((one[cells[j]][1][bs[j]]*mul(one[cells[k]][0][bs[k]] for k in range(n) if k!=j) for j in range(n)),F(0))
          maximum=val if maximum is None else max(maximum,val);minimum=val if minimum is None else min(minimum,val);count+=1
      assert maximum<=1
      # Only x1 faces contribute: t at x1=1 is 1/2, derivative -5/8;
      # others' full B integrals =h and their exterior boundary slopes vanish.
      normal_flux=-c*F(5,8)*h**(n-2);D=1-normal_flux/F(2**n)
      assert D==(F(2401,2400) if n==3 else F(16001,16000))
      # Separate piecewise polynomial volume integral, not a flux reuse.
      integrals=[];second_integrals=[]
      for k in [2]+[0]*(n-1):
        bounds=(F(-1)/h-k,F(1)/h-k);iv=F(0);iw=F(0)
        for piece,a in zip(parts,range(-2,2)):
          lo=max(F(a),bounds[0]);hi=min(F(a+1),bounds[1])
          if lo<hi:
            iv+=h*F(sp.integrate(piece,(t,lo,hi)))
            iw+=F(sp.integrate(sp.diff(piece,t,2),(t,lo,hi)))/h
        integrals.append(iv);second_integrals.append(iw)
      direct=c*sum((second_integrals[j]*mul(integrals[k] for k in range(n) if k!=j) for j in range(n)),F(0))
      assert direct==normal_flux
      # Grid-only false certificate: center .5 in every coordinate, h=.1,
      # support (.3,.7)^n; Delta is zero at every integer lattice point.
      # At (.65,.5,...,.5): B''(1.5)=.5, B(1.5)=1/48, B(0)=2/3, B''(0)=-2.
      bad=100*(F(1,2)*F(2,3)**(n-1)-(n-1)*F(2)*F(1,48)*F(2,3)**(n-2))
      assert bad==(F(50,3) if n==3 else F(250,27)) and bad>1
      rows.append({'n':n,'h':str(h),'coefficient':str(c),'all_support_cells':4**n,'tensor_bernstein_coefficients_checked':count,'laplacian_bernstein_max':str(maximum),'laplacian_bernstein_min':str(minimum),'dual_value_exact':str(D),'dual_value_display':float(D),'grid_only_false_certificate_laplacian':str(bad)})
    costs=[]
    for n,A,b4,b6,outer in [(3,F('.4089479272201547'),F('.005004477192161995'),F('.0009386738283819751'),F(16,5)),(4,F('.18157763153667977'),F('.004'),F('.0005'),F(18,5))]:
      # Independent simplex enumeration of the H6 maximum used below.
      beta=F(15,n+8);gamma=F(30,(n+8)*(n+4));vals=[F(1,k*k)-beta/F(k)+gamma for k in range(1,n+1)];T=2*beta/3
      for p in range(1,n):
        for q in range(1,n-p+1):
          if p!=q:
            a=(1-q*T)/(p-q);b=T-a
            if a>0 and b>0: vals.append(p*a**3+q*b**3-beta*(p*a*a+q*b*b)+gamma)
      assert max(vals)==(F(32,693) if n==3 else F(1,16))
      yy=sp.symbols('z:'+str(n));rr=sum(v*v for v in yy);H6=sum(v**6 for v in yy)-sp.Rational(15,n+8)*rr*sum(v**4 for v in yy)+sp.Rational(30,(n+8)*(n+4))*rr**3
      assert sp.expand(sum(sp.diff(H6,v,2) for v in yy))==0
      r=F(21,10);p4=F(n-1,n+2);p6=F(32,693) if n==3 else F(1,16)
      # Independent positive-core lower and harmonic-mean integral.
      margin=r**(2-n)-A-b4*p4*r**4-b6*p6*r**6
      assert margin>0 and r*r>n and r<outer
      cornerh4=n-F(3,n+2)*n*n;cornerh6=n-F(15,n+8)*n*n+F(30,(n+8)*(n+4))*n**3
      radialupper=F(577351,10**6) if n==3 else F(1,4)
      if n==3:assert radialupper>0 and 3*radialupper**2>1
      cornerupper=radialupper-A-b4*cornerh4-b6*cornerh6
      coremass=(4*F(157,50) if n==3 else 2*F(157,50)**2)*(r*r/2-A*r**n/n)
      lower=coremass/(2**n*cornerupper)
      assert cornerupper>0 and lower>2**n
      costs.append({'n':n,'core_margin_exact':str(margin),'core_mass_lower_exact':str(coremass),'corner_upper_exact':str(cornerupper),'specific_kernel_cost_lower_exact':str(lower),'display_only':float(lower),'scope':'only the two fixed kernels, not C_n or a whole-family obstruction'})
    return {'status':'PASS','source_sha256':sha,'u06_source_sha256':costsha,'spline_checks':rows,'selected_kernel_failures':costs,'shared_libraries':{'Fraction':'Python standard library','SymPy':sp.__version__},'scope':'independent finite exact identities and all support cells for one example per dimension; general correctness separately proved; no LP/search/author imports'}
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--source',type=Path,default=SOURCE);ap.add_argument('--cost-source',type=Path,default=COST);ap.add_argument('--save',action='store_true');args=ap.parse_args()
    out=run(args.source,args.cost_source);out['checker_sha256']=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    if args.save:HERE.joinpath('017_results.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out,indent=2))
