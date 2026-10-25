"""Round9 independent formulas, necessary-condition counterexamples and 2 anchors."""
from pathlib import Path
from fractions import Fraction as F
from itertools import product
import sympy as s
from flint import arb,ctx
import hashlib,json,argparse,time
HERE=Path(__file__).resolve().parent

def A(expr):
    if expr==s.pi:return arb.pi()
    if expr.is_Rational:return arb(int(expr.p))/int(expr.q)
    if expr.func==s.log:return A(expr.args[0]).log()
    if expr.is_Add:return sum((A(z) for z in expr.args),arb(0))
    if expr.is_Mul:
      out=arb(1)
      for z in expr.args:out*=A(z)
      return out
    if expr.is_Pow and expr.exp.is_Integer:return A(expr.base)**int(expr.exp)
    raise ValueError(expr)

def int01(expr,var,lo,hi):return s.integrate(s.cancel(expr),(var,lo,hi))
def toy_jump(alpha,a=s.Rational(99,100)):
    z=s.symbols('z',positive=True);c=s.Rational(1,2**alpha);d=(1-a)/2
    near=((1+z)**(-alpha)+(1-z)**(-alpha)-2)/z**2
    plus=int01(near,z,0,s.Rational(1,2))+int01(((1+z)**(-alpha)+z**(-alpha)-2)/z**2,z,s.Rational(1,2),1)+int01((z**(-alpha)-2+c)/z**2,z,1,2)-(1-c)
    minus=int01(near,z,0,d)+int01(((1+z)**(-alpha)+(a+z)**(-alpha)-2)/z**2,z,d,1)+int01(((a+z)**(-alpha)-2+c)/z**2,z,1,2-a)-2*(1-c)/(2-a)
    both=s.expand_log(s.expand(plus+minus),force=True);value=both/(4*s.pi*(1-c));bound=A(value);assert bound<0
    return {'n':alpha+1,'toy':'G=norm_inf(x)^(-(n-1)), c=2^(-(n-1)), g=1','point':'(1,99/100,0,...)','geometry_only_cost':str(s.Rational(alpha)*2/(1-c)),'L_K_exact':str(value),'L_K_rigorous_ball':bound.str(30),'fails_source_positive_away_origin':True}

def residue_integral(R,z,poles):
    residues=[]
    for pole,order in poles:
      smooth=s.cancel((z-pole)**order*R)
      if order==2:smooth=s.diff(smooth,z)
      residues.append(s.cancel(smooth.subs(z,pole)))
    return s.factor(s.cancel(2*s.pi*s.I*sum(residues)))

def run():
    ctx.prec=128;r,c,g,h,J,E=s.symbols('r c g h J E',positive=True);rows=[]
    for n,gamma in [(3,s.Rational(3,2)),(4,s.Rational(5,2))]:
      alpha=n-1;p=s.Rational(n,alpha);R=(h/c)**s.Rational(1,alpha)
      direct=h*R-c*R**n/n;claimed=s.Rational(alpha,n)*h**p*c**(-s.Rational(1,alpha));assert s.simplify(direct-claimed)==0
      mass=s.Rational(alpha,n)*J*c**(-s.Rational(1,alpha));cost=mass/(2**n*(g-c))
      assert s.simplify(s.diff(s.log(cost),c)-(n*c-g)/(alpha*c*(g-c)))==0
      optimum=s.simplify(cost.subs(c,g/n).subs(J,n*2**n*g**p*E));assert s.simplify(optimum-n**p*E)==0
      q=s.symbols('q',positive=True);angular1=s.integrate(q**(alpha-1-gamma),(q,0,1));angularp=s.integrate(q**(alpha-1-gamma*p),(q,0,1));assert angular1==2 and angularp==s.oo
      rows.append({'n':n,'raw_mass':'(n-1)/n * J * c^(-1/(n-1))','optimal_c':'g/n','optimal_cost':str(optimum),'counterexample_gamma':str(gamma),'counterexample_angular_G_integral':str(angular1),'counterexample_angular_G_power_integral':'infinity'})
    jumps=[toy_jump(2),toy_jump(3)]
    t=s.symbols('t',positive=True);z=s.symbols('z');a=s.Rational(1,2)
    B=s.Matrix([[1,1,1],[1,1,-1],[1,-1,1]]);v4=s.Matrix([1,-1,-1]);assert abs(B.det())==4 and B.inv()*v4==s.Matrix([-1,1,1])
    assert B.inv()*s.Matrix([1,0,0])==s.Matrix([0,a,a]);assert B.inv()*s.Matrix([1,1,0])==s.Matrix([a,a,0])
    axisR=1/((z*z+t*t)**2*((z-a)**2+t*t)**2);edgeR=1/((z*z+t*t)**2*((z-a)**2+t*t)*((z+a)**2+t*t))
    axis=s.factor(t**4/(4*s.pi**4)*residue_integral(axisR,z,[(s.I*t,2),(a+s.I*t,2)]));expected=t*(a*a+20*t*t)/(4*s.pi**3*(a*a+4*t*t)**3);assert s.simplify(axis-expected)==0
    edge=s.factor(t**4/(4*s.pi**4)*residue_integral(edgeR,z,[(s.I*t,2),(a+s.I*t,1),(-a+s.I*t,1)]))
    Gaxis=s.integrate(axis,(t,0,s.oo));Gedge=s.integrate(edge,(t,0,s.oo));assert s.simplify(Gaxis-3/(8*s.pi**3))==0 and s.simplify(Gedge-(4*s.log(2)-1)/(8*s.pi**3))==0
    return {'status':'PASS independent fixed-formula and adversarial audit','mass_and_optimization':rows,'positive_homogeneous_but_illegal_source_toys':jumps,'n3_independent_residue_anchors':{'p_axis':str(axis),'p_edge':str(edge),'G_axis':str(Gaxis),'G_edge':str(Gedge),'scope':'only two fixed points, hence upper bounds on true g; no global face minimum certified'},'analytic_conditions':'For actual transient probability semigroup with pt infinity<=C t^-n: split short-time L1/long-time Linfinity yields weak L^{n/(n-1)}, hence J=n*volume{G>1}<infinity. Semigroup Jensen gives L(G-c)+ >= -delta0 distributionally; complete source theory separately locked/audited in026.','shared_libraries':{'SymPy':s.__version__,'arb':'python-flint, 128 bits','Fraction':'Python standard library'},'scope':'no cone code, no theory implementation, no g/J grid, no new cost or global bound'}
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--save',action='store_true');args=ap.parse_args();start=time.monotonic();out=run();out['elapsed_seconds']=time.monotonic()-start;out['checker_sha256']=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    if args.save:HERE.joinpath('025_results.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out,indent=2))
