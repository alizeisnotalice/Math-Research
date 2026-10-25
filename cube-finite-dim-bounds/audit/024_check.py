"""Independent finite identities supporting U08 audit; no author imports."""
from pathlib import Path
from fractions import Fraction as F
import sympy as s
import hashlib,json,argparse
HERE=Path(__file__).resolve().parent
HASHES={'u08_cauchy_obstacle.md':'c1163c7b40ec5d5723228df0959fc00f7ae91197db49fe0bfc7f132fa6ad7aec','u08_kernel_transfer.md':'98ee459c6d24ad80e4bff39822d7f27746ff046710b0cc5ce3a64183d807f405'}
def run(root):
    for name,h in HASHES.items():assert hashlib.sha256((root/name).read_bytes()).hexdigest()==h
    M,Fm,pi=s.symbols('M F pi',positive=True);potential=[]
    for n in (3,4):
      z=(M/(pi**n*Fm))**s.Rational(1,n);direct=Fm*z+M/pi**n*z**(1-n)/(n-1);claimed=s.Rational(n,n-1)/pi*M**s.Rational(1,n)*Fm**s.Rational(n-1,n)
      assert s.simplify(direct-claimed)==0;potential.append({'n':n,'uniform_potential_integral':str(claimed)})
    x=s.symbols('x',real=True);h=1-s.exp(-x)-(1-s.exp(-1))*x
    assert h.subs(x,0)==h.subs(x,1)==0 and s.diff(h,x,2)==-s.exp(-x)
    eps=s.symbols('epsilon',positive=True);v=s.symbols('v',real=True)
    compact_cutoff_constant=s.integrate((1-v*v)**3,(v,-1,1));assert compact_cutoff_constant==s.Rational(32,35)
    # g_epsilon(x)=x^2(1-x^2/epsilon^2)_+^3, compact C2 diagnostic.
    cutoff_value=compact_cutoff_constant*eps/s.pi
    contact=s.integrate((1-v*v)**2/(2-v)**2,(v,-1,1));contact=s.simplify(contact)
    # On [-1/2,1/2] numerator>=9/16, denominator<=25/4: exact positive lower.
    assert F(9,16)/F(25,4)==F(9,100)
    t,a=s.symbols('t a',positive=True);w=s.exp(-t)/s.sqrt(t);mass=s.integrate(w,(t,0,s.oo));k=s.integrate(w*s.exp(-a*t),(t,0,s.oo));rho=s.simplify(a*k)
    assert mass==s.sqrt(s.pi) and s.simplify(k-s.sqrt(s.pi)/s.sqrt(1+a))==0
    assert s.simplify(rho-s.sqrt(s.pi)*a/s.sqrt(1+a))==0
    assert s.limit(w,t,0,dir='+')==s.oo and s.limit(rho,a,s.oo)==s.oo
    minus_derivative=s.simplify(-s.diff(w,t));assert s.simplify(minus_derivative-s.exp(-t)*(t**s.Rational(-1,2)+s.Rational(1,2)*t**s.Rational(-3,2)))==0
    weighted_small=s.integrate(w,(t,0,1))-w.subs(t,1);assert s.simplify(weighted_small-(s.sqrt(s.pi)*s.erf(1)-s.exp(-1)))==0
    return {'status':'PASS independent finite U08 analytical diagnostics','source_sha256':HASHES,'potential_constants':potential,'spectral_tail_chord_check':'h(0)=h(1)=0 and h_second=-exp(-x)<0; hence 1-exp(-x)>=(1-exp(-1))*x for 0<=x<=1','origin_cutoff_compact_C2_example_L_at_zero':str(cutoff_value),'contact_example':{'u':'product_j (1-x_j^2)_+^2, belongs H1','point':'(2,0,...,0)','L_u_exact':str(contact/s.pi),'strict_positive_lower':'9/(100*pi)','meaning':'nonlocal contact Lu is generally positive, not zero; no new kernel cost'},'fixed_infinite_Levy_source_diagnostic':{'kernel':'K=integral_0^infty exp(-t)/sqrt(t) p_t dt','mass_exact':str(mass),'k_of_psi':str(k),'rho_of_psi':str(rho),'nu_density_mixture_weight':str(minus_derivative),'nu_total_mass':'infinite','Levy_small_time_weighted_integral':str(weighted_small),'meaning':'one fixed analytic example tests that L1 K can require infinite Levy source; no cube majorization or cost optimization'},'libraries':{'SymPy':s.__version__,'Fraction':'Python standard library'},'scope':'finite symbolic/integral auxiliaries only; whole-space limit/CND/LK closure assessed in written proof, not inferred from diagnostics'}
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--upper-dir',type=Path,default=HERE.parent/'upper');ap.add_argument('--save',action='store_true');args=ap.parse_args();out=run(args.upper_dir);out['checker_sha256']=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    if args.save:HERE.joinpath('024_results.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out,indent=2))
