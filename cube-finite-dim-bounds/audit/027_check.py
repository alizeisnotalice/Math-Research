"""Independent fixed single-bump source obstruction, no author/Green imports."""
from pathlib import Path
from fractions import Fraction as F
import json,hashlib,argparse,time
import sympy as s
from flint import arb,acb,ctx
HERE=Path(__file__).resolve().parent

def ends(x):return [str(x.lower().fmpq()),str(x.upper().fmpq())]
def subset(inner,outer):
    il,iu=map(F,ends(inner));ol,ou=map(F,ends(outer));return ol<il<=iu<ou

def simplex_integrate(f,a,b,d):
    return s.integrate(s.integrate(s.integrate(f,(d,0,1-a-b)),(b,0,1-a)),(a,0,1))
def run(root):
    started=time.monotonic();snapshot=json.loads((HERE/'027_source_snapshot.json').read_text())
    for name,h in snapshot.items():assert hashlib.sha256((root/name).read_bytes()).hexdigest()==h,(name,'source changed')
    # Rebuild the Feynman/Jensen constant using scalar integrals and Dirichlet moments.
    aa,bb,dd=s.symbols('alpha0 alpha1 alpha2',real=True);volume=simplex_integrate(s.Integer(1),aa,bb,dd);moment=simplex_integrate(aa*bb,aa,bb,dd)
    assert volume==s.Rational(1,6) and moment==s.Rational(1,120) and moment/volume==s.Rational(1,20)
    z=s.symbols('z',real=True);A,V=s.symbols('A V',positive=True);t=s.symbols('t',positive=True)
    scalar_line=s.integrate((z*z+A)**(-4),(z,-s.oo,s.oo));scalar_time=s.integrate(t**4/(t*t+V)**s.Rational(7,2),(t,0,s.oo))
    assert scalar_line==5*s.pi/(16*A**s.Rational(7,2)) and scalar_time==1/(5*V)
    feynman_coefficient=s.cancel(s.Rational(1,4)/s.pi**4*6*(5*s.pi/16)*s.Rational(1,5));assert feynman_coefficient==3/(32*s.pi**3)
    x,y,w=s.symbols('x y w',real=True);centers=[s.Integer(0),-(y+w)/2,(x-w)/2,(x-y)/2];pair_sum=s.expand(sum((centers[i]-centers[j])**2 for i in range(4) for j in range(i+1,4)));assert pair_sum==x*x+y*y+w*w
    lower=s.cancel(feynman_coefficient*volume*20);assert lower==5/(16*s.pi**3)
    # Independent local b calculation: no PV derivative formula or general Lb code.
    slice_integral=s.integrate((2+3*z)**2,(z,s.Rational(-2,3),0));lb_exact=4*slice_integral/s.pi;assert slice_integral==s.Rational(8,9) and lb_exact==32/(9*s.pi)
    # Reciprocal substitution maps both infinite source tails to finite intervals.
    U,RR=s.symbols('u R',real=True)
    for sign in (-1,1):
      denominator=U*U+2*sign*U+3
      primitive=(1-RR)*U+RR*sign*s.log(denominator)+RR/s.sqrt(2)*s.atan((U+sign)/s.sqrt(2))
      assert s.cancel(s.diff(primitive,U)-(1-RR*U*U/denominator))==0
    ctx.prec=192;pi=arb.pi();two=arb(2);T=4*two.log()-1;c=T/(24*pi**3);kappa=arb(5)/(16*pi**3);R=kappa/c;q=(3*R-2).sqrt();a=(-1-q)/3;b=(-1+q)/3
    assert T>0 and R>3 and kappa/3>c and kappa>c
    pieces=[];closed=arb(0)
    for sign,limit in [(-1,-1/a),(1,1/b)]:
      def f(v,analytic,sgn=sign):return 1-acb(R)*v*v/(v*v+2*sgn*v+3)
      piece=acb.integral(f,0,acb(limit),rel_tol=arb('1e-44'),abs_tol=arb('1e-48'),eval_limit=10000,depth_limit=30,deg_limit=80)
      assert piece.is_finite();pieces.append(piece)
      closed+=(1-R)*limit+R*sign*((limit*limit+2*sign*limit+3)/3).log()+R/two.sqrt()*(((limit+sign)/two.sqrt()).atan()-(arb(sign)/two.sqrt()).atan())
    I=sum((v.real for v in pieces),arb(0));assert I>0 and subset(closed,I)
    stated=(R-1)*(1/a-1/b)+2*R*((-a)/b).log()+R/two.sqrt()*(pi-2*(q/two.sqrt()).atan());assert subset(stated,I)
    mu_expression=4*c*I/pi;lb=arb(32)/(9*pi);eps_bound=mu_expression/lb;candidate=arb(1)/200;negative_expression=mu_expression-candidate*lb
    assert F(ends(eps_bound)[1])<F(3347,10**6)
    assert F(ends(negative_expression)[1])<F(-1871,10**6)
    published=json.loads((root/'results/round009_bump_source_certificate.json').read_text())
    for name,value in [('mu_density_upper_ball',mu_expression),('epsilon_necessary_upper_ball',eps_bound),('corrected_source_density_upper_ball',negative_expression)]:
      outer=arb(published[name].strip('[]'));assert subset(value,outer)
    return {'status':'PASS independent fixed single-bump necessary source certificate','source_sha256':snapshot,'feynman':{'simplex_volume':str(volume),'simplex_pair_moment':str(moment),'normalized_pair_moment':str(moment/volume),'Green_simplex_coefficient':str(feynman_coefficient),'pair_distance_square_sum':str(pair_sum),'G_lower_coefficient':str(lower)},'b_local':{'slice':'b(e1+s*v)=s^2*(2+3s)^2 on [-2/3,0], zero elsewhere','slice_integral':str(slice_integral),'Lb_e1_exact':str(lb_exact),'regularity':'C_c^{1,1}; continuous full directional generator, no sphere atom'},'parameters':{'precision_bits':192,'fixed_c':'(4 log(2)-1)/(24 pi^3)','R_exact_expression':'15/(2*(4 log(2)-1))','e1_inside_Omega':True,'whole_cube_inside_Omega_from_lower':True,'g_equals_edge_assumed':False},'tail_integral':{'route':'s -> -1/u on negative tail, s -> 1/u on positive tail; two finite rational Acb integrals','finite_integral_balls':[v.str(55) for v in pieces],'I_interval':ends(I),'independent_antiderivative_interval':ends(closed),'published_closed_form_interval':ends(stated)},'mu_upper_expression_interval':ends(mu_expression),'mu_upper_expression_ball':mu_expression.str(55),'epsilon_necessary_upper_interval':ends(eps_bound),'epsilon_necessary_upper_ball':eps_bound.str(55),'simple_strict_necessary_condition':'epsilon < 3347/1000000','rejected_epsilon':'1/200','corrected_density_upper_expression_interval':ends(negative_expression),'corrected_density_upper_expression_ball':negative_expression.str(55),'strict_corrected_source_density_upper_below':'-1871/1000000','all_three_expression_intervals_strictly_inside_published_balls':True,'scope':'only fixed 3D c and unit-ball bump; necessary upper bound plus one rejected correction; no C_n, cost/g, joint-bump or full-cb9 certificate','shared_libraries':{'SymPy':s.__version__,'Arb/Acb':'python-flint, 192bit','Fraction':'Python standard library'},'implementation':'no author checker or cone module imports; scalar Feynman factors plus reciprocal finite-tail integration; published expression interval is an upper-bound expression, not enclosure of unknown actual mu','elapsed_seconds':time.monotonic()-started}
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--project',type=Path,default=HERE.parent);ap.add_argument('--save',action='store_true');args=ap.parse_args();out=run(args.project);out['checker_sha256']=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    if args.save:HERE.joinpath('027_results.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({k:v for k,v in out.items() if k!='source_sha256'},indent=2))
