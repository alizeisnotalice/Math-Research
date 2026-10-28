"""R070 independent H2 algebra and rational local-rate checks.

No author helper or prior baseline program is imported. Probability-space
and source/target-pin distinctions are proved in REVIEW.md.
"""
from fractions import Fraction as F
from pathlib import Path
import hashlib,json,os,sys,platform
from datetime import datetime,timezone
import sympy as sp

HERE=Path(__file__).resolve().parent
checks={}
def ok(name,p):
 assert bool(p),name
 checks[name]=True
def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
R,u,rho,q,h,eps,w=sp.symbols('R u rho q h eps w',positive=True)
Q=-sp.log(1-u/(2*R))
T2=2*R*Q**3/3
ok('tail_phi2_over_ell_once',sp.simplify(sp.diff(T2,u)-Q**2/(1-u/(2*R)))==0)
a=sp.symbols('a',positive=True)
K=2-a*((-sp.log(a))**2+2*(-sp.log(a))+2)
ok('source_phi2_full_primitive',sp.simplify(sp.diff(K,a)+sp.log(a)**2)==0)
ok('primitive_endpoint',K.subs(a,1)==0)
Rq=sp.Rational(3,2)*q-h
gq=(sp.Rational(1,2)+eps)*q-h
ok('new_relative_gap',sp.simplify(gq-Rq/3-(eps*q-2*h/3))==0)
ok('new_R_derivative',sp.diff(Rq,q)==sp.Rational(3,2))
ok('new_g_derivative',sp.diff(gq,q)==sp.Rational(1,2)+eps)
t=sp.symbols('t',positive=True)
ok('H2_prefactor_derivative',sp.simplify(sp.diff(Rq/(3*t*q),q)-h/(3*t*q*q))==0)
ok('Q_u_identity',sp.simplify(sp.diff(Q,u)-1/(2*R-u))==0)
ok('Q_R_identity',sp.simplify(sp.diff(Q,R)+u/(R*(2*R-u)))==0)
# The capture-onset jump is derivative of the residual Q^3 term through
# v'=g' versus v'=R', with R/rho prefactor fixed continuously.
Qpin,ellpin=sp.symbols('Qpin ellpin',positive=True)
jump=R/(3*rho)*3*Qpin**2/(2*R*ellpin)*(1-eps)
ok('source_onset_derivative_jump',sp.simplify(jump-(1-eps)*Qpin**2/(2*rho*ellpin))==0)
ok('zero_residual_C1',sp.limit(Q**2/(1-u/(2*R)),u,0,dir='+')==0)

hc=F(1,2**20); s=F(3,2); qmin=F(2,3)
# Independent log enclosure: alternating log1p even12/odd13 sums.
def log_bounds(a):
 z=1/a-1
 assert 0<z<1
 lo=sum(((-1)**(k+1)*z**k/F(k) for k in range(1,13)),F(0))
 hi=lo+z**13/F(13)
 return lo,hi
def K_bounds(a):
 lo,hi=log_bounds(a)
 return 2-a*(hi**2+2*hi+2),2-a*(lo**2+2*lo+2)
amin=F(5,6)
cmin=s-hc/qmin
amax=(1+1/cmin)/2
Kmin,Kmin_hi=K_bounds(amax)
Kmax_lo,Kmax=K_bounds(amin)
free_lo,free_hi=2*cmin*Kmin,2*s*Kmax
ok('log_interval_domain',0<amin<amax<1 and cmin>1)
ok('true_chain_free_mean_lower',free_lo>F(528,100000))
ok('true_chain_free_mean_upper',free_hi<F(5295,1000000))
ok('relative_domain_alone_insufficient',amin<F(999,1000)<1)
# A separate finite n guard pays the omitted free coordinate, while the
# target pin is nonnegative for the lower and bounded for the upper.
n0=10**6
full_lo=(1-F(1,n0))*free_lo
full_hi=free_hi+F(1,250*qmin*n0)
ok('finite_target_mean_lower_n1e6',full_lo>F(528,100000))
ok('finite_target_mean_upper_n1e6',full_hi<F(5295,1000000))
ok('good_target_lower_gap',full_lo-F(527,100000)>F(1,100000))
ok('good_target_upper_gap',F(53,10000)-full_hi>F(1,200000))

ok('strong_derivative_coefficient',F(21,250)+F(2,250)*hc<F(1,10))
ok('derivative_cubic_chain',F(1,3)*3*F(1,5)**2*(F(9,10)+F(6,5))==F(21,250))
ok('prefactor_derivative_payment',2*F(1,5)**3/(3*qmin)==F(2,250))
ok('bridge_from_width',F(1,10)/qmin*(2/s)==F(1,5))
ok('source_anypin_second_profile',F(3,2)/3*2*F(1,5)**3==F(1,125))
ok('target_truepin_half',F(3,2)/3*F(1,5)**3==F(1,250))

source_mean,source_second,J=sp.symbols('H1 H2 J',positive=True)
mass=sp.symbols('mass',positive=True)
# Under normalized inverse capture, raw first and second moments are H1/J,H2/J.
variance=source_second/J-(source_mean/J)**2
ok('centered_variance_identity',sp.simplify(variance*J**2-(source_second*J-source_mean**2))==0)
ok('matched_lower_transport',F(527,100000)-hc/5>F(52,10000))
ok('matched_upper_transport',F(53,10000)+hc/5<F(531,100000))
ok('per_target_failure_safe_constant',F(1,125)**2/(qmin**2*F(1,200000)**2)==5760000)

out={'status':'PASS','scope':'independent LOCAL H2 algebra/rates; finite n>=10^6 used only for these independent two-sided target-mean bounds; no operator or cap-conditioned law claim',
 'utc':datetime.now(timezone.utc).isoformat(),'pid':os.getpid(),'python':sys.executable,'python_version':platform.python_version(),'sympy':sp.__version__,'code_sha256':sha(Path(__file__)),
 'checks':checks,'check_count':len(checks),'recipe':'SymPy exact primitives/chain/jump/centralization; Fraction even12/odd13 log1p enclosures for new H2 rate only; hmax=2^-20',
 'values':{k:{'rational':str(v),'display_only':float(v)} for k,v in {'free_mean_lower':free_lo,'free_mean_upper':free_hi,'full_target_mean_lower_n1e6':full_lo,'full_target_mean_upper_n1e6':full_hi,'max_derivative_coefficient':(F(21,250)+F(2,250)*hc)}.items()},
 'no_author_or_old_baseline_program':True}
(HERE/'verification.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps({'status':'PASS','checks':len(checks),'values':{k:v['display_only'] for k,v in out['values'].items()}},indent=2))
