"""R070 new second-profile identities and explicit finite local guards.

Low-resource single process; no old baseline or production-grid replay.
"""
from pathlib import Path
from fractions import Fraction as F
import hashlib,json
import sympy as sp

HERE=Path(__file__).resolve().parent
R,u,rho,a,q,h,eps,t,w=sp.symbols('R u rho a q h eps t w',positive=True)
def Q(R,u):
    return -sp.log(1-u/(2*R))
raw=2*R*Q(R,u)**3/3
assert sp.simplify(sp.diff(raw,u)-Q(R,u)**2/(1-u/(2*R)))==0
K2=2-a*((-sp.log(a))**2+2*(-sp.log(a))+2)
assert sp.simplify(sp.diff(K2,a)+sp.log(a)**2)==0
assert K2.subs(a,1)==0

s=sp.Rational(3,2);Rq=s*q-h;gq=(s-1+eps)*q-h;aq=(1-eps)*q;rh=t*q
ellpin=1-gq/(2*Rq)
central=2*Rq*Q(Rq,gq)**3/(3*rh)
edge=Rq*(Q(Rq,gq)**3+Q(Rq,Rq-w)**3)/(3*rh)
outside=Rq*Q(Rq,gq)**3/(3*rh)
assert sp.simplify((edge-central).subs(w,aq))==0
assert sp.simplify((sp.diff(edge,q)-sp.diff(central,q)).subs(w,aq)-(1-eps)*Q(Rq,gq)**2/(2*rh*ellpin))==0
assert sp.simplify((edge-outside).subs(w,Rq))==0
assert sp.simplify((sp.diff(edge,q)-sp.diff(outside,q)).subs(w,Rq))==0
assert sp.simplify(sp.diff(Rq,q)-Rq/q-h/q)==0

hmax=F(1,2**20)
assert F(3,25)*(F(9,10)+F(6,5))==F(63,250)
derivative_coefficient=(21+2*hmax)/250
assert derivative_coefficient<F(1,10)<1
assert F(1,10)*F(3,2)*F(4,3)==F(1,5)
assert hmax/5<F(1,1000000)
assert hmax<F(1,100000) # first-profile shift <=x h still fits its budget.

# New mu2: accurate lower/upper log(6/5) from four positive atanh terms.
z=F(1,11);K=4
log_lo=2*sum((z**(2*k+1)/F(2*k+1) for k in range(K)),F(0))
log_tail=2*z**(2*K+1)/(F(2*K+1)*(1-z*z))
log_hi=log_lo+log_tail
mu_lo=1-5*log_hi-F(5,2)*log_hi**2
mu_hi=1-5*log_lo-F(5,2)*log_lo**2
assert F(5289,1000000)<mu_lo<mu_hi<F(529,100000)

n0=2000;qmin=F(2,3)
target_mean_lo=(1-F(1,n0))*(mu_lo-hmax*(mu_hi+F(1,25)))
target_mean_hi=mu_hi+F(1,250*n0)/qmin
assert target_mean_lo>F(528,100000)
assert target_mean_hi<F(1059,200000)
assert F(528,100000)-F(527,100000)>=F(1,200000)
assert F(53,10000)-F(1059,200000)==F(1,200000)
bad_coefficient=F(200000**2,125**2)/qmin**2
assert bad_coefficient==5760000

# Exact conditional-face variance bounds, removing arbitrary source pin.
eta=F(1,1024);C1=F(3,50);C2=F(1,125)
match_lo=F(5269,1000000);match_hi=F(5301,1000000)
xnmax=F(1,1000)
var_lo=(match_lo-C2/(n0*qmin))/(1+eta)-C1*C1*xnmax/(qmin*qmin*(1-eta)**2)
var_hi=match_hi/(1-eta)
assert var_lo>F(13,2500)
assert var_hi<F(531,100000)

# Finite local guard implications for both widths at L>=4096.
L0=4096
assert F(1,256*L0)==hmax
assert F(1,512*L0)<=hmax
assert F(L0**2,2)>=n0
assert F(2,L0)<xnmax
assert F(L0,3)>6
assert L0**4>=6144
assert F(3,L0**2)<F(1,2**20) # stronger J point normalization if wanted.
profile_bad_at_L0=bad_coefficient*24/F(L0**4)
assert profile_bad_at_L0<F(1,2**20)

result={'status':'PASS new local second-profile identities and rational bounds',
        'sympy':sp.__version__,
        'symbolic_checks':['phi^2/ell tail primitive','unweighted source phi^2 / baseline target primitive','new onset second-profile derivative jump','C1 at |w|=R','q-prefactor chain'],
        'derivative_coefficient_at_hmax':str(derivative_coefficient),
        'safe_derivative':'<1/(10rho), hence <=1/rho',
        'matching_shift':'<=x h/5 < x/1000000 when h<=2^-20',
        'log6over5_lower':str(log_lo),'log6over5_upper':str(log_hi),
        'mu2_lower':str(mu_lo),'mu2_upper':str(mu_hi),
        'mu2_display_interval':['0.005289','0.005290'],
        'baseline_target_mean_coefficient_lower':str(target_mean_lo),
        'baseline_target_mean_coefficient_upper':str(target_mean_hi),
        'baseline_second_profile_bad_probability':'<=5760000/n for n>=2000,h<=2^-20',
        'conditional_source_variance_coefficient_lower':str(var_lo),
        'conditional_source_variance_coefficient_upper':str(var_hi),
        'conditional_variance_target':['13/2500','531/100000'],
        'L0_local_guards':{'L0':L0,'candidate_hc':['1/256','1/512'],'profile_bad_upper':str(profile_bad_at_L0),'uses':'e^L>=L^2/2 and e^L>=L^4/24; x in [L/2,L]'},
        'verifier_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'scope':'New local H2/profile and full inverse captured-face variance only; remove source pin and subtract coordinate mean squares. Not mixture variance, actual-mu randomness, Cantelli/BE or a Cn theorem. No old pass/production/MATLAB/Git work.'}
(HERE/'verification.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
