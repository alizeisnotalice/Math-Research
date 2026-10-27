"""New R066 first-profile symbolic identities and exact rational budgets.

This does not replay old J baselines or a production grid.
"""
from pathlib import Path
from fractions import Fraction as F
import hashlib,json
import sympy as sp

HERE=Path(__file__).resolve().parent
R,u,rho,a,q,h,eps,t,w=sp.symbols('R u rho a q h eps t w',positive=True)
def Q(R,u):
    return -sp.log(1-u/(2*R))
def H(R,u):
    return R*Q(R,u)**2
assert sp.simplify(sp.diff(H(R,u),u)-Q(R,u)/(1-u/(2*R)))==0
mean=2*R/rho*(1-a+a*sp.log(a))
assert sp.simplify(sp.diff(mean,a)-2*R/rho*sp.log(a))==0
assert mean.subs(a,1)==0
assert sp.diff(1-a+a*sp.log(a),a)==sp.log(a)

s=sp.Rational(3,2);Rq=s*q-h;gq=(s-1+eps)*q-h;aq=(1-eps)*q;rh=t*q
ellpin=1-gq/(2*Rq)
central=H(Rq,gq)/rh
edge=(H(Rq,gq)+H(Rq,Rq-w))/(2*rh)
outside=H(Rq,gq)/(2*rh)
assert sp.simplify((edge-central).subs(w,aq))==0
assert sp.simplify((sp.diff(edge,q)-sp.diff(central,q)).subs(w,aq)-(1-eps)*Q(Rq,gq)/(2*rh*ellpin))==0
assert sp.simplify((edge-outside).subs(w,Rq))==0
assert sp.simplify((sp.diff(edge,q)-sp.diff(outside,q)).subs(w,Rq))==0
assert sp.simplify(sp.diff(Rq,q)-Rq/q-h/q)==0

# New first-profile derivative constants, including the q-dependent prefactor.
assert F(3,2)*F(1,5)+F(3,5)==F(9,10)
assert F(3,2)*(F(1,5)+F(3,5))==F(6,5)
assert F(2,5)*(F(9,10)+F(6,5))==F(21,25)
assert (F(21)+3*F(1,1000))/50<F(1,2)<8
assert (2*F(27,50)+F(3,25))/2==F(3,5)<8
assert 8*F(3,2)*F(4,3)==16
assert 16*F(1,2**80)==F(1,2**76)<F(1,100000)

# To upper-bound 1/2-(5/2)log(6/5), use a LOWER log bound.
z=F(1,5)
log_lower=sum(((-1)**(k+1)*z**k/k for k in range(1,9)),F(0))
mu_upper=F(1,2)-F(5,2)*log_lower
assert 0<mu_upper<F(221,5000)
margin=F(221,5000)-mu_upper
target_event=F(221,5000)+F(1,100000)
transported=target_event+F(1,100000)
assert transported==F(2211,50000)
eJ=F(1,2**20)
source_rate=transported/(1-eJ)
assert source_rate<F(223,5000)
Q5=sum(((-1)**(k+1)*z**k/k for k in range(1,6)),F(0))
variance_rate=Q5**3/(1-eJ)
assert variance_rate<F(1,150)
assert F(9,200)-F(223,5000)-F(1,5000)==F(1,5000)
assert F(5000**2,150)==F(500000,3)

# Uniform target-profile Chebyshev bookkeeping: one deterministic pin.
C=F(3,50);Cpin=F(3,100);qmin=F(2,3)
assert Cpin/(9000*qmin)==F(1,200000)
bad_coefficient=C*C*200000**2/qmin**2
assert bad_coefficient==324000000

result={'status':'PASS new first-profile identities and exact budgets',
        'sympy':sp.__version__,
        'symbolic_checks':['tail primitive H=R Q^2','target Fubini-mean primitive','new Jphi onset derivative jump','C1 at |w|=R','R/rho prefactor derivative cancellation'],
        'derivative_bound':'<(21+3h)/(50 rho)<1/(2rho), hence <=8/rho',
        'profile_shift':'<=16 L h=2^-76; also <=16 x h using q>=2/3',
        'even_eight_log_lower':str(log_lower),
        'target_mean_upper':str(mu_upper),'strict_margin_to_221_over_5000':str(margin),
        'target_event_rate':str(target_event),'transported_rate':str(transported),
        'source_mean_rate_after_J':str(source_rate),'variance_rate_upper':str(variance_rate),
        'source_chebyshev_bound':'500000/(3x) once pin/x<=1/5000',
        'target_profile_bad_probability':'<=324000000/n for n>=9000 under stated local guards',
        'verifier_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'scope':'New Jphi and first-profile interfaces only. Target full-face law versus source inverse-face law are distinct; mixing uses conditional probability bounds. No old passed baseline/grid, full Cn, MATLAB, n7, lower-bound or Git work.'}
(HERE/'verification.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
