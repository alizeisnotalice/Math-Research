"""New R064 relative-coordinate chain and exact rate checks.

No replay of old baseline integrals or production numerical batches.
Usage: python verify.py [current PROOF.md]
"""
from pathlib import Path
from fractions import Fraction as Q
import hashlib,json,sys
import sympy as sp

HERE=Path(__file__).resolve().parent
q,h,eps,t,w=sp.symbols('q h eps t w',positive=True)
s=sp.Rational(3,2)
R=s*q-h;g=(s-1+eps)*q-h;a0=(1-eps)*q;rho=t*q
ellpin=(R+a0)/(2*R)
assert sp.diff(R,q)==s
assert sp.diff(g,q)==sp.Rational(1,2)+eps
assert sp.simplify(g-R/3-(eps*q-2*h/3))==0
assert sp.simplify(sp.diff(ellpin,q)+(1-eps)*h/(2*R*R))==0

def FF(r,u):
    return -2*r*sp.log(1-u/(2*r))-u
RR,U=sp.symbols('RR U',positive=True)
FR=sp.diff(FF(RR,U),RR);FU=sp.diff(FF(RR,U),U)
chain=s*FR.subs({RR:R,U:g})+(s-1+eps)*FU.subs({RR:R,U:g})
assert sp.simplify(sp.diff(FF(R,g),q)-chain)==0
central=1+FF(R,g)/rho
edge=1+(FF(R,g)+FF(R,R-w)-(w-a0))/(2*rho)
outside=1+(FF(R,g)-(w-a0))/(2*rho)
assert sp.simplify((edge-central).subs(w,a0))==0
assert sp.simplify((sp.diff(edge,q)-sp.diff(central,q)).subs(w,a0)-(1-eps)/(2*rho*ellpin))==0
assert sp.simplify((sp.diff(edge,q)-sp.diff(outside,q)).subs(w,R))==0

# New relative-domain constants, not old absolute-d endpoint checks.
umax=Q(1,6)
assert umax/(1-umax)==Q(1,5)
assert 2*umax**2/(1-umax)==Q(1,15)
assert Q(3,2)*Q(1,15)+Q(1,5)==Q(3,10)
assert Q(3,10)+Q(3,2)*(Q(1,15)+Q(1,5))+1==Q(17,10)
assert Q(3,2)*(Q(1,15)+Q(1,3))+Q(3,1000)==Q(603,1000)<1
assert 1-Q(1,24)-Q(1,4000)==Q(11497,12000)>Q(1,2)
assert Q(17,10)+Q(603,1000)<3
assert 256*Q(3,2)*Q(4,3)==512
assert 512*Q(1,2**64)==Q(1,2**55)

# Odd alternating log Taylor upper bound at log(6/5).
z=Q(1,5)
Q5=sum(((-1)**(k+1)*z**k/k for k in range(1,6)),Q(0))
e=Q(1,1024)
meanrate=Q(3,2)*Q5**2/(1-e)
varrate=Q5**3/(1-e)
assert meanrate<Q(999,20000)
assert varrate<Q(1,150)

record={'status':'PASS new multiplicative chain, kink and exact rate checks',
        'sympy':sp.__version__,
        'relative_guards':'s=3/2; active R>=1; q>=2/3; 1<=r<3/2; rho>=4r; 0<h<=1/1000; 0<=epsilon and epsilon*q<=h/s; g>0',
        'new_symbolic_checks':['Rprime=3/2 and gprime=1/2+epsilon','g<=R/3 equivalent to epsilon*q<=2h/3','ellpin derivative negative','full F chain with new derivatives','onset one-sided derivative jump','C1 at |w|=R'],
        'derivative_bound':'strictly <3/rho a.e.; hence 256/rho',
        'product_log_shift':'<=512 L h=2^-55 for h=2^-64/L',
        'Q5_log_upper':str(Q5),'mean_rate_upper':str(meanrate),'variance_rate_upper':str(varrate),
        'verified_mean_target':'999/20000','verified_variance_target':'1/150',
        'verifier_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'scope':'Local J/Psi/phi and conditional-face estimates only. Baseline target W is uniform [-R,R]; phi fullsource Z is uniform [-rho,rho]. No global Cn, target count, radial column, old/production/MATLAB/n7/Git work.'}
if len(sys.argv)>1:
    source=Path(sys.argv[1]);record['source_sha256']=hashlib.sha256(source.read_bytes()).hexdigest()
(HERE/'verification.json').write_text(json.dumps(record,indent=2)+'\n')
print(json.dumps(record,indent=2))
