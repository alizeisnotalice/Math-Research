"""R056 new symbolic phi moments and enlarged-domain constant checks.

No historical baseline replay, production quadrature, or MATLAB.
Usage: python verify.py [current PROOF.md]
"""
from fractions import Fraction as F
from pathlib import Path
import hashlib,json,sys
import sympy as sp

HERE=Path(__file__).resolve().parent
R,rho,a=sp.symbols('R rho a',positive=True)
I1=R/rho*sp.log(a)**2
I2=-2*R/(3*rho)*sp.log(a)**3
assert sp.simplify(sp.diff(I1,a)-2*R/rho*sp.log(a)/a)==0
assert sp.simplify(sp.diff(I2,a)+2*R/rho*sp.log(a)**2/a)==0
assert I1.subs(a,1)==I2.subs(a,1)==0

d=F(1,3)
# R>=1, g<=d imply u=g/(2R)<=1/6.
assert 1/(2-d)==F(3,5)<1
assert 1/(2*(2-d))==F(3,10)
assert 1+d+3*d*d==F(5,3)<2
assert F(6,5)*d*d+F(5,2)*d==F(29,30)<1
assert 1-F(5,32)*d==F(91,96)>F(1,2)
assert F(2)+F(3,5)<3
assert F(2,9)*F(1,2)+F(8,27)==F(11,27)<3
assert 1280*F(1,2**24)==F(5,65536)
assert F(3,10)>=F(1,4) # K<=n<=4nd for the candidate d, not every d<=1/3.

# Exact identity for the mixture correction, distinct from conditional variance.
m0,m1,m2,p0,p1,c=sp.symbols('m0 m1 m2 p0 p1 c',real=True)
weights=[p0,p1,1-p0-p1]
mus=[m0,m1,m2]
means=[c+sum(mus)-m for m in mus]
avg=lambda values:sum(p*v for p,v in zip(weights,values))
mixvar=sp.expand(avg([v*v for v in means])-avg(means)**2)
excludedvar=sp.expand(avg([v*v for v in mus])-avg(mus)**2)
assert sp.simplify(mixvar-excludedvar)==0
positive_correction=(F(0)**2+F(1,7)**2)/2-((F(0)+F(1,7))/2)**2
assert positive_correction==F(1,196)>0

record={'status':'PASS new R056 local identities and enlarged-domain constants',
        'sympy':sp.__version__,
        'new_checks':['exact full-source E(phi/ell) primitive','exact full-source E(phi^2/ell) primitive','d<=1/3 derivative constants at monotone endpoint','baseline second-moment coefficient 11/27<3','full matching shift 5/65536','face-mixture mean-variance correction identity'],
        'nonzero_mixture_mean_variance_example':'1/196',
        'verifier_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'scope':'Local J/Psi/phi only; target W uniform [-R,R] for baseline, source Z uniform [-rho,rho] for phi moments. No full Cn theorem audit or old/production/MATLAB/Git work.'}
if len(sys.argv)>1:
    source=Path(sys.argv[1])
    record['source_sha256']=hashlib.sha256(source.read_bytes()).hexdigest()
(HERE/'verification.json').write_text(json.dumps(record,indent=2)+'\n')
print(json.dumps(record,indent=2))
