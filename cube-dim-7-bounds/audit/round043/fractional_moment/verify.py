"""New R043 symbolic/Arb checks, not the author's ten identity checks.

Usage: python verify.py FINAL_PROOF.md [ORIGINAL_PROOF.md]
Runtime: existing local SymPy/python-flint environment; no MATLAB/production.
"""
from fractions import Fraction as F
from pathlib import Path
import hashlib,json,re,sys
import sympy as sp
from flint import arb,fmpq,ctx

HERE=Path(__file__).resolve().parent
a,k,p,s,A,d,q=sp.symbols('a k p s A d q',positive=True)
# k=2r/rho, a=1-g/(2r). Real-p antiderivative, checked by differentiation
# and the zero-gap boundary rather than integer-moment enumeration.
moment=1-k*(1-a)+k*(1-a**(p+1))/(p+1)
assert sp.simplify(sp.diff(moment,a)-k*(1-a**p))==0
assert sp.simplify(moment.subs(a,1)-1)==0
mean=1-k*(1-a)**2/2
assert sp.simplify(moment.subs(p,1)-mean)==0
variance=sp.expand(moment.subs(p,2)-mean**2)
assert sp.simplify(variance-(k*(1-a)**3/3-k*k*(1-a)**4/4))==0
raw_entropy=k*((a*a-1)/4-a*a*sp.log(a)/2)
assert sp.simplify(sp.diff(moment,p).subs(p,1)-raw_entropy)==0
poly=sp.series(moment.subs(a,1-s),s,0,5).removeO()
logratio=sp.series(sp.log(poly),s,0,5).removeO()-p*sp.series(sp.log(mean.subs(a,1-s)),s,0,5).removeO()
coeff=sp.expand(logratio)
assert all(sp.simplify(coeff.coeff(s,j))==0 for j in (0,1,2))
assert sp.simplify(coeff.coeff(s,3)-k*p*(p-1)/6)==0
assert sp.simplify(coeff.coeff(s,4)+k*p*(p-1)*(p-2)/24+k*k*p*(p-1)/8)==0
# Independent derivation of the tilted-variance and Chebyshev constants.
assert F(64,27)*F(1,12)==F(16,81)<1
assert F(1,12)-F(1,32)==F(5,96)>F(1,24)
assert sp.simplify(A/(A/64-A/128)**2-16384/A)==0
assert sp.simplify(sp.diff(d*sp.exp(-q*d**3),d)/sp.exp(-q*d**3)-(1-3*q*d**3))==0

# One different admissible exact geometry point, with Arb-directed logs.
ctx.prec=512
def ab(v):
    return arb(fmpq(v.numerator,v.denominator))
r,rho=F(5),F(13)
aa=F(31,32)**3
g=2*r*(1-aa)
assert 0<g<=min(r/2,rho/2)
m=1-g*g/(4*r*rho)
entropy=(ab(r*(aa*aa-1)/(2*rho))-ab(r*aa*aa/rho)*ab(aa).log())/ab(m)-ab(m).log()
lower=g**3/(64*r*r*rho)
assert entropy>ab(lower)
entropy_mid=str(entropy)

final=Path(sys.argv[1])
receipt={'status':'PASS independent symbolic real-p identities and new Arb point',
         'source_sha256':hashlib.sha256(final.read_bytes()).hexdigest(),
         'sympy':sp.__version__,'arb_precision_bits':ctx.prec,
         'symbolic_checks':['real-p moment antiderivative and boundary','mean/variance exact','entropy equals p derivative at one','cubic and quartic log-normalized coefficients','tilted variance coefficient 16/81','Chebyshev constant 16384','uniform-cap scalar derivative'],
         'new_arb_point':{'r':str(r),'rho':str(rho),'g':str(g),'entropy_enclosure':entropy_mid,'entropy_lower':str(lower),'strict_positive_margin':True},
         'verifier_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
         'scope':'Only the frozen independent-coordinate product/envelope law. No actual C_A energy necessity from C_A<=Y, no Cn barrier, no author ten-check replay, no production/n7/MATLAB/Git mutation.'}
if len(sys.argv)>2:
    old=Path(sys.argv[2])
    blocks=lambda t:re.findall(r'\\\[(.*?)\\\]',t,re.S)
    old_blocks,new_blocks=blocks(old.read_text()),blocks(final.read_text())
    assert len(old_blocks)>=15 and old_blocks==new_blocks
    receipt['original_source_sha256']=hashlib.sha256(old.read_bytes()).hexdigest()
    receipt['all_displayed_math_blocks_unchanged']=True
    receipt['unchanged_displayed_math_block_count']=len(old_blocks)
(HERE/'verification.json').write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps(receipt,indent=2))
