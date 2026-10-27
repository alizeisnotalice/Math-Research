"""R048 bounded symbolic/rational transcription checks. No production or old-baseline replay."""
from pathlib import Path
from fractions import Fraction as F
import hashlib,json
import sympy as sp
R,g,T=sp.symbols('R g T',positive=True)
Psi=lambda s:-2*R*sp.log(1-s/(2*R))-s
u=g/(2*R)
assert sp.simplify(sp.diff(Psi(g),g)-g/(2*R-g))==0
assert sp.simplify(sp.diff(Psi(g),R)-(-2*sp.log(1-u)-2*u/(1-u)))==0
primitive=g*g/2-(2*R-g)*Psi(g)
assert sp.simplify(sp.diff(primitive,g)-Psi(g))==0
assert sp.simplify(primitive.subs(g,0))==0
A1=2*R*sp.log(1+T/(2*R-g))
stableA=T+Psi(g)-Psi(g-T)
assert sp.simplify(sp.diff(A1-stableA,T))==0
assert sp.simplify((A1-stableA).subs(T,0))==0
checks={
 'Psi_s_exact':True,'Psi_R_exact':True,'baseline_primitive_derivative_and_zero':True,
 'exact_J_piecewise_equivalence_by_derivative_and_zero':True,
 'bridge_1280epsilon':1280*F(1,2**20)==F(5,4096),
}
d=F(1,2**24);sigma=32*d*d;Delta=d/65536
checks.update({'endpoint_sigma':sigma==F(1,2**43),'endpoint_actual_budget':Delta==F(1,2**40),'endpoint_buffer_equals_quarter_delta':64*d*d==Delta/4,'cap_exponent_positive':12*d<F(1,64),'half_delta_display':Delta/2==F(1,2**41)})
L0=2**240
E=F(2,L0**2)+F(2,1)/(sigma*sigma*L0*L0)+F(16,1)/(sigma*sigma*2**360)
checks['full_endpoint_scalar_margin']=E<Delta/2-2*sigma
checks['R042_full_cost_lower_gap']=F(4,L0)+F(1,2**111)<F(1,2**41)
checks['d30_full_cost_lower_gap']=F(4,L0)+F(1,2**46)<F(1,2**41)
assert all(checks.values())
root=Path(__file__).resolve().parent
receipt={'status':'PASS','sympy_version':sp.__version__,'checks':checks,'scope':'Only R045 local formulas and R048 endpoint/scalar comparison transcription; symbolic derivative identities are interpreted on the positive domains proved in the manuscript. L0 is not an operator threshold. No n=7, old-baseline replay, production calculation or optimization.','script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
(root/'transcription_receipt.json').write_text(json.dumps(receipt,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(receipt,ensure_ascii=False,indent=2))
