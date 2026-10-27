"""Small symbolic/exact transcription checks, not an operator proof or threshold search."""
from fractions import Fraction as F
from pathlib import Path
import hashlib
import json
import sympy as sp
r,rho,g,y,q,d,h,zeta=sp.symbols('r rho g y q d h zeta', positive=True)
exact_integral=(rho-g+sp.integrate((1-y/(2*r))**2,(y,0,g)))/rho
claimed=1-g**2/(2*r*rho)+g**3/(12*r*r*rho)
assert sp.simplify(exact_integral-claimed)==0
rminus=q+d-h
pinned=(1+(q-zeta*rho)/rminus)/2
assert sp.simplify(pinned-(1-(d-h+zeta*rho)/(2*rminus)))==0
s,L=sp.symbols('s L',positive=True)
assert sp.simplify((s**3-s)/s**2-(s-1/s))==0
cE=F(1,2**10);H=2**20
checks={
 'free_coordinate_second_moment_exact':True,
 'pinned_coordinate_factor_exact':True,
 'new_scalar_xd2_equals_s_minus_d':True,
 'column_cutoff_integer_comparison':H>=4*144/cE,
 'row_upper_below_cutoff':H>=36,
 'deficit_constant_exact':cE/(2*H)==F(1,2**31),
 'raw_energy_constant':F(1,64*9)>cE,
 'moment_loss_coefficient_exact':(F(1,50)+F(3,4))/2==F(77,200),
 'row_log_error_limit_strict':F(1,4)+F(13,200)<2,
}
assert all(checks.values())
root=Path(__file__).resolve().parent
receipt={'status':'PASS','sympy_version':sp.__version__,'checks':checks,'scope':'Exact one-coordinate integration, pinned-face algebra and rational constant transcription only. No n=7 calculation, optimization, full-operator re-certification or explicit-threshold claim.','source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
(root/'second_moment_transcription_receipt.json').write_text(json.dumps(receipt,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(receipt,ensure_ascii=False,indent=2))
