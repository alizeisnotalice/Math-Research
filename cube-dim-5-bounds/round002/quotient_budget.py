"""Finite outward arithmetic budget, conditional on geometric volume audit."""
from flint import arb, ctx
from fractions import Fraction
import hashlib,json
from pathlib import Path
ctx.prec=256
old=Fraction(256,3); new=Fraction(14455,31104)
a=lambda q:arb(q.numerator)/q.denominator
ratio=old/new
factor=a(ratio)**(arb(1)/5)
coefficient=arb(5)/6*a(new)**(arb(1)/5)*arb.pi()**(-arb(9)/5)
rows=[{'epsilon':f'1/{d}','positive_face_short_time_norm_upper_ball':str(coefficient*(arb(1)/d)**(arb(6)/5))} for d in [32,64,128,256,512,1024]]
payload={'status':'conditional_budget_not_C5_bound','precision_bits':ctx.prec,'old_over_new_face_constant_exact':str(ratio),'norm_improvement_factor':str(factor),'equal_budget_epsilon_increase_factor':str(a(ratio)**(arb(1)/6)),'new_positive_face_norm_coefficient':str(coefficient),'rows':rows,'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
raw=json.dumps(payload,indent=2)+'\n'
Path(__file__).with_suffix('.json').write_text(raw)
print(raw)
