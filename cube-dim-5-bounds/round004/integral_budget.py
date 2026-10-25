from flint import arb,ctx,fmpq
from fractions import Fraction
from pathlib import Path
import json,hashlib
ctx.prec=256
A=lambda q:arb(fmpq(q.numerator,q.denominator))
q=Fraction(5,4)
# C_model=(5*||G||_q/g)^q. Minkowski: ||G||_q <= ||H_eps||_q+E_eps.
rows=[]
for g in [Fraction(105,10**6),Fraction(11,10**5),Fraction(115,10**6)]:
 for den in [256,512,1024]:
  eps=Fraction(1,den)
  E=arb(5)/6*(arb(14455)/31104)**(arb(1)/5)*arb.pi()**(-arb(9)/5)*A(eps)**(arb(6)/5)
  norm_limit=A(g)/5*arb(32)**(arb(4)/5)
  h_limit=norm_limit-E
  rows.append({'g':str(g),'eps':str(eps),'shorttime_norm_budget':str(E),'G_norm_strict_limit':str(norm_limit),'H_integral_strict_limit':str(h_limit**(arb(5)/4)) if h_limit>0 else None,'positive_remaining_norm':bool(h_limit>0)})
out={'status':'conditional_certification_budget_not_a_bound','interpretation':'A rigorous integral of H_eps^(5/4) strictly below the listed threshold, together with valid global g, suffices for model cost<32. No g or integral certified here.','rows':rows,'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
Path(__file__).with_suffix('.json').write_text(json.dumps(out,indent=2)+'\n')
for r in rows:print(r['g'],r['eps'],r['H_integral_strict_limit'])
