"""New theorem changes B; six finite acceptance budgets, not an I certificate."""
from pathlib import Path
from flint import arb,ctx
from fractions import Fraction
import json,hashlib,time
p=Path(__file__).resolve().parents[1]; started=time.monotonic();ctx.prec=256
source=p/'upper/round010_alpha_strata.md'
assert hashlib.sha256(source.read_bytes()).hexdigest()=='8cdf321d9692fbb8d71180671cb090698e937c3bca3cc3ba32c84ff1c91f99ba'
A=2*(arb(2)/3).gamma()/(3*arb.pi());K=(2/arb.pi()).sqrt()/arb(3)**(arb(1)/4)
P=(arb(14455)/191102976)*(arb(11)/3).gamma()*arb(15)**(arb(4)/3)/(2*arb.pi())**4
B=arb(5)/6*K*(arb(16)**(arb(7)/2)*P)**(arb(3)/10)
I0=arb(11)/1000000;q=arb(10)/7
rows=[]
for denom in [1024,4096,16384]:
 for T in [1,4]:
  e=arb(1)/denom;short=B*e**(arb(6)/5);long=arb(3)/7*A**5/48*arb(T)**(-arb(7)/3)
  rem=I0**(arb(7)/10)-short-long;assert rem>0
  gate=rem**q
  rows.append({'epsilon':f'1/{denom}','T':T,'short_norm_upper_ball':short.str(65,more=True),'long_norm_upper_ball':long.str(65,more=True),'middle_gate_ball_take_lower_endpoint':gate.str(65,more=True)})
weakg=arb(1)/65000;required=32*weakg**q/(arb(10)/3)**q
out={'status':'six_conditional_budgets_not_an_operator_bound','source_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),'B':B.str(65,more=True),'rows':rows,'I_required_with_existing_g_1_over_65000':required.str(65,more=True),'elapsed_seconds':time.monotonic()-started,'precision_bits':ctx.prec}
(p/'round006/new_shorttime_budget.json').write_text(json.dumps(out,indent=2)+'\n')
print('B',B);print('eps1/4096,T4 gate',rows[3]['middle_gate_ball_take_lower_endpoint']);print('required I with weak g',required)
