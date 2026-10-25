from flint import arb,ctx
from fractions import Fraction as F
from pathlib import Path
import json,hashlib
ctx.prec=256
A=lambda x:arb(x.numerator)/x.denominator
pi=arb.pi(); C0=arb(520066141)/(arb(3121246686412800)*pi**6)
Cf=arb(32)/(243*pi**6)
assert C0/Cf < A(F(1,256))**2
# For t>=1/256>tau, the analytic budget integrand is decreasing.
# Right endpoint rectangles provide a LOWER bound on the budget B,
# NOT a lower bound for the actual Green norm.
rows=[]; S=arb(0)
for j in range(-8,0):
 a=F(2)**j;b=2*a
 term=A(b-a)*C0**(arb(1)/6)*A(b)**(-arb(1)/6)/(pi*(1+1024*A(b)**2))**(arb(5)/6)
 S+=term;rows.append({'a':str(a),'b':str(b),'lower_rectangle':str(term)})
U=A(F(349,10**7)) # 0.0000349, above independently certified exact point
source=Path(__file__).resolve().parents[2]/'round002/green_point/laurent_points_results.json'
point=json.loads(source.read_text())['rows'][0]['G_exact']
coeff=F(point.split(')/pi^5')[0][1:])
assert A(coeff)/pi**5 < U
assert S>arb(16)/3*U
out={'status':'PASS_budget_obstruction_only_NOT_model_exclusion','B_lower':str(S),'16_over3_U':str(arb(16)/3*U),'U':'349/10000000','rectangles':rows,'point_input_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
Path(__file__).with_suffix('.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
