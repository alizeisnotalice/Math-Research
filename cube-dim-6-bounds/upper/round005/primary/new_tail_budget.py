from fractions import Fraction as Q
from pathlib import Path
import json,hashlib
from flint import arb,ctx
ctx.prec=384
N=32;assert Q(N*N,12**5)==Q(1,243)
K=(arb(6)/7)*(arb(1)/243)**(arb(1)/6)*arb.pi()**(-arb(11)/6)
q=arb(6)/5
A=lambda f:arb(f.numerator)/f.denominator
rows=[]
for e,g,H in [(Q(1,1000),Q(1,60000),Q(1,16000)),(Q(3,2000),Q(1,60000),Q(1,15000)),(Q(1,500),Q(1,60000),Q(1,17000)),(Q(1,400),Q(1,50000),Q(1,16000))]:
 tail=K*A(e)**(arb(7)/6);cost=6**q*((A(H)+tail)/A(g))**q
 rows.append({'epsilon_I':str(e),'requires_g_ge':str(g),'requires_H_norm_le':str(H),'tail':str(tail),'conditional_cost':str(cost),'strictly_below64':bool(cost<64),'global_premises_certified':False})
out={'status':'PASS_exact_new_tail_constant_and_conditional_arithmetic','Cf_coefficient_of_pi_minus6':'1/243','K':str(K),'Knew_over_Kold':str(arb(32)**(-arb(1)/6)),'conditional_targets':rows,'scope':'all global g/H premises remain unproved; no new C6 bound','source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
Path(__file__).with_suffix('.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
