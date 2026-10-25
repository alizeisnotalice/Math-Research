"""Symbolic knot-limit identities only, no extra spatial sample points."""
from pathlib import Path
import sympy as s,json
OUT=Path(__file__).resolve().parent
e=s.symbols('epsilon',positive=True);b=s.symbols('b',real=True,nonzero=True)
Lb=s.log(b*b)/2
onezero=(s.log(e)-s.log(e*e+b*b)/2+b*b/(e*e+b*b))/(b*b)
twozero=(s.log(e*e+b*b)/2-s.log(e))/(b*b)
triplenonzero=(e*e-b*b)/(2*(e*e+b*b)**2)
assert s.simplify(s.limit(onezero-s.log(e)/(b*b),e,0,dir='+')-(1-Lb)/(b*b))==0
assert s.simplify(s.limit(twozero+s.log(e)/(b*b),e,0,dir='+')-Lb/(b*b))==0
assert s.limit(triplenonzero,e,0,dir='+')==-1/(2*b*b)
c=s.symbols('c',real=True,nonzero=True)
assert s.factor(1/(b*c)+1/(b*(b-c))+1/(c*(c-b)))==0
r={'status':'PASS exact coalescent/zero-knot algebra','L_b':'log|b|=log(b^2)/2','rows':[{'nodes':'0,b,b','log_epsilon_coefficient':'1/b^2','finite_real_constant':'(1-L_b)/b^2'},{'nodes':'0,0,b','log_epsilon_coefficient':'-1/b^2','finite_real_constant':'L_b/b^2'},{'nodes':'b,b,b; b!=0','log_epsilon_coefficient':'0','finite_real_constant':'-1/(2b^2)'},{'nodes':'0,b,c; b,c nonzero distinct','log_epsilon_coefficient':'1/(bc)','finite_real_constant':'L_b/[b(b-c)]+L_c/[c(c-b)]'}],'scope':'formal identities, no additional tau/spatial numerical points beyond preregistered36'}
(OUT/'cb11_knot_limits.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(r,indent=2))
