"""Exact certificate for a restricted kernel-family obstruction, not for C_2."""
from fractions import Fraction as F
from pathlib import Path
import json
from upper_certificate import I, logq

r=F(86425249,10**8)
T=2*(2-r)*(1-logq(2)+logq(r))-2*(1-r)*logq(1-r)
Tp=2*logq(2)-2*logq(r)+2*logq(1-r)+4/r-2
b=2-2*logq(2)
endpoint=b-r**3*T
slope=F(27,2)*r-4*Tp
assert T.a>0 and endpoint.a>0 and slope.a>0
R=F(4746319,2500000)
s=2*r/R
J=I(F(1,2)-2/s)-4/s**2*logq(1-s/2)
floor=R**2/(R-1)-J/(6*R*(R-1)*T)
incumbent=F(27226955153507,7019746875000)
gap=incumbent-floor
out={'r0':str(r),'T':T.pair(),'endpoint_margin':endpoint.pair(),
     'negative_Hprime_margin':slope.pair(),'fixed_R':str(R),
     'restricted_family_cost_floor':floor.pair(),
     'incumbent_minus_floor':gap.pair(),
     'gap_upper_display':float(gap.b),
     'status':'exact_certificate_for_nonnegative_scaled_bump_sufficient_condition_model_only'}
(Path(__file__).resolve().parents[1]/'results/round4_obstruction.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps({k:out[k] for k in ('status','gap_upper_display')},indent=2))
