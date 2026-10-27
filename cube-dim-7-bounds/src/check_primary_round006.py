#!/usr/bin/env python3
"""Primary exact saved bound interfaces and four-moment algebra; read only."""
from pathlib import Path
from fractions import Fraction as Q
import json,hashlib
import sympy as s
P=Path(__file__).resolve().parents[1]
def dy(x):return Q(int(x[0]))*Q(2)**int(x[1])
def main():
 c=json.loads((P/'upper/round006/radau4_mellin/certificate.json').read_text());lo,hi=dy(c['I_upper']['lower']),dy(c['I_upper']['upper']);cap=Q(7723627498344,10**18)
 assert 0<lo<=hi<cap<dy(c['R005_exact_budget_expression_lower']['lower'])
 n=json.loads((P/'lower/round006/results.json').read_text());old=json.loads((P/'lower/round003/results.json').read_text());lower,upper=map(Q,n['full_construction_comparison_interval']);finite=Q(n['finite_L1']['exact_ratio_lower']);oldhi=Q(old['full_construction_comparison_interval'][1])
 assert Q(n['complete_round003_upper'])==oldhi and upper>=lower>finite>oldhi
 assert finite>Q(2136578727666,10**12)
 state=json.loads((P/'lower/round006/continuation_state.json').read_text());assert state['this_round_evaluated_indices']==[52,38]and state['remaining_new_events']==57
 m1,m2,m3,m4=s.symbols('m1 m2 m3 m4');D=m1*m3-m2*m2;alpha=(m2*m3-m1*m4)/D;beta=(m2*m4-m3*m3)/D
 assert s.cancel(m3+alpha*m2+beta*m1)==0 and s.cancel(m4+alpha*m3+beta*m2)==0
 w=s.symbols('w');assert s.expand(s.chebyshevt(8,1+2*w*(w+2))-s.chebyshevt(16,1+w))==0
 r1,r2,r3,r4=s.symbols('r1 r2 r3 r4');x=s.symbols('x')
 R=[s.Integer(1),r1,r2,r3,r4];raw=[sum(s.binomial(j,k)*R[abs(j-2*k)]for k in range(j+1))/2**j for j in range(5)]
 moments=[s.expand(sum(s.binomial(j,k)*(-1)**(j-k)*raw[k]for k in range(j+1)))for j in range(5)]
 assert moments[4]==r4/8-r3+7*r2/2-7*r1+s.Rational(35,8)
 g=Q(1,219000);cost6=Q(7)**7*hi**6/g**7;assert cost6<128**6
 print(json.dumps(dict(status='PASS primary exact bound comparisons and orthogonality/Chebyshev/four-moment identities',I_finer_upper=str(hi),I_rational_upper=str(cap),new_strict_C7_lower='2136578727666/1000000000000',finite_L1_margin_vs_complete_R003=str(finite-oldhi),conditional_g_target=str(g),conditional_cost_sixth_power_exact=str(cost6),global_g_proved=False,new_C7_upper=False,lower_certificate_sha256=hashlib.sha256((P/'lower/round006/results.json').read_bytes()).hexdigest(),upper_certificate_sha256=hashlib.sha256((P/'upper/round006/radau4_mellin/certificate.json').read_bytes()).hexdigest()),indent=2))
if __name__=='__main__':main()
