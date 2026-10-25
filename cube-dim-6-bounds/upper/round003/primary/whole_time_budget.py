import mpmath as m
import json,hashlib
from pathlib import Path
m.mp.dps=60
C0=m.mpf(520066141)/(m.mpf(3121246686412800)*m.pi**6)
Cf=m.mpf(32)/(243*m.pi**6)
tau=m.sqrt(C0/Cf)
q=m.mpf(6)/5
# Positive orthant factor32^-5/6 exactly cancels marginal's32^5/6.
def integrand(t):
    peak=min(C0/t**6,Cf/t**4)
    return peak**(m.mpf(1)/6)*(t/(m.pi*(1+1024*t*t)))**(m.mpf(5)/6)
parts=[m.quad(integrand,[0,tau]),m.quad(integrand,[tau,m.mpf(1)/32]),m.quad(integrand,[m.mpf(1)/32,m.inf])]
B=sum(parts)
U=m.mpf('0.000034897358069585267267294466548812617323309243852628340274667865970046219614897077')
out={'status':'numerical_evaluation_of_proved_analytic_norm_upper_bound_not_interval_certificate','C_peak':str(C0),'C_face':str(Cf),'tau':str(tau),'pieces':list(map(str,parts)),'B_norm_upper_expression_numeric':str(B),'conditional_cost_using_point_U_NOT_a_valid_C6_bound':str(6**q*(B/U)**q),'required_g_for_this_budget_to_beat64':str(m.mpf(3)/16*B),'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
Path(__file__).with_suffix('.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
