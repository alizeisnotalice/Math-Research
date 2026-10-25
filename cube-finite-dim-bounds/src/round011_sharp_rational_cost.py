"""Exact final-cost arithmetic only; analytic g/I certificates are separate.

No nodes, optimization, baseline recomputation or ledger mutation.
Machin's identity and alternating sums certify the rational pi upper bound.
"""
from fractions import Fraction as F
import json


def atan_sum(x, terms):
    return sum(((-1)**k * x**(2*k+1)/F(2*k+1)
                for k in range(terms)), F(0))


# tan(4 atan(1/5))=120/119, and the subtraction tangent is 1.
# 0 < 4 atan(1/5)-atan(1/239) < 4/5 < pi/2 fixes the branch.
t2 = 2*F(1,5)/(1-F(1,5)**2)
t4 = 2*t2/(1-t2**2)
assert (t4-F(1,239))/(1+t4*F(1,239)) == 1
pi_upper_machin = 16*atan_sum(F(1,5),13)-4*atan_sum(F(1,239),4)
pi_upper = F(355,113)
assert pi_upper_machin < pi_upper

g3, i3, u3 = F(2621323989,500000000000), F(23129127,40000000000), F(1979,250)
i4, u4 = F(92072351,500000000000), F(14371,1000)
power3 = 27*i3**2/g3**3
# g4=17/(640*pi^3); bounding pi from above bounds the cost cube above.
power4 = 256*i4**3*(F(640,17)*pi_upper**3)**4
assert power3 < u3**2 < 8**2
assert power4 < u4**3 < 16**3
print(json.dumps({
    'status': 'PASS exact conditional cost comparisons; analytic certification separate',
    'pi_upper': str(pi_upper),
    'machin_to_pi_upper_gap': str(pi_upper-pi_upper_machin),
    'n3': {'candidate_U':str(u3),'cost_squared':str(power3),'strict_power_gap':str(u3**2-power3)},
    'n4': {'candidate_U':str(u4),'cost_cube_rational_upper':str(power4),'strict_power_gap':str(u4**3-power4)},
}, indent=2))
