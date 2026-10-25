"""Exact scalar checks quoted in the n4 paper; no search or new cells."""
from fractions import Fraction as Q
import json

I = Q(92072351, 500000000000)
coefficient = 256 * I**3 * Q(640, 17)**4
assert coefficient == Q(6244212653911047061276408, 1944624818861484527587890625)
def atan_sum(n, x):
    return sum(((-1)**k * x**(2*k+1) / Q(2*k+1) for k in range(n)), Q())
pi_upper = 16*atan_sum(13, Q(1,5)) - 4*atan_sum(4, Q(1,239))
assert pi_upper < Q(355,113)
gap = Q(14371,1000)**3 - coefficient*Q(355,113)**12
assert gap == Q(2439293356294278395819615382243120026281259439,
               17676938664604974026861767968798828125000000000) > 0
assert 32*Q(1,8)+24*Q(1,16) == Q(11,2)
assert 4*8**2*4*Q(11,2) == 5632
assert Q(1,36)+Q(1,12)+Q(1,18) == Q(1,6)
assert min(Q(1,24), Q(17,640), Q(482191,13088544)) == Q(17,640)
print(json.dumps({
    'status': 'PASS actual Python Fraction paper scalar checks',
    'I0': str(I), 'coefficient_pi12': str(coefficient),
    'machin_upper': str(pi_upper), 'cost_gap': str(gap),
    'cap_heat_coefficient': 5632, 'model_minimum_coefficient': '17/640',
    'new_nodes': 0, 'new_integral_cells': 0,
    'scope': 'Scalar check only; continuous certificates require 035/037 and second review.'
}, indent=2))
