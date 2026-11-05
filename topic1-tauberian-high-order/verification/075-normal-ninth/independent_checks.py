"""Finite exact bookkeeping checks, not model coefficients or geometric proof."""
from fractions import Fraction as F
from itertools import product
import json

grid = (F(0), F(1, 4), F(1, 2), F(1), F(2))
checks = 0
for top, normal, bad, erosion in product(grid, repeat=4):
    clip = lambda x: min(top, max(F(0), x))
    nminus = clip(normal - erosion)
    aminus = clip(normal + bad - erosion)
    assert 0 <= aminus - nminus <= bad
    assert 0 <= clip(normal) - nminus <= erosion * (normal > 0)
    checks += 1

clip = lambda x: min(F(1), max(F(0), x))
normal_loss = clip(F(1)) - clip(F(1)-F(1, 2))
full_loss = clip(F(2)) - clip(F(2)-F(1, 2))
assert normal_loss == F(1, 2) and full_loss == 0

# All physical errors below carry one common factor a.
# Keys are (power of C, power of M, power of R, power of v).
def plus(a, b):
    out = a.copy()
    for key, value in b.items():
        out[key] = out.get(key, F(0)) + value
    return out

def multiply_lambda(a):
    return {(c+1, m, r, v+1): value for (c, m, r, v), value in a.items()}

def shift_v(a, power):
    return {(c, m, r, v+power): value for (c, m, r, v), value in a.items()}

rho = {(0, 0, 1, 8): F(1)}
mass = {(0, 1, 0, 5): F(1)}
e0 = plus(rho, multiply_lambda(mass))
e1 = plus(rho, multiply_lambda(e0))
e2 = plus(rho, multiply_lambda(e1))
bdy3 = multiply_lambda(e2)
source = shift_v(bdy3, 1)
refund = shift_v(e2, 2)
assert min(key[3] for key in e2) == 8
assert min(key[3] for key in bdy3) == 9
assert min(key[3] for key in source) == 10
assert min(key[3] for key in refund) == 10
assert source != refund
assert shift_v(refund, 0) != source
# Source and target error polynomials differ by C, not by a power of v.
assert source == {(c+1, m, r, v): value for (c,m,r,v),value in refund.items()}

def describe(poly):
    return [{'C_power': c, 'M_power': m, 'R_power': r, 'v_power': v,
             'coefficient': str(value)} for (c,m,r,v),value in sorted(poly.items(), key=lambda x:x[0][3])]

result = {
    'all_pass': True,
    'scalar_exact_cases': checks,
    'clip_same_eroder_bound': True,
    'normal_loss_bounded_by_integrated_eroder': True,
    'full_clipped_loss_does_not_bound_normal_loss':
        {'top': '1', 'normal': '1', 'bad': '1', 'eroder': '1/2',
         'normal_loss': str(normal_loss), 'full_loss': str(full_loss)},
    'D_minus_Y2_after_a': describe(e2),
    'B_D_minus_Y3_after_a': describe(bdy3),
    'source_price': describe(source),
    'refund_price': describe(refund),
    'source_and_refund_same_order_not_same_polynomial': True,
    'bad_mass_to_root_weight_orders': {'M_bad': 9, 'edge_lower_cost': 6, 'root_weight': 3},
    'bad_layer_refund_order': 7+3,
    'bad_raw_refund_order': 9+1,
    'clip_order': 9+3,
    'scope': 'Finite scalar and error-order checks only; no model coefficients and no C_n ninth coefficient.'
}
print(json.dumps(result, ensure_ascii=False, indent=2))
