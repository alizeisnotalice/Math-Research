"""Exact rational threshold bookkeeping for the analytic cross-review."""
from fractions import Fraction as F
import json

z = F(1, 3)
N = 4
log2_lower = 2 * sum(z ** (2 * k + 1) / (2 * k + 1) for k in range(N))
log2_upper = log2_lower + 2 * z ** (2 * N + 1) / ((2 * N + 1) * (1 - z*z))
assert F(1, 2) < log2_lower < log2_upper < F(7, 10)
L = F(32)
logL_upper = 5 * F(7, 10)
checks = {
    'half_denominator_margin_at_L32': L/2 - 3*logL_upper,
    'x0_below_2L_margin_at_L32': L - 4*logL_upper,
    'derivative_half_denominator_lower_bound': F(1, 2) - 3/L,
    'derivative_x0_comparison_lower_bound': 1 - 4/L,
    'derivative_expL_over_L2_lower_bound': 1 - 2/L,
    'exp32_lower_bound_minus_2L2': F(2**32) - 2*L*L,
    'low_x_region_surplus_per_n_lower_bound': F(1, 2) - 4/L,
    'L3_minus_e_upper_bound': L**3 - 3,
    'sqrt2_upper_square_margin': F(3, 2)**2 - 2,
    'small_parameter_domain_exp32_minus_2L': F(2**32) - 2*L,
}
assert all(v > 0 for v in checks.values())
print(json.dumps({
    'status': 'all exact rational checks passed',
    'log2_lower': str(log2_lower), 'log2_upper': str(log2_upper),
    'strict_positive_margins': {k: str(v) for k, v in checks.items()},
    'scope': 'threshold arithmetic only; monotonicity and arbitrary-n estimates are proved analytically in cross-review'
}, indent=2))
