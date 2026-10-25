#!/usr/bin/env python3
"""No Arb or Gamma: recheck the recorded integer-power threshold certificates."""
from fractions import Fraction as F
from pathlib import Path
import hashlib
import json

HERE = Path(__file__).resolve().parent


def V(t):
    assert 0 <= t <= F(1, 3)
    return 27*t**3*(4-9*t)


def root_interval(q, degree, digits=40):
    """Exact rational enclosure of a positive root; no floating point."""
    scale = 10**digits
    z = (q*scale**degree).__floor__()
    lo, hi = 0, 1 << ((z.bit_length()+degree-1)//degree)
    assert hi**degree > z
    while hi-lo > 1:
        mid = (lo+hi)//2
        if mid**degree <= z:
            lo = mid
        else:
            hi = mid
    lower, upper = F(lo,scale), F(hi,scale)
    assert lower**degree <= q < upper**degree
    return {'lower_fraction':str(lower), 'upper_fraction':str(upper)}


def main():
    data = json.loads((HERE/'results.json').read_text())
    assert data['precision_bits'] == 384
    assert data['old_r_3over2_recomputed'] is False
    rows = data['rows']
    assert [r['r'] for r in rows] == ['8/5','13/8','33/20']
    summary = []
    for row in rows:
        r = F(row['r']); gamma = (r-1)/(2*r); k = 5*gamma-F(1,2)
        assert gamma == F(row['gamma']) and k == F(row['trace_k'])
        assert F(5,4) < r < F(5,3) and 0 < k < F(1,2)
        a,b,d = row['integer_power_test_A_B_D']
        assert F(a,d) == F(5,4)/r and F(b,d) == 1-F(a,d)
        cap, low = map(F, [row['strict_rational_K_cap'],row['strict_rational_K_lower']])
        assert 0 < low <= F(row['K']['lower_dyadic'])
        assert F(row['K']['upper_dyadic']) < cap
        delta = F(row['delta']); n = row['maximum_certified_dyadic_radius_exponent']
        assert delta == F(1,2**n) and n > 2
        assert V(delta) == F(row['tube_volume'])
        rhs = F(1, 10**(6*d)); assert rhs == F(row['target_rhs'])
        pass_ratio = cap**a*V(delta)**b/rhs
        fail_ratio = cap**a*V(2*delta)**b/rhs
        lower_fail_ratio = low**a*V(2*delta)**b/rhs
        assert pass_ratio <= 1 < lower_fail_ratio <= fail_ratio
        assert pass_ratio == F(row['exact_pass_ratio_lhs_over_rhs'])
        assert fail_ratio == F(row['exact_previous_fail_ratio_lhs_over_rhs'])
        assert lower_fail_ratio == F(row['exact_previous_lower_fail_ratio'])
        assert delta/F(data['old_delta']) == F(row['improvement_factor_over_locked_2pow_minus23'])
        # Every saved binary comparison is rechecked from its integer exponent.
        for trial in row['binary_exponent_trials']:
            assert (cap**a*V(F(1,2**trial['n']))**b <= rhs) is trial['passes']
        error = root_interval(cap**a*V(delta)**b, d)
        previous_lower_error = root_interval(low**a*V(2*delta)**b, d)
        assert F(error['upper_fraction']) < F(1,10**6)
        assert F(previous_lower_error['lower_fraction']) > F(1,10**6)
        summary.append({'r':row['r'], 'delta':str(delta),
                        'envelope_error_at_delta':error,
                        'envelope_error_at_previous_delta_using_K_lower':previous_lower_error})
    # Both first candidates permit the same dyadic radius. Compare their
    # envelopes at that radius without extracting any roots.
    x = F(rows[0]['strict_rational_K_cap'])**25 * V(F(rows[0]['delta']))**7
    y = F(rows[1]['strict_rational_K_cap'])**10 * V(F(rows[1]['delta']))**3
    first_better_at_shared_delta = x**13 < y**32
    out = {'status':'all_three_exact_rational_threshold_certificates_verified',
           'results_sha256':hashlib.sha256((HERE/'results.json').read_bytes()).hexdigest(),
           'verifier_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
           'Gamma_or_Arb_calls':0, 'rows':summary,
           'r_8over5_has_smaller_error_than_13over8_at_shared_delta':first_better_at_shared_delta}
    (HERE/'rational_verification.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out, indent=2))


if __name__ == '__main__':
    main()
