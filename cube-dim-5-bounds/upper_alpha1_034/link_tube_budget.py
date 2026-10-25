#!/usr/bin/env python3
"""Read locked037 cap; exact-rational radius/Holder error, no moments replayed."""
from fractions import Fraction as F
from pathlib import Path
import hashlib
import json

from check_obstruction import root, encode

HERE = Path(__file__).resolve().parent


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def volume(delta):
    assert 0 <= delta <= F(1, 3)
    return 108*delta**3 - 243*delta**4


def main():
    reg = json.loads((HERE / 'tube_link_registration.json').read_text())
    for row in reg['sources']:
        assert sha(Path(row['source'])) == row['sha256']
    source = json.loads(Path(reg['sources'][1]['source']).read_text())
    assert source['parameters']['r'] == '3/2'
    assert source['parameters']['gamma'] == '1/6'
    assert source['parameters']['trace_k'] == '1/3'
    assert source['rational_K_checks']['1/3840']
    cap = F(source['strict_rational_K_upper'])
    assert cap == F(1, 3840)
    allowance = F(reg['fixed_error'])
    k = 2
    while cap**5 * volume(F(1, 2**k)) > allowance**6:
        k += 1
    delta = F(1, 2**k)
    v = volume(delta)
    assert cap**5 * v <= allowance**6
    assert cap**5 * volume(2*delta) > allowance**6
    error = root(cap**5 * v, 6, digits=100)
    assert error[1] < allowance
    out = {
        'status': 'locked037_to_exact15tube_Holder_budget_passed',
        'script_sha256': sha(Path(__file__)),
        'registration_sha256': sha(HERE / 'tube_link_registration.json'),
        'consumed_037_as_external_proved_input_not_independently_recomputed': True,
        'strict_L_3over2_integral_cap': str(cap),
        'target_tube_Iq_error': str(allowance),
        'chosen_dyadic_exponent': k,
        'chosen_delta': str(delta),
        'exact_tube_union_volume': str(v),
        'tube_Iq_error_upper_by_rational_cap': encode(error),
        'previous_dyadic_radius_fails_this_budget_envelope': True,
        'interpretation': 'A sufficient full-tube error budget, not true tube mass or a necessary geometric radius.',
        'new_Green_or_moment_evaluations': 0,
        'new_global_I_or_C5_bound': False,
        'external_sources': reg['sources'],
    }
    (HERE / 'tube_link_results.json').write_text(json.dumps(out, indent=2)+'\n')
    print(json.dumps({'status': out['status'], 'delta': str(delta),
                      'volume': str(v), 'error_upper': out['tube_Iq_error_upper_by_rational_cap']['decimal_outward']}, indent=2))


if __name__ == '__main__':
    main()
