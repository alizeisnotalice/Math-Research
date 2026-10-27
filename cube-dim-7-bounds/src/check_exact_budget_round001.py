#!/usr/bin/env python3
"""Outward rational endpoints for the full analytic I budget; no g proof."""
import hashlib
import json
import time
from fractions import Fraction as Q
from pathlib import Path

from flint import arb, ctx


def main():
    start = time.perf_counter()
    ctx.prec = 256
    pi = arb.pi()
    moment = arb(1)/4*((arb(65)/128).gamma()/(arb(1)/2).gamma())**64
    coefficient = 1/(105*arb(2)**(arb(25)/6)*pi**(arb(13)/12))
    integral = coefficient*moment**(arb(7)/6)
    ml, mu = map(lambda a: str(a.fmpq()), (moment.lower(), moment.upper()))
    il, iu = map(lambda a: str(a.fmpq()), (integral.lower(), integral.upper()))
    assert 0 < Q(ml) < Q(mu) and 0 < Q(il) < Q(iu)
    target = Q(1, 179000)
    cost6 = 7**7*Q(iu)**6/target**7
    assert cost6 < 128**6
    print(json.dumps({
        'status': 'PASS complete analytic I upper budget only',
        'arb_bits': ctx.prec,
        'M_Gamma_exact_expression': '(Gamma(65/128)/Gamma(1/2))^64/4',
        'M_Gamma_expression_outward_interval': [ml, mu],
        'I_upper_expression': 'M_Gamma^(7/6)/(105*2^(25/6)*pi^(13/12))',
        'I_upper_expression_outward_interval': [il, iu],
        'I_certified_upper': iu,
        'conditional_g_target': str(target),
        'conditional_cost_sixth_power_exact': str(cost6),
        'conditional_gap_below_128_sixth_power': str(128**6-cost6),
        'global_g_target_proved': False,
        'new_C7_upper_bound': False,
        'new_Green_nodes': 0, 'new_random_samples': 0,
        'checker_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'seconds': time.perf_counter()-start,
    }, indent=2))


if __name__ == '__main__':
    main()
