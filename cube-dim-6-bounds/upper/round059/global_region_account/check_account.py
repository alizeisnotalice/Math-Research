"""Portable, saved-rational budget arithmetic. No kernels or private paths."""
import argparse
import datetime
import hashlib
import itertools
import json
import math
import platform
import time
from fractions import Fraction as Q
from pathlib import Path


def enclosure(x, n, denominator=10**60):
    assert x >= 0 and n > 0
    target = x.numerator * denominator**n
    divisor = x.denominator
    lo, hi = 0, 1
    while hi**n * divisor <= target:
        hi *= 2
    while hi - lo > 1:
        mid = (lo + hi) // 2
        if mid**n * divisor <= target:
            lo = mid
        else:
            hi = mid
    assert lo**n * divisor <= target < hi**n * divisor
    return Q(lo, denominator), Q(hi, denominator)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output-dir', type=Path, required=True)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=False)
    tic = time.monotonic()
    start = datetime.datetime.now(datetime.timezone.utc).isoformat()
    certificates = []

    def root(x, n, label):
        lower, upper = enclosure(x, n)
        certificates.append(dict(label=label, radicand=str(x), degree=n,
                                 lower=str(lower), upper=str(upper)))
        return lower, upper

    counts = [math.comb(5, k-1) for k in range(1, 6)]
    assert counts == [1, 5, 10, 10, 5] and sum(counts) == 31
    assignments = []
    for c, e, d, h in itertools.product([False, True], repeat=4):
        labels = [c, e and not c, d and not(c or e),
                  h and not(c or e or d), not(c or e or d or h)]
        assert sum(labels) == 1
        assignments.append(dict(C=c, E=e, D=d, H=h, part=labels.index(True)))
    assert 0 < Q(1, 4) < Q(3, 4) < 1

    b_c = Q(2118205441, 10**15)
    b_e = Q(13404268290817, 10**18)
    b_d_short = Q(4689838525803, 10**18)
    c_i = root(b_c**6, 5, 'cap_B_to_I')
    e_i = root(b_e**6, 5, 'low_B_to_I')
    paid_i = c_i[1] + e_i[1]
    paid_b = root(paid_i**5, 6, 'paid_disjoint_regions_B')
    b_crit = Q(4, 48750)
    crit_i = root(b_crit**6, 5, 'conditional_target_I')
    room = crit_i[0] - paid_i
    assert room > 0
    room_b = root(room**5, 6, 'sufficient_remaining_B_allowance')
    # Exact algebraic equivalence of the two cost comparison input conventions.
    assert (6*b_crit)**6 == 64**5 * Q(1, 65000)**6

    result = dict(
        status='PASS_PARTITION_AND_PARTIAL_BUDGET_ARITHMETIC',
        configurations=5, counts=counts, boolean_partition_cases=assignments,
        exact_region_proof='theorem.md', root_certificates=certificates,
        paid_cap_B=str(b_c), paid_low_B=str(b_e), middle_only_short_B=str(b_d_short),
        paid_I_conservative_upper=str(paid_i),
        paid_disjoint_B_upper=str(paid_b[1]),
        conditional_Bcrit=str(b_crit), conditional_global_g='1/65000',
        sufficient_unpaid_I_allowance=str(room),
        sufficient_unpaid_B_allowance_lower=str(room_b[0]),
        unpaid_terms=['middle long-time norm L_D', 'high-endpoint integral I_H',
                      'far-region integral I_F', 'global Green floor g'],
        global_I_certified=False, global_g_certified=False, new_C6_bound=False,
        kernel_or_native_calls=0, source_production_replays=0,
        start_utc=start, end_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
        elapsed_seconds=time.monotonic()-tic, python=platform.python_version(),
        checker_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    (args.output_dir/'results.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(dict(status=result['status'], configurations=5,
                         paid_B_display=float(paid_b[1]),
                         remaining_B_sufficient_display=float(room_b[0]),
                         global_I=False, global_g=False, new_C6=False)))


if __name__ == '__main__':
    main()
