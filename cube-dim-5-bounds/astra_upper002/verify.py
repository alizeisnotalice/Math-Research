"""Exactly three registered shape-moment certificates; no stochastic estimates."""
from pathlib import Path
from fractions import Fraction as F
from itertools import combinations
from math import comb
import hashlib
import importlib.util
import json
import sys
import time
from flint import arb, ctx

ROOT = Path(__file__).resolve().parent
SOURCE = ROOT.parent / 'astra_upper001' / 'verify.py'
spec = importlib.util.spec_from_file_location('upper001_verified', SOURCE)
old = importlib.util.module_from_spec(spec)
spec.loader.exec_module(old)


def ar(value):
    value = F(value)
    return arb(value.numerator) / value.denominator


def cofactor_moments():
    supersets = [{} for _ in range(5)]
    total = 0
    nonzero = 0
    for indices in combinations(range(16), 4):
        rows = [old.V[i] for i in indices] + [(1, 1, 1, 1, 1)]
        weight = old.bareiss_det(rows) ** 2
        if not weight:
            continue
        nonzero += 1
        total += weight
        for j in range(5):
            for subset in combinations(indices, j):
                supersets[j][subset] = supersets[j].get(subset, 0) + weight
    assert total == 5 * 16 ** 4
    sums = [sum(w ** 2 for w in row.values()) for row in supersets]
    overlaps = [sum((-1) ** (j - k) * comb(j, k) * sums[j] for j in range(k, 5)) for k in range(5)]
    assert all(x >= 0 for x in overlaps)
    assert sum(overlaps) == total ** 2
    moments = [old.laplace_moment([2] * k + [1] * (8 - 2 * k)) for k in range(5)]
    mean = total * old.laplace_moment([1] * 4)
    second = sum((w * m for w, m in zip(overlaps, moments)), F(0))
    return mean, second, {
        'four_subsets_enumerated': comb(16, 4),
        'nonzero_b_x_corner_terms': nonzero,
        'coefficient_sum': total,
        'ordered_pair_overlap_weights': overlaps,
        'overlap_weight_moments': [str(x) for x in moments],
        'E_b_corner': str(mean), 'E_b_corner_squared': str(second),
    }


def main():
    began = time.monotonic()
    ctx.prec = 256
    m1 = 16 ** 5 * old.laplace_moment([1] * 5)
    assert m1 == F(648, 1463)
    receipt001 = json.loads((ROOT.parent / 'astra_upper001' / 'verification.json').read_text())
    m2 = F(receipt001['finite_determinant_minorant']['E_detM_squared'])
    assert m2 == F(126404864, 474619145)
    b1, b2, cofactor = cofactor_moments()
    c = 3 / (2 ** (arb(15) / 2) * arb.pi() ** (arb(5) / 2))
    tau = m2 / m1
    trace = c * 5 ** (-arb(7) / 4) * ar(m1) ** (arb(5) / 4)
    holder1 = c * ar(m1) ** (arb(17) / 4) / ar(m2) ** (arb(3) / 2) / ar(b1) ** (arb(7) / 4)
    holder2 = c * ar(m1) ** (arb(5) / 2) / ar(m2) ** (arb(5) / 8) / ar(b2) ** (arb(7) / 8)
    common_ceiling = c * ar(m1) ** (arb(5) / 4) / ar(b1) ** (arb(7) / 4)
    assert common_ceiling < arb(1) / 10000
    results = []
    for name, value in [('trace_d_mean', trace), ('holder_b_first_moment', holder1), ('holder_b_second_moment', holder2)]:
        results.append({'name': name, 'g_ball': str(value), 'strictly_above_1_over_65000': bool(value > arb(1) / 65000)})
    assert holder2 > arb(1) / 27000
    assert holder2 > holder1 > trace > arb(1) / 65000
    out = {
        'status': 'three_finite_global_Green_lower_bound_certificates',
        'new_C5_upper_bound': False,
        'exact_inputs': {'E_d': str(m1), 'E_d_squared': str(m2), 'tau': str(tau)},
        'cofactor_receipt': cofactor,
        'candidate_results': results,
        'certified_rational_g': '1/27000',
        'holder_family_ceiling': str(common_ceiling),
        'family_can_reach_g_1e_minus4': False,
        'source_sha256': {
            str(p): hashlib.sha256(p.read_bytes()).hexdigest()
            for p in [SOURCE, ROOT.parent / 'astra_upper001' / 'verification.json', ROOT / 'registration.md', Path(__file__)]
        },
        'runtime': {'python': sys.executable, 'arb_bits': ctx.prec, 'seconds': time.monotonic() - began},
    }
    (ROOT / 'verification.json').write_text(json.dumps(out, indent=2) + '\n')
    print(json.dumps(out, indent=2))


if __name__ == '__main__':
    main()
