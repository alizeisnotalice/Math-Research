#!/usr/bin/env python3
"""Finite new face-tube geometry and one conditional budget; no Green calls."""
from fractions import Fraction as F
from itertools import combinations, product
from pathlib import Path
import hashlib
import json

from flint import arb, ctx
import sympy as sp

HERE = Path(__file__).resolve().parent
ctx.prec = 256


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def af(x):
    x = F(x)
    return arb(x.numerator) / x.denominator


def saved_ball_upper(s):
    center, radius = map(F, s[1:-1].split(' +/- '))
    return center + radius


def main():
    reg = json.loads((HERE / 'cover_registration.json').read_text())
    for row in reg['sources']:
        assert sha(Path(row['source'])) == row['sha256']
    signs = list(product((-1, 1), repeat=4))
    assert len(signs) == 16
    classes = {'negative_fixed_coordinate': 0, 'positive_diagonal_line': 0,
               'mixed_variable_signs_contained_in_diagonal_tube': 0}
    dimension_counts = {d: {'diagonal': 0, 'mixed': 0, 'negative_fixed': 0}
                        for d in range(1, 5)}
    for v, w in combinations(signs, 2):
        D = [j for j in range(4) if v[j] != w[j]]
        fixed = [j for j in range(4) if v[j] == w[j]]
        assert D
        if any(v[j] == -1 for j in fixed):
            classes['negative_fixed_coordinate'] += 1
            dimension_counts[len(D)]['negative_fixed'] += 1
        elif len({v[j] for j in D}) == 1:
            classes['positive_diagonal_line'] += 1
            dimension_counts[len(D)]['diagonal'] += 1
        else:
            assert {v[j] for j in D} == {-1, 1}
            classes['mixed_variable_signs_contained_in_diagonal_tube'] += 1
            dimension_counts[len(D)]['mixed'] += 1
    assert classes == {'negative_fixed_coordinate': 80,
                       'positive_diagonal_line': 15,
                       'mixed_variable_signs_contained_in_diagonal_tube': 25}

    # Exact range-volume formula in a k-dimensional interval of length L:
    # vol(range <= a) = k*L*a**(k-1) - (k-1)*a**k, 0<=a<=L.
    # In the union of 15 diagonal tubes, classify coordinates below 1-delta.
    delta = sp.symbols('delta', nonnegative=True)
    strata = {'0': delta**4}
    for k in range(1, 5):
        strata[str(k)] = sp.expand(
            sp.binomial(4, k) * delta**(4-k) *
            (k*(1-delta)*(2*delta)**(k-1) - (k-1)*(2*delta)**k))
    union = sp.expand(sum(strata.values()))
    assert union == 108*delta**3 - 243*delta**4
    individual_sum = sp.expand(sum(
        sp.binomial(4, k)*delta**(4-k)*
        (k*(2*delta)**(k-1) - (k-1)*(2*delta)**k)
        for k in range(1, 5)))
    assert individual_sum == 108*delta**3 - 136*delta**4
    assert union.subs(delta, sp.Rational(1, 3)) == 1

    # Read a locked complete old short-time enclosure, not recomputed Gamma.
    old = json.loads(Path(reg['sources'][0]['source']).read_text())
    candidates = [r for r in old['rows'] if r['eps'] == '1/512']
    assert len({r['shorttime_norm_budget'] for r in candidates}) == 1
    e_short = saved_ball_upper(candidates[0]['shorttime_norm_budget'])
    goal, e_tube, h_cap = F(1, 20000), F(1, 1000000), F(1, 25000)
    q, iq = af(F(5, 4)), af(F(4, 5))
    resulting_I = (af(h_cap)**iq + af(e_short))**q + af(e_tube)
    h_acceptance = ((af(goal-e_tube))**iq - af(e_short))**q
    assert resulting_I < af(goal)
    assert af(h_cap) < h_acceptance
    assert af(goal) < af(F(53, 1000000))
    out = {
        'status': 'new_exact_face_tube_geometry_and_conditional_budget_passed',
        'script_sha256': sha(Path(__file__)),
        'registration_sha256': sha(HERE / 'cover_registration.json'),
        'pair_count': 120, 'classes': classes, 'dimension_counts': dimension_counts,
        'cover': '15 lines L_D, D nonempty; y_j=t on D and y_j=1 off D',
        'tube_metric': 'linfinity distance in the affine face x1=1',
        'exact_union_volume_delta_le_1_over3': str(union),
        'individual_volume_union_bound_delta_le_1_over2': str(individual_sum),
        'disjoint_order_statistic_strata': {k: str(v) for k, v in strata.items()},
        'fullspace_distance_relation': 'outside affineface delta tube implies dist_infinity((1,y),S2)>delta/2',
        'higher_moment_required': 'A >= integral_[0,1]^4 G(1,y)^(3/2), to be independently supplied by037',
        'tube_Iq_bound': 'A^(5/6)*(108*delta^3-243*delta^4)^(1/6), delta<=1/3',
        'sufficient_delta_for_tube_allowance': 'delta<=min(1/3,(e_tube^6/(108*A^5))^(1/3))',
        'one_future_certificate_budget': {
            'eps': '1/512', 'I_goal': str(goal), 'tube_I_allowance': str(e_tube),
            'shorttime_norm_upper_fraction_from_locked_ball': str(e_short),
            'regular_H_eps_integral_cap': str(h_cap),
            'combined_conditional_I_upper_ball': str(resulting_I),
            'maximal_regular_H_eps_integral_acceptance_ball': str(h_acceptance),
            'strictly_below_goal_if_prerequisites_hold': True,
        },
        'not_done': ['No037 moment recomputation', 'No new Green point',
                     'No Taylor boxes integrated', 'No continuum quadrature error certified',
                     'No new I or C5 bound'],
    }
    (HERE / 'cover_results.json').write_text(json.dumps(out, indent=2) + '\n')
    print(json.dumps({'status': out['status'], 'classes': classes,
                      'volume': str(union), 'conditional_I_upper': str(resulting_I),
                      'H_integral_acceptance': str(h_acceptance)}, indent=2))


if __name__ == '__main__':
    main()
