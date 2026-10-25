#!/usr/bin/env python3
"""Independent staggered-position certification, using E's own rational helpers.

Positions: x_i=a_i k_i, x_j=a_j(k_j+delta*(-1)^k_i), other x_h=a_h k_h.
Only the i and other-four coordinate reflections are used (32, not 64).
"""
from argparse import ArgumentParser
from fractions import Fraction as Q
from itertools import product
from math import prod
from pathlib import Path
import hashlib
import json
import sys
import time

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from correlated_fraction_checker import (N, Budget, BudgetStop, root6, scale,
    default_spacing, floorq, ceilq, prune, union_integer, finite_transfer)


def physical_extrema(start, lengths, axes, delta):
    """Dimensionless extrema across *all* indices in the Cartesian index block."""
    i, j = axes
    extremes = [(Q(q), Q(q+ell-1)) for q, ell in zip(start, lengths)]
    if lengths[i] == 1:
        shift = delta*(1-2*(start[i] % 2))
        shift_lo = shift_hi = shift
    else:
        shift_lo, shift_hi = -abs(delta), abs(delta)
    extremes[j] = (Q(start[j])+shift_lo, Q(start[j]+lengths[j]-1)+shift_hi)
    return tuple(extremes)


def geometry(candidate, denominator, grid, budget):
    axes = tuple(candidate['axes'])
    assert len(axes) == 2 and len(set(axes)) == 2
    assert all(type(x) is int and 0 <= x < N for x in axes)
    i, j = axes
    delta, level = Q(candidate['delta']), Q(candidate.get('level', '1'))
    assert abs(delta) <= Q(1, 2) and level > 0
    if candidate.get('spacing') == 'unit_example':
        sp = ((Q(1), Q(1)),)*N
    else:
        assert candidate.get('spacing', 'frozen_R6_003') == 'frozen_R6_003'
        sp = default_spacing(denominator)
    p = tuple(2 if h == i else 1 for h in range(N))
    # H_i=[0,a_i], H_j=[-a_j/2,a_j/2], H_other=[0,a_h/2].
    domain = [((Q(0), Q(0)), scale(sp[h], Q(p[h], 2))) for h in range(N)]
    domain[j] = (scale(sp[j], Q(-1, 2)), scale(sp[j], Q(1, 2)))
    inner, outer = [], []
    margin = None
    collar = [0]*N
    for block in candidate['blocks']:
        budget.check()
        q, ell = tuple(block['start']), tuple(block['lengths'])
        assert len(q) == len(ell) == N
        assert all(type(x) is int for x in q)
        assert all(type(x) is int and 1 <= x <= 4 for x in ell)
        mass = Q(prod(ell))
        r0, r1 = scale(root6(mass/level, denominator), Q(1, 2))
        assert level*(2*r0)**6 <= mass < level*(2*r1)**6
        extremes = physical_extrema(q, ell, axes, delta)
        ib, ob = [], []
        physical = []
        for h, ((x, y), a, (domlo, domhi)) in enumerate(zip(extremes, sp, domain)):
            src_lo, src_hi = scale(a, x), scale(a, y)
            physical.append((src_lo, src_hi))
            ilo, ihi = max(domlo[1], src_hi[1]-r0), min(domhi[0], src_lo[0]+r0)
            olo, ohi = max(domlo[0], src_hi[0]-r1), min(domhi[1], src_lo[1]+r1)
            ib.extend((ceilq(ilo*grid)+1, floorq(ihi*grid)-1))
            ob.extend((floorq(olo*grid), ceilq(ohi*grid)))
        if all(ob[2*h] < ob[2*h+1] for h in range(N)):
            outer.append(tuple(ob))
        if not all(ib[2*h] < ib[2*h+1] for h in range(N)):
            continue
        local_margin = None
        for h, ((src_lo, src_hi), (domlo, domhi)) in enumerate(zip(physical, domain)):
            lo, hi = Q(ib[2*h], grid), Q(ib[2*h+1], grid)
            assert domlo[1] < lo < hi < domhi[0]
            gap = min(r0+lo-src_hi[1], r0-hi+src_lo[0])
            assert gap > 0
            local_margin = gap if local_margin is None else min(local_margin, gap)
            s = max(abs(q[h]), abs(q[h]+ell[h]-1))
            collar[h] = max(collar[h], (s+p[h]-1)//p[h])
        margin = local_margin if margin is None else min(margin, local_margin)
        inner.append(tuple(ib))
    return {'inner': prune(inner, budget), 'outer': prune(outer, budget),
            'margin': margin, 'collar': tuple(collar), 'level': level,
            'periods': p, 'spacing': sp}


def certify(candidate, target, denominator=10**50, grid=10**35):
    budget = Budget(nodes=10**6, seconds=20*60)
    geo = geometry(candidate, denominator, grid, budget)
    vi = Q(union_integer(geo['inner'], budget), grid**N)
    vo = Q(union_integer(geo['outer'], budget), grid**N)
    raw_lo, raw_hi = 16*geo['level']*vi, 16*geo['level']*vo
    assert 0 <= raw_lo <= raw_hi
    out = {'status': 'certified_listed_staggered_union',
           'scope': 'listed index blocks only; no complete maximal-set claim',
           'root_denominator': str(denominator), 'grid_denominator': str(grid),
           'input_block_count': len(candidate['blocks']),
           'inner_retained_boxes': len(geo['inner']), 'outer_retained_boxes': len(geo['outer']),
           'listed_union_ratio_interval': [str(raw_lo), str(raw_hi)],
           'target': str(target), 'ratio_lower_minus_target': str(raw_lo-target),
           'geometric_cube_margin': str(geo['margin']),
           'coarse_collar': geo['collar'], 'union_nodes': budget.nodes,
           'node_cap': budget.node_cap, 'wall_seconds_cap': budget.seconds,
           'elapsed_seconds': time.monotonic()-budget.started,
           'above_target': raw_lo > target,
           'outer_proves_listed_candidate_below_target': raw_hi < target}
    if raw_lo > target:
        out['finite_L1_witness'] = finite_transfer(raw_lo, geo['level'], target,
            geo['collar'], geo['margin'], period_mass=2, budget=budget)
    return out


def selfcheck():
    extrema_cases = 0
    reflection_cases = 0
    for axes in ((0, 1), (2, 5)):
        i, j = axes
        for delta in (Q(-1, 3), Q(0), Q(1, 3)):
            for seed in range(12):
                q = tuple((seed+2*h) % 5-2 for h in range(N))
                ell = tuple((seed+h) % 3+1 for h in range(N))
                atoms = []
                for k in product(*(range(q[h], q[h]+ell[h]) for h in range(N))):
                    x = [Q(kh) for kh in k]
                    x[j] += delta*(1-2*(k[i] % 2))
                    atoms.append(x)
                exact = physical_extrema(q, ell, axes, delta)
                assert exact == tuple((min(x[h] for x in atoms), max(x[h] for x in atoms))
                                      for h in range(N))
                extrema_cases += 1
            # Check every allowed coordinate reflection on negative/even/odd columns.
            reflected = tuple(h for h in range(N) if h != j)
            for ki in (-2, -1, 0, 1, 2):
                k = [0]*N
                k[i], k[j] = ki, -1
                x = [Q(kh) for kh in k]
                x[j] += delta*(1-2*(k[i] % 2))
                for sign_bits in product((-1, 1), repeat=5):
                    kp, xp = list(k), list(x)
                    for h, s in zip(reflected, sign_bits):
                        kp[h] *= s
                        xp[h] *= s
                    expected = [Q(kh) for kh in kp]
                    expected[j] += delta*(1-2*(kp[i] % 2))
                    assert xp == expected
                    reflection_cases += 1
            # Translating a full coarse period gives the corresponding physical shift.
            k = [-1]*N
            p = tuple(2 if h == i else 1 for h in range(N))
            base = [Q(x) for x in k]
            base[j] += delta*(1-2*(k[i] % 2))
            for h in range(N):
                kp = list(k)
                kp[h] += p[h]
                expected = [Q(x) for x in kp]
                expected[j] += delta*(1-2*(kp[i] % 2))
                translated = list(base)
                translated[h] += p[h]
                assert expected == translated
    example = {'axes': [0, 1], 'delta': '1/3', 'level': '1', 'spacing': 'unit_example',
               'blocks': [{'start': [0]*N, 'lengths': [1]*N},
                          {'start': [1, 0, 0, 0, 0, 0], 'lengths': [1]*N}]}
    interface = certify(example, Q(1, 10), denominator=10**24, grid=10**12)
    assert interface['above_target'] and interface['finite_L1_witness']['period_mass'] == 2
    assert interface['finite_L1_witness']['strictly_above_target']
    # All 32 allowed reflected blocks, at both extreme receiver translates.
    collar_checks = 0
    p = (2, 1, 1, 1, 1, 1)
    reflected = (0, 2, 3, 4, 5)
    for block in example['blocks']:
        for signs in product((-1, 1), repeat=5):
            q = list(block['start'])
            for h, s in zip(reflected, signs):
                q[h] *= s
            for h, s in enumerate(q):
                ch = interface['coarse_collar'][h]
                for t in (-3, 3):
                    assert p[h]*(-3-ch) <= p[h]*t+s <= p[h]*(3+ch+1)-1
                    collar_checks += 1
    # Explicit K=0,C=1 source count, p=(2,1,1,1,1,1).
    coarse_source_count = prod((2 if h == 0 else 1)*3 for h in range(N))
    assert coarse_source_count == 2*3**6
    return {'status': 'pass', 'new_C6_bound': False,
            'scope': 'staggered extrema/reflection/period/collar/finite-L1 interface example only',
            'extrema_direct_enumeration_cases': extrema_cases,
            'allowed_reflection_cases': reflection_cases, 'period_translation_cases': 36,
            'reflected_endpoint_collar_checks': collar_checks,
            'explicit_coarse_source_mass': str(coarse_source_count),
            'unit_spacing_example': interface}


def main():
    parser = ArgumentParser(description=__doc__)
    parser.add_argument('--selfcheck', action='store_true')
    parser.add_argument('--candidate', type=Path)
    parser.add_argument('--frozen', type=Path)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    assert args.selfcheck != bool(args.candidate), 'choose exactly one mode'
    here = Path(__file__).resolve()
    try:
        if args.selfcheck:
            result = selfcheck()
        else:
            assert args.frozen is not None
            candidate = json.loads(args.candidate.read_text())
            assert candidate.get('spacing', 'frozen_R6_003') == 'frozen_R6_003'
            frozen = json.loads(args.frozen.read_text())
            result = certify(candidate, Q(frozen['lower']['exact_interval'][1]))
            result['candidate_sha256'] = hashlib.sha256(args.candidate.read_bytes()).hexdigest()
            result['frozen_sha256'] = hashlib.sha256(args.frozen.read_bytes()).hexdigest()
            result['new_C6_bound'] = bool(result['above_target'])
            result['candidate_only_status'] = 'primary must integrate theory and independent review'
    except BudgetStop as exc:
        result = {'status': 'budget_stop', 'reason': str(exc), 'new_C6_bound': False}
    result['code_sha256'] = hashlib.sha256(here.read_bytes()).hexdigest()
    helper = here.parent.parent/'correlated_fraction_checker.py'
    result['helper_sha256'] = hashlib.sha256(helper.read_bytes()).hexdigest()
    output = args.output or here.with_name('interface_selfcheck_results.json' if args.selfcheck
                                         else 'candidate_certification_results.json')
    assert output.resolve().parent == here.parent, 'write only in owned audit/lower006/staggered'
    output.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
