#!/usr/bin/env python3
"""Independent Fraction checker for *listed* two-axis parity witness blocks.

No author modules, floating point geometry, algebraic-number package or MATLAB.
The selfcheck is an interface example, not a new C6-bound experiment.
"""
from argparse import ArgumentParser
from fractions import Fraction as Q
from functools import lru_cache
from itertools import product
from math import prod
from pathlib import Path
import hashlib
import json
import time

N = 6


class BudgetStop(Exception):
    pass


class Budget:
    def __init__(self, nodes=10**6, seconds=20*60):
        self.node_cap, self.seconds = nodes, seconds
        self.nodes, self.started = 0, time.monotonic()

    def check(self, node=False):
        self.nodes += int(node)
        if self.nodes > self.node_cap:
            raise BudgetStop('union node cap')
        if time.monotonic() - self.started > self.seconds:
            raise BudgetStop('wall time cap')


def floorq(x):
    return x.numerator // x.denominator


def ceilq(x):
    return -floorq(-x)


@lru_cache(None)
def root6(q, denominator):
    """l^6 <= q < h^6, h-l=1/denominator, proved by integer powers."""
    q = Q(q)
    assert q >= 0 and denominator > 0
    z = q.numerator * denominator**6 // q.denominator
    lo, hi = 0, 1
    while hi**6 <= z:
        hi *= 2
    while hi-lo > 1:
        mid = (lo+hi)//2
        if mid**6 <= z:
            lo = mid
        else:
            hi = mid
    a, b = Q(lo, denominator), Q(lo+1, denominator)
    assert a**6 <= q < b**6
    return a, b


def scale(interval, q):
    a, b = interval[0]*q, interval[1]*q
    return min(a, b), max(a, b)


def default_spacing(denominator):
    """Enclose the frozen R6-003 spacing, using only positive operations."""
    u, v = root6(Q(2), denominator), root6(Q(3), denominator)
    return (
        (1+u[0]**4/2, 1+u[1]**4/2),
        ((u[0]**3+u[0]**4*v[0]**3)/4,
         (u[1]**3+u[1]**4*v[1]**3)/4),
        *((((1+u[0])/2, (1+u[1])/2),)*4),
    )


def parity_sum(start, length):
    # Avoid Python's floating return value for (-1)**negative_integer.
    return (1-2*(start % 2)) * (length % 2)


def block_mass(start, lengths, axes, delta):
    i, j = axes
    return (Q(prod(lengths)) + delta * parity_sum(start[i], lengths[i])
            * parity_sum(start[j], lengths[j])
            * prod(lengths[k] for k in range(N) if k not in axes))


def prune(boxes, budget):
    boxes = tuple(sorted(set(boxes)))
    out = []
    for k, box in enumerate(boxes):
        budget.check()
        if not any(other != box and all(
                other[2*j] <= box[2*j] and box[2*j+1] <= other[2*j+1]
                for j in range(N)) for other in boxes):
            out.append(box)
    return tuple(out)


def union_integer(boxes, budget):
    """Exact integer-coordinate slicing; cache changes only running time."""
    @lru_cache(None)
    def rec(active):
        budget.check(node=True)
        if not active:
            return 0
        if len(active[0]) == 2:
            spans = sorted(active)
            a, b = spans[0]
            total = 0
            for c, d in spans[1:]:
                if c > b:
                    total += b-a
                    a, b = c, d
                else:
                    b = max(b, d)
            return total+b-a
        cuts = sorted({x for b in active for x in b[:2]})
        total = 0
        for a, b in zip(cuts, cuts[1:]):
            child = tuple(sorted({x[2:] for x in active
                                  if x[0] <= a and b <= x[1]}))
            total += (b-a)*rec(child)
        return total
    return rec(tuple(sorted(set(boxes))))


def inward_outer(candidate, denominator, grid, budget):
    axes = tuple(candidate['axes'])
    assert len(axes) == 2 and len(set(axes)) == 2
    assert all(type(x) is int and 0 <= x < N for x in axes)
    delta, level = Q(candidate['delta']), Q(candidate.get('level', '1'))
    assert -1 <= delta <= 1 and level > 0
    if candidate.get('spacing') == 'unit_example':
        sp = ((Q(1), Q(1)),)*N
    else:
        assert candidate.get('spacing', 'frozen_R6_003') == 'frozen_R6_003'
        sp = default_spacing(denominator)
    periods = tuple(2 if k in axes else 1 for k in range(N))
    inner, outer, labels = [], [], []
    margin = None
    collar = [0]*N
    for tag, block in enumerate(candidate['blocks']):
        budget.check()
        start, lengths = tuple(block['start']), tuple(block['lengths'])
        assert len(start) == len(lengths) == N
        assert all(type(k) is int for k in start)
        assert all(type(k) is int and 1 <= k <= 4 for k in lengths)
        mass = block_mass(start, lengths, axes, delta)
        assert mass >= 0
        if not mass:
            continue
        r0, r1 = scale(root6(mass/level, denominator), Q(1, 2))
        assert (2*r0)**6 * level <= mass < (2*r1)**6 * level
        ib, ob = [], []
        for k, (q, ell, a, p) in enumerate(zip(start, lengths, sp, periods)):
            first, last = scale(a, q), scale(a, q+ell-1)
            half = scale(a, Q(p, 2))
            ilo, ihi = max(Q(0), last[1]-r0), min(half[0], first[0]+r0)
            olo, ohi = max(Q(0), last[0]-r1), min(half[1], first[1]+r1)
            ib.extend((ceilq(ilo*grid)+1, floorq(ihi*grid)-1))
            ob.extend((floorq(olo*grid), ceilq(ohi*grid)))
        if all(ob[2*k] < ob[2*k+1] for k in range(N)):
            outer.append(tuple(ob))
        if not all(ib[2*k] < ib[2*k+1] for k in range(N)):
            continue
        local_margin = None
        for k, (q, ell, a, p) in enumerate(zip(start, lengths, sp, periods)):
            lo, hi = Q(ib[2*k], grid), Q(ib[2*k+1], grid)
            assert 0 < lo < hi < p*a[0]/2
            for source_index in (q, q+ell-1):
                al, ah = scale(a, source_index)
                gap = min(r0+lo-ah, r0-hi+al)
                assert gap > 0
                local_margin = gap if local_margin is None else min(local_margin, gap)
            # Reflections send all source indices into [-S_k,S_k].
            s = max(abs(q), abs(q+ell-1))
            collar[k] = max(collar[k], (s+p-1)//p)
        margin = local_margin if margin is None else min(margin, local_margin)
        inner.append(tuple(ib))
        labels.append({'block': tag, 'start': start, 'lengths': lengths,
                       'mass': str(mass), 'rational_radius': str(r0),
                       'cube_boundary_margin': str(local_margin)})
    inner, outer = prune(inner, budget), prune(outer, budget)
    return {'inner': inner, 'outer': outer, 'labels': labels, 'margin': margin,
            'collar': tuple(collar), 'level': level, 'axes': axes,
            'delta': delta, 'periods': periods, 'spacing': sp}


def finite_transfer(raw_ratio, level, target, collar, margin, period_mass=4, budget=None):
    """An explicit finite L1 input, conditional only on raw_ratio > target."""
    assert raw_ratio > target >= 0 and level > 0 and margin > 0
    beta = level*(raw_ratio+target)/(2*raw_ratio)
    assert 0 < beta < level
    k = 1
    while True:
        if budget is not None:
            budget.check()
        loss = prod((Q(2*k+1, 2*k+2*c+1) for c in collar), start=Q(1))
        value = (beta/level)*raw_ratio*loss
        if value > target:
            break
        k *= 2
        if k.bit_length() > 4096:
            raise BudgetStop('finite K bit-length cap')
    source_mass = period_mass*prod(2*k+2*c+1 for c in collar)
    return {'receiver_coarse_index_radius_K': k, 'coarse_collar': collar,
            'source_mass_exact': str(source_mass), 'period_mass': period_mass,
            'threshold_beta': str(beta),
            'bump_halfwidth_eta': str(margin/2),
            'boundary_loss_exact': str(loss), 'finite_ratio': str(value),
            'finite_ratio_minus_target': str(value-target),
            'bump': '(2 eta)^(-6) 1_[-eta,eta]^6',
            'source_coarse_cells': 't_j in [-K-C_j,K+C_j], residues 0..p_j-1',
            'strictly_above_target': True}


def certify(candidate, target, denominator=10**50, grid=10**35,
            node_cap=10**6, seconds=20*60):
    budget = Budget(node_cap, seconds)
    geo = inward_outer(candidate, denominator, grid, budget)
    vi = Q(union_integer(geo['inner'], budget), grid**N)
    vo = Q(union_integer(geo['outer'], budget), grid**N)
    raw_lo, raw_hi = 16*geo['level']*vi, 16*geo['level']*vo
    assert 0 <= raw_lo <= raw_hi
    out = {'status': 'certified_listed_union', 'scope': 'listed blocks only; no claim of complete maximal set',
           'root_denominator': str(denominator), 'grid_denominator': str(grid),
           'input_block_count': len(candidate['blocks']),
           'inner_retained_boxes': len(geo['inner']), 'outer_retained_boxes': len(geo['outer']),
           'listed_union_ratio_interval': [str(raw_lo), str(raw_hi)],
           'target': str(target), 'ratio_lower_minus_target': str(raw_lo-target),
           'geometric_cube_margin': str(geo['margin']),
           'coarse_collar': geo['collar'], 'union_nodes': budget.nodes,
           'node_cap': node_cap, 'wall_seconds_cap': seconds,
           'elapsed_seconds': time.monotonic()-budget.started,
           'above_target': raw_lo > target,
           'outer_proves_listed_candidate_below_target': raw_hi < target}
    if raw_lo > target:
        out['finite_L1_witness'] = finite_transfer(raw_lo, geo['level'], target,
                                                  geo['collar'], geo['margin'], budget=budget)
    return out


def selfcheck():
    root_cases = (Q(0), Q(1), Q(1, 64), Q(729, 64), Q(17, 5), Q(4096))
    for x in root_cases:
        a, b = root6(x, 10**24)
        assert a**6 <= x < b**6
    cases = 0
    for axes in ((0, 1), (2, 5)):
        for delta in (Q(-1), Q(-1, 3), Q(0), Q(1, 3), Q(1)):
            for seed in range(12):
                q = tuple((seed+2*j) % 5-2 for j in range(N))
                ell = tuple((seed+j) % 3+1 for j in range(N))
                exact = sum((1+delta*(1-2*((k[axes[0]]+k[axes[1]]) % 2))
                             for k in product(*(range(q[j], q[j]+ell[j]) for j in range(N)))), Q(0))
                assert exact == block_mass(q, ell, axes, delta)
                cases += 1
            residues = product(*(range(2 if j in axes else 1) for j in range(N)))
            assert sum((1+delta*(1-2*((k[axes[0]]+k[axes[1]]) % 2))
                        for k in residues), Q(0)) == 4
    # Explicit coarse-source mass check; K=0,C=1 gives 2916 atoms.
    axes, delta, k, c = (0, 1), Q(1, 3), 0, 1
    p = (2, 2, 1, 1, 1, 1)
    ranges = [range(pj*(-k-c), pj*(k+c+1)) for pj in p]
    coarse_mass = sum((1+delta*(1-2*((x[0]+x[1]) % 2)) for x in product(*ranges)), Q(0))
    assert coarse_mass == 4*(2*k+2*c+1)**6
    # Independent known-volume oracle for exact slicing in dimension six.
    a = (0, 2)*6
    b = (1, 3)*6
    ub = Budget()
    assert union_integer((a, b), ub) == 2*2**6-1
    assert prune((a, (0, 1)*6, b), ub) == tuple(sorted((a, b)))
    example = {'axes': [0, 1], 'delta': '1/3', 'level': '1', 'spacing': 'unit_example',
               'blocks': [{'start': [0]*6, 'lengths': [1]*6},
                          {'start': [1, 1, 0, 0, 0, 0], 'lengths': [1]*6}]}
    interface = certify(example, Q(1, 10), denominator=10**24, grid=10**12)
    assert interface['above_target'] and interface['finite_L1_witness']['strictly_above_target']
    # Check the stated collar on every sign reflection of every source endpoint.
    for block in example['blocks']:
        for signs in product((-1, 1), repeat=N):
            for j in range(N):
                s = signs[j]*block['start'][j]
                cj = interface['coarse_collar'][j]
                for t in (-3, 3):
                    assert p[j]*(-3-cj) <= p[j]*t+s <= p[j]*(3+cj+1)-1
    return {'status': 'pass', 'new_C6_bound': False,
            'scope': 'mass/parity/root/union/collar/finite-L1 interface example only',
            'root_cases': len(root_cases), 'direct_mass_cases': cases,
            'period_mass_cases': 10, 'coarse_mass': str(coarse_mass),
            'union_known_volume': '127', 'unit_spacing_example': interface}


def main():
    parser = ArgumentParser(description=__doc__)
    parser.add_argument('--candidate', type=Path)
    parser.add_argument('--frozen', type=Path)
    parser.add_argument('--output', type=Path)
    parser.add_argument('--selfcheck', action='store_true')
    args = parser.parse_args()
    here = Path(__file__).resolve()
    assert args.selfcheck != bool(args.candidate), 'choose exactly one mode'
    frozen = None
    if args.frozen:
        frozen = json.loads(args.frozen.read_text())
    try:
        if args.selfcheck:
            result = selfcheck()
            if frozen:
                result['frozen_threshold_readonly'] = {
                    'id': frozen['lower']['id'],
                    'exact_upper': frozen['lower']['exact_interval'][1],
                    'sha256': hashlib.sha256(args.frozen.read_bytes()).hexdigest(),
                    'example_is_new_C6_candidate': False}
        else:
            assert frozen is not None, 'real candidate requires --frozen exact threshold file'
            candidate = json.loads(args.candidate.read_text())
            assert candidate.get('spacing', 'frozen_R6_003') == 'frozen_R6_003'
            target = Q(frozen['lower']['exact_interval'][1])
            result = certify(candidate, target)
            result['candidate_sha256'] = hashlib.sha256(args.candidate.read_bytes()).hexdigest()
            result['frozen_sha256'] = hashlib.sha256(args.frozen.read_bytes()).hexdigest()
            result['new_C6_bound'] = bool(result['above_target'])
            result['candidate_only_status'] = 'primary must integrate theory and independent review'
    except BudgetStop as exc:
        result = {'status': 'budget_stop', 'reason': str(exc), 'new_C6_bound': False}
    result['code_sha256'] = hashlib.sha256(here.read_bytes()).hexdigest()
    output = args.output or here.with_name('interface_selfcheck_results.json' if args.selfcheck
                                         else 'candidate_certification_results.json')
    assert output.resolve().parent == here.parent, 'write only in owned audit/lower006'
    output.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
