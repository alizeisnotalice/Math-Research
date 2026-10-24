#!/usr/bin/env python3
"""One authorized finite rational 2D certification attempt; no LP or grid."""
from __future__ import annotations
from fractions import Fraction as Q
from functools import lru_cache
from hashlib import sha256
import importlib.util
import json
import math
from pathlib import Path
import sys
from time import monotonic

sys.dont_write_bytecode = True
ROOT = Path(__file__).resolve().parent
ARCHIVE = ROOT.parent.parent/'Math-Research-lowdim-agent/cube-lowdim-bounds/experiments'
R = Q(15, 8)
E = {2: Q(69, 100), 6: Q(19, 50), 7: Q(11, 10), 8: Q(3, 4)}
OLD_UPPER = Q(411881, 106500)
DELTA = Q(1, 100)
BOX_LIMIT = 30000
SECONDS = 285.
WIDTH_FLOOR = Q(1, 2**60)


def load_core():
    for source in [ROOT.parent/'round008/certify.py', ARCHIVE/'round008/certify.py']:
        if source.exists():
            spec = importlib.util.spec_from_file_location('round016_fraction_core', source)
            c = importlib.util.module_from_spec(spec)
            sys.modules[spec.name] = c
            spec.loader.exec_module(c)
            c.R = R
            return c, source
    raise FileNotFoundError('round008/certify.py not found')


def make_lower(c):
    # W_n lower enclosures converge at positive arguments by natural interval
    # multiplication, and at zero by a proved monotone envelope. Do not use
    # core.xlogx, whose fixed critical-point floor need not converge.
    def w_lower(lo, hi, n):
        if hi == 0:
            return Q(0)
        if lo == 0:
            if hi <= Q(1, 4):
                return hi**n*c.log_point(hi).lo
            return -Q(1, 2*n)
        ll, lh = c.log_point(lo).lo, c.log_point(hi).hi
        return min(lo**n*ll, lo**n*lh, hi**n*ll, hi**n*lh)

    def negative_log_upper(coefficient_lo, argument_hi):
        assert coefficient_lo >= 0 and 0 <= argument_hi <= 1
        if coefficient_lo == 0 or argument_hi == 0:
            return Q(0)  # A safe upper bound, never 0 times an infinity.
        return coefficient_lo*min(Q(0), c.log_point(argument_hi).hi)

    @lru_cache(maxsize=None)
    def f_lower(x):
        return sum((x**n/((n+2)*R**(n+1)) for n in range(29)), Q(0))

    @lru_cache(maxsize=None)
    def even_pair_lower(x):
        return 2*sum((x**(2*k)/((2*k+2)*R**(2*k+1)) for k in range(15)), Q(0))

    def lower(box):
        rl, rh, ql, qh = box
        assert DELTA <= rl <= rh <= 1 and 0 <= ql <= qh <= 1
        dl, dh = rl*ql, rh*qh
        ul, uh = rl*(1+ql)/2, rh*(1+qh)/2
        vl, vh = rl*(1-qh)/2, rh*(1-ql)/2
        sl, sh = 1-rh, 1-rl
        ah, bh = 1-vl, 1-ul
        zpl, zph, zml, zmh = 1+dl, 1+dh, 1-dh, 1-dl
        tail = 2*(2*f_lower(rl)+even_pair_lower(dl))/(R-1)
        upper = {}
        for p in E:
            n = p-1
            poly = sum(((ah**(n-k)*(zph**k+sh**k)
                         +bh**(n-k)*(zmh**k+sh**k))/Q(n-k) for k in range(n)), Q(0))
            cu, cv = zpl**n+sl**n, zml**n+sl**n
            bp = (poly+negative_log_upper(cu, uh)+negative_log_upper(cv, vh)
                  -w_lower(zpl, zph, n)-w_lower(zml, zmh, n)-2*w_lower(sl, sh, n))
            upper[p] = 2*p*bp
        total = tail-sum((E[p]*upper[p] for p in E), Q(0))
        return total, tail, upper
    return lower


def auxiliary(c):
    harm7 = sum((Q(1, k) for k in range(1, 8)), Q(0))
    B = (2*Q(101, 100)**6*harm7+c.log_point(Q(1, 101))
         +Q(99, 100)**7*c.log_point(Q(1, 99)))
    expected_B_hi = Q(-237533260236026023244253198287469777005146586609,
                       70000000000000000000000000000000000000000000000)
    endpoints = {}
    for p in E:
        a = Q(2, p-1)
        for k in range(p+1):
            if k != 1:
                a += 2*math.comb(p, k)*Q(2)**(p-k)*(-1)**k*(Q(2)**(k-1)-1)/Q(k-1)
        endpoints[p] = a-p*2**p*c.log_point(Q(2))
    D1 = 2*(2*c.tail_series(c.I.point(1))+c.I.point(1)/R)/(R-1)
    outer = D1-sum((E[p]*endpoints[p] for p in E), c.I.point(0))
    dominance = R-(R-1)*sum(E.values())/4
    cost = R*R/(R-1)-sum((2*e/Q((p+1)*(p+2)) for p, e in E.items()), Q(0))
    return {'small_radius': {'range': '0<r<=1/100, open quadrants', 'B_interval': B.pair(),
                             'B_upper_negative': B.hi < 0,
                             'B_upper_matches_audit': B.hi == expected_B_hi},
            'exterior': {'D_R1_interval': D1.pair(), 'Hp_intervals': {str(p): h.pair() for p, h in endpoints.items()},
                         'new_slack_interval': outer.pair(), 'new_slack_positive': outer.lo > 0},
            'kernel_dominance': str(dominance), 'kernel_dominance_positive': dominance > 0,
            'cost': str(cost), 'cost_display': float(cost), 'previous_strict_upper': str(OLD_UPPER),
            'improvement': str(OLD_UPPER-cost), 'strict_cost_improvement': cost < OLD_UPPER}


def box_row(box):
    return {'r_lo': str(box[0]), 'r_hi': str(box[1]), 'q_lo': str(box[2]), 'q_hi': str(box[3])}


def main():
    output = ROOT/'certificate.json'
    if output.exists():
        raise SystemExit('Certificate attempt file already exists; no automatic retry permitted.')
    output.write_text(json.dumps({'round': 16, 'status': 'started', 'candidate_attempts': 1})+'\n')
    c, core_source = load_core()
    lower = make_lower(c)
    aux = auxiliary(c)
    started = monotonic()
    root = (DELTA, Q(1), Q(0), Q(1))
    nodes = [{'id': 0, 'parent': None, **box_row(root), 'status': 'pending'}]
    todo = [(0, root)]
    leaves, unresolved = [], []
    evaluated, splits = 0, 0
    stop = 'all_boxes_closed'
    minimum = None
    worst_evaluated = None
    while todo:
        # Hard limits are checked inside the loop before every evaluated box.
        if evaluated >= BOX_LIMIT:
            stop = 'box_evaluation_budget_exhausted'
            break
        if monotonic()-started >= SECONDS:
            stop = 'time_budget_exhausted'
            break
        node_id, box = todo.pop()
        evaluated += 1
        value, tail, upper = lower(box)
        if worst_evaluated is None or value < worst_evaluated[0]:
            worst_evaluated = (value, box)
        if value >= 0:
            nodes[node_id]['status'] = 'accepted'
            leaves.append({'node_id': node_id, **box_row(box), 'tail_lower': str(tail),
                           'bump_upper_by_power': {str(p): str(x) for p, x in upper.items()},
                           'generator_lower': str(value)})
            minimum = value if minimum is None else min(minimum, value)
        else:
            rl, rh, ql, qh = box
            if max(rh-rl, qh-ql) <= WIDTH_FLOOR:
                nodes[node_id]['status'] = 'unresolved_width_floor'
                unresolved.append({'node_id': node_id, **box_row(box), 'last_generator_lower': str(value),
                                   'reason': 'width_floor'})
            else:
                axis = 'r' if (rh-rl)/(1-DELTA) >= qh-ql else 'q'
                if axis == 'r':
                    mid = (rl+rh)/2
                    children = [(rl, mid, ql, qh), (mid, rh, ql, qh)]
                else:
                    mid = (ql+qh)/2
                    children = [(rl, rh, ql, mid), (rl, rh, mid, qh)]
                child_ids = []
                for child in children:
                    child_id = len(nodes)
                    child_ids.append(child_id)
                    nodes.append({'id': child_id, 'parent': node_id, **box_row(child), 'status': 'pending'})
                    todo.append((child_id, child))
                nodes[node_id].update({'status': 'split', 'axis': axis, 'midpoint': str(mid), 'children': child_ids})
                splits += 1
        if evaluated % 500 == 0:
            print(json.dumps({'evaluated_boxes': evaluated, 'accepted_leaves': len(leaves),
                              'pending_boxes': len(todo), 'elapsed_seconds': monotonic()-started}), flush=True)
    numeric_elapsed = monotonic()-started
    for node_id, box in todo:
        unresolved.append({'node_id': node_id, **box_row(box), 'last_generator_lower': None, 'reason': stop})
    cover = not unresolved and all(n['status'] in ['accepted', 'split'] for n in nodes)
    aux_passed = (aux['small_radius']['B_upper_negative'] and aux['small_radius']['B_upper_matches_audit']
                  and aux['exterior']['new_slack_positive'] and aux['kernel_dominance_positive']
                  and aux['strict_cost_improvement'])
    success = cover and aux_passed
    archive = ARCHIVE if (ARCHIVE/'round013/outer_analysis.md').exists() else ROOT.parent
    dependencies = {'round013/outer_analysis.md': archive/'round013/outer_analysis.md',
                    'round016/continuous_interface.md': ROOT/'continuous_interface.md',
                    'round016/certification_plan.md': ROOT/'certification_plan.md'}
    result = {'round': 16, 'status': 'closed_pending_independent_review' if success else 'not_closed',
              'candidate_attempts': 1, 'parameters': {'R': str(R), 'e_by_power': {str(p): str(e) for p, e in E.items()}},
              'root_domain': box_row(root), 'limits': {'actual_box_evaluations': BOX_LIMIT, 'seconds': SECONDS,
                                                     'width_floor': str(WIDTH_FLOOR), 'tail_n_max': 28, 'even_pair_k_max': 14},
              'inner': {'actual_evaluated_boxes': evaluated, 'accepted_leaves': len(leaves), 'splits': splits,
                        'unresolved_boxes': len(unresolved), 'stop_reason': stop, 'complete_tree_cover': cover,
                        'minimum_accepted_lower': str(minimum) if minimum is not None else None,
                        'minimum_accepted_lower_display': float(minimum) if minimum is not None else None,
                        'elapsed_seconds': numeric_elapsed, 'time_checked_before_every_box': True,
                        'worst_evaluated_box': {'lower': str(worst_evaluated[0]), **box_row(worst_evaluated[1])} if worst_evaluated else None},
              'auxiliary': aux, 'accepted_leaf_data': leaves, 'unresolved_queue': unresolved, 'partition_tree': nodes,
              'accepted_leaf_sha256': sha256(json.dumps(leaves, sort_keys=True, separators=(',', ':')).encode()).hexdigest(),
              'partition_tree_sha256': sha256(json.dumps(nodes, sort_keys=True, separators=(',', ':')).encode()).hexdigest(),
              'core_dependency': 'round008/certify.py', 'core_sha256': sha256(core_source.read_bytes()).hexdigest(),
              'certifier_sha256': sha256(Path(__file__).read_bytes()).hexdigest(),
              'analytic_dependencies': {name: sha256(path.read_bytes()).hexdigest() for name, path in dependencies.items()},
              'method': 'Exact Fraction one-sided weighted logarithm envelopes, positive full-angular tail truncations, finite 2D binary partition',
              'scope': 'One authorized attempt. Only a fully closed partition plus independent review and transfer audit can update the rigorous upper bound.'}
    output.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({k: v for k, v in result.items() if k not in ['accepted_leaf_data', 'unresolved_queue', 'partition_tree']}, indent=2), flush=True)


if __name__ == '__main__':
    main()
