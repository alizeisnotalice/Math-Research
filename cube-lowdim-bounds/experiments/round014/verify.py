#!/usr/bin/env python3
"""Independently recompute the saved round014 rational leaves; no search/writes.

Import only the previously reviewed round008 Fraction/log/tail/b2 primitives.
Never import the round014 certifier or run an adaptive subdivision or LP.
"""
from __future__ import annotations
import argparse
from fractions import Fraction as Q
from hashlib import sha256
import importlib.util
import json
from pathlib import Path
import sys

sys.dont_write_bytecode = True
ROOT = Path(__file__).resolve().parent
PASS = 'certified_under_reviewed_angular_and_transfer'


def require(value, label):
    if not value:
        raise ValueError(label)


def find_core(explicit):
    if explicit is not None:
        return explicit.resolve()
    for path in (ROOT.parent/'round008/certify.py',
                 ROOT.parent.parent/'Math-Research-lowdim-agent/cube-lowdim-bounds/experiments/round008/certify.py'):
        if path.is_file():
            return path.resolve()
    raise FileNotFoundError('Pass --core <round008/certify.py> after relocation.')


def verify(path, dependencies_root=None):
    data = json.loads((ROOT/'certificate.json').read_text())
    require(data['round'] == 14 and data['status'] == PASS, 'round and successful status')
    require(sha256(path.read_bytes()).hexdigest() == data['core_sha256'], 'reviewed core hash')
    require(sha256((ROOT/'certify.py').read_bytes()).hexdigest() == data['certifier_sha256'], 'certifier hash')
    dependencies = {'round012/quartic_angular.md', 'round012/transfer_extension.md',
                    'round012/small_radius.json', 'round013/outer_analysis.md', 'round013/outer_audit.md'}
    require(set(data['analytic_dependencies']) == dependencies, 'complete analytic dependencies')
    roots = ([dependencies_root.resolve()] if dependencies_root is not None else
             [ROOT.parent, ROOT.parent.parent/'Math-Research-lowdim-agent/cube-lowdim-bounds/experiments'])
    matching = [root for root in roots if all((root/name).is_file() and
                sha256((root/name).read_bytes()).hexdigest() == digest
                for name, digest in data['analytic_dependencies'].items())]
    require(bool(matching), 'one complete matching analytic dependency root; pass --dependencies-root if needed')
    dependency_root = matching[0]
    require(data['limits'] == {'max_candidate_attempts': 2, 'max_boxes_per_attempt': 100000,
                              'seconds_per_attempt': 285.0, 'tail_terms': 28, 'log_terms': 72}, 'prescribed budget')
    attempts = data['attempts']
    require(1 <= len(attempts) <= 2, 'candidate count')
    intended = [('657/1000', '14/25'), ('16/25', '27/50')]
    for index, attempt in enumerate(attempts):
        e2_text, e4_text = intended[index]
        require(attempt['attempt'] == index+1, 'attempt sequence')
        require(attempt['parameters'] == {'R': '151/80', 'e2': e2_text, 'e3': '0', 'e4': e4_text},
                'preselected candidate sequence')
        require(index == len(attempts)-1 or attempt['status'] != PASS, 'stop after first success')
    require(attempts[-1]['status'] == PASS, 'last attempt is successful')
    saved = attempts[-1]
    R, e2, e4 = (Q(saved['parameters'][key]) for key in ('R', 'e2', 'e4'))
    require(R > 1 and e2 >= 0 and e4 >= 0, 'rational parameter domain')
    spec = importlib.util.spec_from_file_location('round014_independent_saved_leaf_core', path)
    core = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = core
    spec.loader.exec_module(core)
    core.R = R
    I = core.I

    def quartic_diagonal(r):
        s = 1-r
        a = 1-r/2
        return 8*(a*r.square()+Q(22, 3)*(a*a*a)
                  +2*(1+s*s*s)*core.log_i(r/2)-2*s.square()*core.xlogx(s))

    def angular_cap(r):
        if r.hi <= Q(2, 3):
            q = -Q(5, 2)*r-3+18/r-4/r.square()
        elif r.lo >= Q(2, 3):
            q = (165-101*r)/8+1/(2*r.square())
        else:
            raise ValueError('leaf straddles the exact 2/3 curvature split')
        return 12*core.log_i(r/2)+q

    def lower(r):
        diagonal = (2*(2*core.tail_series(r)+I.point(1)/R)/(R-1)
                    -e2*core.b2_diag(r)-e4*quartic_diagonal(r))
        curvature = (e2*64*(1-r).square()/(2-r).square()/(2-r)
                     +I.point(1)/(R**3*(R-1))-4*e4*angular_cap(r))
        # d^2 <= r^2 <= r.hi^2, so a negative lower curvature uses r.hi^2.
        return (diagonal+min(curvature.lo, Q(0))*r.hi*r.hi).lo

    leaves = saved['accepted_leaves']
    require(bool(leaves), 'nonempty complete leaf list')
    require(Q(leaves[0]['lo']) == Q(1, 10) and Q(leaves[-1]['hi']) == 1, 'cover endpoints')
    values = []
    for index, leaf in enumerate(leaves):
        lo, hi, recorded = (Q(leaf[key]) for key in ('lo', 'hi', 'lower'))
        require(Q(1, 10) <= lo < hi <= 1, 'positive leaf width and domain')
        if index:
            require(Q(leaves[index-1]['hi']) == lo, 'exact adjacency without gaps/overlaps')
        value = lower(I(lo, hi))
        require(value >= 0 and value == recorded, 'exactly reproduced nonnegative rational leaf')
        values.append(value)
    inner = saved['inner']
    require(inner['accepted_boxes'] == len(leaves) and inner['splits']+2 == len(leaves), 'binary leaf/split counts')
    require(inner['status'] == 'certified' and inner['cover_verified'], 'coverage metadata')
    require(inner['unresolved_boxes'] == 0 and inner['unresolved_intervals'] == [], 'zero unresolved boxes')
    require(len(leaves)+inner['splits'] <= 100000 and 0 <= inner['elapsed_seconds'] <= 285, 'recorded finite budget')
    minimum = min(values)
    require(minimum == Q(inner['minimum_accepted_lower']), 'exact saved minimum')
    text = '\n'.join(f"{leaf['lo']}/{leaf['hi']}" for leaf in leaves)
    require(sha256(text.encode()).hexdigest() == inner['accepted_leaf_sha256'], 'complete leaf digest')

    h4_small = 8*(Q(553, 150)+Q(729, 1000)*core.log_point(Q(1, 9))+core.log_point(Q(1, 11)))
    t2 = 3-5*core.log_point(2)
    small = saved['small_radius']
    supplied = json.loads((dependency_root/'round012/small_radius.json').read_text())
    require(h4_small.hi < 0 and t2.hi < 0, 'strict small-radius signs')
    require(h4_small.pair() == small['h4_upper_interval'] == supplied['h4_uniform_upper_interval'], 'quartic small interval')
    require(t2.pair() == small['quadratic_T_half_interval'], 'quadratic small interval')
    require(small['h4_strictly_negative'] and small['matches_round012_interval'] and small['quadratic_T_half_negative'],
            'small-radius metadata')

    log2 = core.log_point(2)
    h2 = 8-8*log2
    h4 = Q(136, 3)-64*log2
    d1 = 2*(2*core.tail_series(I.point(1))+I.point(1)/R)/(R-1)
    exterior = d1-e2*h2-e4*h4
    dominance = R-(e2+e4)*(R-1)/4
    cost = R*R/(R-1)-e2/6-e4/15
    previous = Q(173331, 44750)
    checks = saved['outer_and_cost']
    require(h2.pair() == checks['H2_interval'] and h4.pair() == checks['H4_interval'], 'endpoint intervals')
    require(d1.pair() == checks['D_R1_interval'], 'tail lower interval')
    require(exterior.pair() == checks['new_exterior_slack_interval'] and exterior.lo > 0, 'strict new exterior condition')
    require(dominance == Q(checks['kernel_dominance_margin']) and dominance > 0, 'kernel dominance')
    require(previous == Q(checks['previous_U12']) and cost == Q(checks['cost']) and cost < previous, 'strict exact cost')
    require(previous-cost == Q(checks['U12_minus_cost']) and checks['strict'], 'exact improvement and metadata')
    if saved['attempt'] == 1:
        require(cost == Q(411881, 106500) and previous-cost == Q(112307, 19063500), 'intended final first-candidate bound')
    return {'status': 'all_saved_rational_leaves_independently_verified',
            'method': 'no adaptive search or LP; independent formula; reviewed Fraction primitives only',
            'leaf_count': len(leaves), 'exact_cover': ['1/10', '1'],
            'all_recorded_lowers_exactly_reproduced': True,
            'minimum_lower_display': float(minimum), 'unresolved_boxes': 0,
            'small_h4_upper_display': float(h4_small.hi),
            'new_exterior_slack_lower_display': float(exterior.lo),
            'kernel_dominance_margin': str(dominance), 'cost': str(cost),
            'improvement': str(previous-cost), 'source_and_analytic_hashes_match': True,
            'scope': 'scalar certificate verified; angular, small-radius and general-input transfer remain explicit reviewed analytic dependencies'}


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--core', type=Path)
    parser.add_argument('--dependencies-root', type=Path)
    args = parser.parse_args()
    print(json.dumps(verify(find_core(args.core), args.dependencies_root), indent=2))
