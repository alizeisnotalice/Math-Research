#!/usr/bin/env python3
"""Round14 bounded rational certification; run only after approved candidates.

Imports round008 interval/log/tail/quadratic primitives. The quartic formula,
angular correction and new round013 exterior row are explicit here.
"""
from __future__ import annotations
import argparse
from fractions import Fraction as Q
from hashlib import sha256
import importlib.util
import json
from pathlib import Path
import sys
from time import monotonic

sys.dont_write_bytecode = True
ROOT = Path(__file__).resolve().parent
R = Q(151, 80)
INCUMBENT = Q(173331, 44750)
CANDIDATES = [(Q(657, 1000), Q(14, 25)), (Q(16, 25), Q(27, 50))]
MAX_BOXES = 100000
SECONDS = 285.0
LOCAL_ARCHIVE = ROOT.parent.parent / 'Math-Research-lowdim-agent/cube-lowdim-bounds/experiments'


def find_core(explicit=None):
    if explicit is not None:
        return explicit.resolve()
    for path in [ROOT.parent/'round008/certify.py', LOCAL_ARCHIVE/'round008/certify.py']:
        if path.exists():
            return path.resolve()
    raise FileNotFoundError('Pass --core <round008/certify.py> after relocation.')


def proof_root():
    # During local continuation use the archived authoritative proofs. After
    # archiving this round, the same proofs live beside round014.
    return LOCAL_ARCHIVE if (LOCAL_ARCHIVE/'round013/outer_analysis.md').exists() else ROOT.parent


def load_core(path):
    spec = importlib.util.spec_from_file_location('round014_fraction_core', path)
    c = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = c
    spec.loader.exec_module(c)
    c.R = R
    return c


def formulas(c, e2, e4):
    I = c.I
    def cube(x):
        return x*x*x
    def b4_diag(r):
        s = 1-r
        a = 1-r/2
        return 8*(a*r.square()+Q(22, 3)*cube(a)
                  +2*(1+cube(s))*c.log_i(r/2)-2*s.square()*c.xlogx(s))
    def c4(r):
        if r.hi <= Q(2, 3):
            q = -Q(5, 2)*r-3+18/r-4/r.square()
        elif r.lo >= Q(2, 3):
            q = (165-101*r)/8+1/(2*r.square())
        else:
            raise ValueError('C4 interval straddles 2/3')
        return 12*c.log_i(r/2)+q
    def lower(r):
        tail = 2*(2*c.tail_series(r)+I.point(1)/R)/(R-1)
        diag = tail-e2*c.b2_diag(r)-e4*b4_diag(r)
        k = (e2*64*(1-r).square()/(2-r).square()/(2-r)
             +I.point(1)/(R**3*(R-1))-4*e4*c4(r))
        return diag+min(k.lo, Q(0))*r.hi*r.hi
    return lower


def auxiliary(c, e2, e4):
    # Independent exact checks for the omitted small-radius strip.
    h4 = 8*(Q(553, 150)+Q(729, 1000)*c.log_point(Q(1, 9))+c.log_point(Q(1, 11)))
    t2 = 3-5*c.log_point(Q(2))
    supplied = json.loads((proof_root()/'round012/small_radius.json').read_text())
    small = {'range': '0<r<=1/10', 'h4_upper_interval': h4.pair(),
             'h4_strictly_negative': h4.hi < 0,
             'matches_round012_interval': h4.pair() == supplied['h4_uniform_upper_interval'],
             'quadratic_T_half_interval': t2.pair(), 'quadratic_T_half_negative': t2.hi < 0,
             'analytic_dependency': 'round012/transfer_extension.md; t=u/(1-v)<=r; quadratic T monotonicity'}
    log2 = c.log_point(Q(2))
    h2 = 8-8*log2
    h4end = Q(136, 3)-64*log2
    d1 = 2*(2*c.tail_series(c.I.point(1))+c.I.point(1)/R)/(R-1)
    exterior = d1-e2*h2-e4*h4end
    dominance = R-(e2+e4)*(R-1)/4
    cost = R*R/(R-1)-e2/6-e4/15
    outer = {'range': '1<=r<R; r>=R follows from K>=0 and K=0',
             'H2_interval': h2.pair(), 'H4_interval': h4end.pair(),
             'D_R1_interval': d1.pair(), 'new_exterior_slack_interval': exterior.pair(),
             'new_exterior_slack_display': [float(exterior.lo), float(exterior.hi)],
             'kernel_dominance_margin': str(dominance),
             'cost': str(cost), 'cost_display': float(cost),
             'previous_U12': str(INCUMBENT), 'U12_minus_cost': str(INCUMBENT-cost),
             'U12_minus_cost_display': float(INCUMBENT-cost),
             'strict': e2 >= 0 and e4 >= 0 and exterior.lo > 0 and dominance > 0 and cost < INCUMBENT}
    return small, outer


def certify(c, e2, e4):
    lower = formulas(c, e2, e4)
    start = monotonic()
    todo = [c.I(Q(1, 10), Q(2, 3)), c.I(Q(2, 3), Q(1))]
    leaves = []
    unresolved = []
    splits = 0
    minimum = None
    while todo:
        if len(leaves)+splits >= MAX_BOXES or monotonic()-start > SECONDS:
            unresolved.extend(todo)
            break
        box = todo.pop()
        bound = lower(box)
        if bound.lo >= 0:
            leaves.append({'lo': str(box.lo), 'hi': str(box.hi), 'lower': str(bound.lo)})
            minimum = bound.lo if minimum is None else min(minimum, bound.lo)
        elif box.hi-box.lo <= Q(1, 2**80):
            unresolved.append(box)
        else:
            mid = (box.lo+box.hi)/2
            todo.extend([c.I(box.lo, mid), c.I(mid, box.hi)])
            splits += 1
        if (len(leaves)+splits) % 500 == 0:
            print(f'accepted={len(leaves)} splits={splits} pending={len(todo)}', flush=True)
    leaves.sort(key=lambda x: Q(x['lo']))
    cover = bool(leaves) and not unresolved and Q(leaves[0]['lo']) == Q(1, 10) and Q(leaves[-1]['hi']) == 1
    cover = cover and all(Q(a['hi']) == Q(b['lo']) for a, b in zip(leaves, leaves[1:]))
    leaf_text = '\n'.join(f"{x['lo']}/{x['hi']}" for x in leaves)
    inner = {'status': 'certified' if cover else 'unresolved', 'accepted_boxes': len(leaves),
             'splits': splits, 'unresolved_boxes': len(unresolved), 'cover_verified': cover,
             'minimum_accepted_lower': str(minimum) if minimum is not None else None,
             'minimum_accepted_lower_display': float(minimum) if minimum is not None else None,
             'accepted_leaf_sha256': sha256(leaf_text.encode()).hexdigest(),
             'unresolved_intervals': [[str(x.lo), str(x.hi)] for x in unresolved],
             'elapsed_seconds': monotonic()-start}
    small, outer = auxiliary(c, e2, e4)
    passed = (cover and small['h4_strictly_negative'] and small['matches_round012_interval']
              and small['quadratic_T_half_negative'] and outer['strict'])
    return {'status': 'certified_under_reviewed_angular_and_transfer' if passed else 'not_closed',
            'parameters': {'R': str(R), 'e2': str(e2), 'e3': '0', 'e4': str(e4)},
            'inner': inner, 'small_radius': small, 'outer_and_cost': outer, 'accepted_leaves': leaves}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--core', type=Path)
    args = ap.parse_args()
    path = ROOT/'certificate.json'
    data = json.loads(path.read_text()) if path.exists() else {'round': 14, 'attempts': []}
    if any(a.get('status') == 'certified_under_reviewed_angular_and_transfer' for a in data['attempts']):
        raise SystemExit('A successful certificate exists; success-stop rule prevents another attempt.')
    if len(data['attempts']) >= len(CANDIDATES):
        raise SystemExit('Two-candidate budget exhausted.')
    number = len(data['attempts'])+1
    e2, e4 = CANDIDATES[number-1]
    data['attempts'].append({'attempt': number, 'status': 'started',
                             'parameters': {'R': str(R), 'e2': str(e2), 'e3': '0', 'e4': str(e4)}})
    path.write_text(json.dumps(data, indent=2)+'\n')
    source = find_core(args.core)
    c = load_core(source)
    result = certify(c, e2, e4)
    result['attempt'] = number
    data['attempts'][-1] = result
    archive = proof_root()
    dependencies = [archive/'round012/quartic_angular.md', archive/'round012/transfer_extension.md',
                    archive/'round012/small_radius.json', archive/'round013/outer_analysis.md',
                    archive/'round013/outer_audit.md']
    data.update({'status': result['status'], 'core_dependency': 'round008/certify.py',
                 'core_sha256': sha256(source.read_bytes()).hexdigest(),
                 'certifier_sha256': sha256(Path(__file__).read_bytes()).hexdigest(),
                 'analytic_dependencies': {str(p.relative_to(archive)): sha256(p.read_bytes()).hexdigest() for p in dependencies},
                 'limits': {'max_candidate_attempts': 2, 'max_boxes_per_attempt': MAX_BOXES,
                            'seconds_per_attempt': SECONDS, 'tail_terms': 28, 'log_terms': 72},
                 'method': 'Fraction intervals; full-r quartic angular lower bound; exact initial split at 2/3; complete accepted leaves saved',
                 'scope': 'Independent accepted-leaf recomputation by another reviewer remains required'})
    path.write_text(json.dumps(data, indent=2)+'\n')
    print(json.dumps({k: v for k, v in result.items() if k != 'accepted_leaves'}, indent=2), flush=True)
    if result['status'] == 'not_closed':
        raise SystemExit(1)


if __name__ == '__main__':
    main()
