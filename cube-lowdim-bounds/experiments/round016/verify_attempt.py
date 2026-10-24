#!/usr/bin/env python3
"""Review the failed attempt without optimization or full leaf recomputation.

Default: verify the archived summary, its auxiliary scalars and exactly three
selected leaves. --full FILE additionally checks the local full binary tree.
The two modes have distinct status/output files; neither proves global G>=0.
"""
from pathlib import Path
from fractions import Fraction as Q
from hashlib import sha256
import argparse
import importlib.util
import json
import sys

sys.dont_write_bytecode = True
ROOT = Path(__file__).resolve().parent
ARCHIVE = ROOT.parent.parent/'Math-Research-lowdim-agent/cube-lowdim-bounds/experiments'
R = Q(15, 8)
E = {2: Q(69, 100), 6: Q(19, 50), 7: Q(11, 10), 8: Q(3, 4)}


def require(value, label):
    if not value:
        raise ValueError(label)


def box(row):
    return tuple(Q(row[key]) for key in ('r_lo', 'r_hi', 'q_lo', 'q_hi'))


def digest(data):
    return sha256(json.dumps(data, sort_keys=True, separators=(',', ':')).encode()).hexdigest()


def load_core(data, explicit):
    candidates = [explicit] if explicit else [ROOT.parent/'round008/certify.py', ARCHIVE/'round008/certify.py']
    matches = [p for p in candidates if p.is_file() and sha256(p.read_bytes()).hexdigest() == data['core_sha256']]
    require(bool(matches), 'reviewed core source/hash; use --core if relocated')
    spec = importlib.util.spec_from_file_location('round016_independent_attempt_log', matches[0])
    core = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = core
    spec.loader.exec_module(core)
    core.R = R
    return core


def verify_auxiliary(core, data):
    aux = data['auxiliary']
    harmonic = sum((Q(1, k) for k in range(1, 8)), Q(0))
    small = 2*Q(101, 100)**6*harmonic+core.log_point(Q(1, 101))+Q(99, 100)**7*core.log_point(Q(1, 99))
    require(small.pair() == aux['small_radius']['B_interval'] and small.hi < 0, 'exact common small-radius interval')
    log2 = core.log_point(2)
    hp = {2: 8-8*log2, 6: Q(1334, 5)-384*log2,
          7: Q(3108, 5)-896*log2, 8: Q(149104, 105)-2048*log2}
    d1 = 2*(2*core.tail_series(core.I.point(1))+core.I.point(1)/R)/(R-1)
    exterior = d1-sum((E[p]*hp[p] for p in E), core.I.point(0))
    require({str(p): h.pair() for p, h in hp.items()} == aux['exterior']['Hp_intervals'], 'endpoint constant intervals')
    require(d1.pair() == aux['exterior']['D_R1_interval'], 'tail endpoint interval')
    require(exterior.pair() == aux['exterior']['new_slack_interval'] and exterior.lo > 0, 'new exterior condition')
    dominance = R-(R-1)*sum(E.values())/4
    cost = R*R/(R-1)-sum((2*e/Q((p+1)*(p+2)) for p, e in E.items()), Q(0))
    old = Q(411881, 106500)
    require(dominance == Q(aux['kernel_dominance']) == Q(989, 800), 'kernel dominance')
    require(cost == Q(aux['cost']) == Q(4841, 1260), 'candidate exact cost')
    require(old == Q(aux['previous_strict_upper']) and old-cost == Q(aux['improvement']) == Q(28363, 1118250), 'candidate improvement')
    return {'exact_auxiliary_intervals_reproduced': True, 'small_B_upper_display': float(small.hi),
            'exterior_lower_display': float(exterior.lo), 'candidate_cost_only': str(cost),
            'rigorous_upper_unchanged': str(old)}


def recompute_three(core, selected):
    require([x['accepted_leaf_index'] for x in selected] == [0, 7494, 14988], 'fixed first/middle/last selection')
    require([x['node_id'] for x in selected] == [60, 15034, 30022], 'fixed selected nodes')
    def w_lower(l, h, n):
        if h == 0:
            return Q(0)
        if l == 0:
            return h**n*core.log_point(h).lo if h <= Q(1, 4) else -Q(1, 2*n)
        return min(x**n*y for x in (l, h) for y in (core.log_point(l).lo, core.log_point(h).hi))
    def negative_log_upper(c, x):
        return Q(0) if c == 0 or x == 0 else c*min(Q(0), core.log_point(x).hi)
    checks = []
    for leaf in selected:
        l, h, a, b = box(leaf)
        require(Q(1, 100) <= l < h <= 1 and 0 <= a < b <= 1, 'selected leaf domain')
        dl, dh = l*a, h*b
        ul, uh = l*(1+a)/2, h*(1+b)/2
        vl, vh = l*(1-b)/2, h*(1-a)/2
        sl, sh = 1-h, 1-l
        aa, bb = 1-vl, 1-ul
        zp, zP, zm, zM = 1+dl, 1+dh, 1-dh, 1-dl
        tail = 2*(2*sum((l**j/((j+2)*R**(j+1)) for j in range(29)), Q(0))
                  +2*sum((dl**(2*j)/((2*j+2)*R**(2*j+1)) for j in range(15)), Q(0)))/(R-1)
        uppers = {}
        for p in E:
            n = p-1
            polynomial = sum(((aa**(n-k)*(zP**k+sh**k)+bb**(n-k)*(zM**k+sh**k))/Q(n-k)
                              for k in range(n)), Q(0))
            uppers[p] = 2*p*(polynomial+negative_log_upper(zp**n+sl**n, uh)
                             +negative_log_upper(zm**n+sl**n, vh)-w_lower(zp, zP, n)
                             -w_lower(zm, zM, n)-2*w_lower(sl, sh, n))
        value = tail-sum((E[p]*uppers[p] for p in E), Q(0))
        require(str(tail) == leaf['tail_lower'], 'selected exact tail lower')
        require({str(p): str(v) for p, v in uppers.items()} == leaf['bump_upper_by_power'], 'selected exact bump uppers')
        require(str(value) == leaf['generator_lower'] and value >= 0, 'selected exact nonnegative total')
        checks.append({'leaf_index': leaf['accepted_leaf_index'], 'node_id': leaf['node_id'],
                       'box': [str(x) for x in box(leaf)], 'lower_display': float(value),
                       'exact_record_reproduced': True})
    return checks


def verify_full_tree(data):
    nodes, leaves, queue = data['partition_tree'], data['accepted_leaf_data'], data['unresolved_queue']
    require(box(nodes[0]) == (Q(1, 100), Q(1), Q(0), Q(1)) and nodes[0]['parent'] is None, 'exact full-tree root')
    seen, todo, terminal, splits = set(), [0], [], 0
    while todo:
        i = todo.pop()
        require(i not in seen and nodes[i]['id'] == i, 'unique reachable node identity')
        seen.add(i)
        node = nodes[i]
        l, h, a, b = box(node)
        require(Q(1, 100) <= l < h <= 1 and 0 <= a < b <= 1, 'positive tree box domain')
        if node['status'] != 'split':
            terminal.append(i)
            continue
        splits += 1
        ids, mid = node['children'], Q(node['midpoint'])
        require(len(ids) == 2 and ids[0] != ids[1], 'two distinct children')
        if node['axis'] == 'r':
            require(mid == (l+h)/2, 'exact r midpoint')
            children = [(l, mid, a, b), (mid, h, a, b)]
        else:
            require(node['axis'] == 'q' and mid == (a+b)/2, 'exact q midpoint')
            children = [(l, h, a, mid), (l, h, mid, b)]
        for child, expected in zip(ids, children):
            require(nodes[child]['parent'] == i and box(nodes[child]) == expected, 'exact parent/child binary cover')
            todo.append(child)
    require(len(seen) == len(nodes) == 1+2*splits, 'all nodes reachable and binary count')
    accepted = {x['node_id'] for x in leaves}
    pending = {x['node_id'] for x in queue}
    require(len(accepted) == len(leaves) and len(pending) == len(queue) and not accepted & pending, 'unique disjoint terminal lists')
    require(accepted | pending == set(terminal), 'every terminal accounted for')
    for leaf in leaves:
        require(nodes[leaf['node_id']]['status'] == 'accepted' and box(nodes[leaf['node_id']]) == box(leaf), 'accepted leaf identity')
        require(Q(leaf['generator_lower']) >= 0, 'recorded accepted sign, not recomputation')
    for leaf in queue:
        require(box(nodes[leaf['node_id']]) == box(leaf), 'unresolved terminal identity')
    inner = data['inner']
    require(len(leaves) == inner['accepted_leaves'] and len(queue) == inner['unresolved_boxes'] and splits == inner['splits'], 'saved counts')
    evaluated = len(leaves)+splits+sum(nodes[i]['status'] == 'unresolved_width_floor' for i in pending)
    require(evaluated == inner['actual_evaluated_boxes'] == 30000, 'actual evaluation accounting')
    require(digest(leaves) == data['accepted_leaf_sha256'] and digest(nodes) == data['partition_tree_sha256'], 'full leaf and tree digests')
    def area(rows):
        return sum(((x[1]-x[0])*(x[3]-x[2]) for x in map(box, rows)), Q(0))
    accepted_area, pending_area = area(leaves), area(queue)
    require(accepted_area+pending_area == Q(99, 100), 'exact root area ledger')
    require(min(Q(x['generator_lower']) for x in leaves) == Q(inner['minimum_accepted_lower']), 'saved minimum from all recorded signs')
    return {'full_tree_structure_checked_this_run': True, 'tree_nodes': len(nodes),
            'accepted_leaves': len(leaves), 'splits': splits, 'unresolved_terminals': len(queue),
            'exact_root_cover_including_unresolved': True,
            'accepted_parameter_area': str(accepted_area), 'unresolved_parameter_area': str(pending_area),
            'accepted_fraction_of_root': str(accepted_area/Q(99, 100)),
            'full_leaf_formula_recomputation': False}


def main():
    parser = argparse.ArgumentParser()
    inputs = parser.add_mutually_exclusive_group()
    inputs.add_argument('--summary', type=Path, default=None)
    inputs.add_argument('--full', type=Path)
    parser.add_argument('--core', type=Path)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    source = args.full or args.summary or ROOT/'certificate_summary.json'
    data = json.loads(source.read_text())
    require(data['round'] == 16 and data['status'] == 'not_closed' and data['candidate_attempts'] == 1, 'failed single-attempt status')
    require(data['parameters'] == {'R': '15/8', 'e_by_power': {str(p): str(e) for p, e in E.items()}}, 'unique authorized rational candidate')
    require(data['inner']['stop_reason'] == 'box_evaluation_budget_exhausted' and not data['inner']['complete_tree_cover'], 'budget stop, no certified cover')
    require(0 <= data['inner']['elapsed_seconds'] < 285 and data['inner']['actual_evaluated_boxes'] == 30000, 'recorded finite limits')
    require(sha256((ROOT/'certify.py').read_bytes()).hexdigest() == data['certifier_sha256'], 'certifier hash')
    for name, wanted in data['analytic_dependencies'].items():
        require(any(p.is_file() and sha256(p.read_bytes()).hexdigest() == wanted
                    for p in (ROOT.parent/name, ARCHIVE/name)), 'analytic dependency hash: '+name)
    core = load_core(data, args.core)
    auxiliary = verify_auxiliary(core, data)
    if args.full:
        structure = verify_full_tree(data)
        leaves = data['accepted_leaf_data']
        selected = [{'accepted_leaf_index': i, **leaves[i]} for i in [0, len(leaves)//2, len(leaves)-1]]
        status = 'failed_attempt_full_tree_structure_and_three_leaves_verified'
    else:
        require('partition_tree' not in data and 'accepted_leaf_data' not in data, 'summary mode does not ingest a full tree')
        selected = data['selected_accepted_leaves']
        structure = {'full_tree_structure_checked_this_run': False,
                     'reason': 'Only the compact summary is available in this mode; full-tree evidence requires --full FILE.',
                     'full_leaf_formula_recomputation': False}
        status = 'failed_attempt_summary_three_leaves_and_auxiliary_verified'
    checks = recompute_three(core, selected)
    result = {'status': status, 'mode': 'full_tree_plus_three_leaves' if args.full else 'summary_only',
              'candidate_status': 'not_closed', 'new_LPs': 0, 'new_adaptive_certification_runs': 0,
              **structure, **auxiliary, 'independently_recomputed_leaf_count': 3, 'checks': checks,
              'saved_full_leaf_sha256': data['accepted_leaf_sha256'],
              'saved_partition_tree_sha256': data['partition_tree_sha256'],
              'input_sha256': sha256(source.read_bytes()).hexdigest(),
              'verifier_sha256': sha256(Path(__file__).read_bytes()).hexdigest(),
              'scope': 'A failed finite attempt; three accepted formulas checked. No complete leaf revalidation or new general upper bound.'}
    output = args.output or ROOT/('attempt_verification.json' if args.full else 'summary_verification.json')
    output.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
