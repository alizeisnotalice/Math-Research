#!/usr/bin/env python3
"""Finite exact audit; no import or evaluation of a heat-kernel implementation."""
from fractions import Fraction as F
from itertools import product
from pathlib import Path
import hashlib
import json
import os
import platform
import time

HERE = Path(__file__).resolve().parent
WORK = Path('/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/upper_g105089')


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def digest(obj):
    return hashlib.sha256(json.dumps(obj, sort_keys=True, separators=(',', ':')).encode()).hexdigest()


def read(path):
    return json.loads(Path(path).read_text())


def key(node):
    root, path = node
    assert isinstance(root, int) and 0 <= root < 52360
    assert isinstance(path, str) and set(path) <= {'0', '1'}
    return root, path


def corners_min(eps, delta, gradient, bounds, anchor):
    eps = F(eps)
    assert eps > 0
    penalty = 1 / (16 * eps * eps)
    margins = []
    for bits in product((0, 1), repeat=4):
        u = [F(bounds[j][bits[j]]) - F(anchor[j]) for j in range(4)]
        margins.append(delta + sum(min(F(a)*x, F(b)*x)
                                  for x, (a, b) in zip(u, gradient))
                       - penalty * sum(x*x for x in u))
    return min(margins)


def bound_key(binding):
    """Illustrative required binding, not a production098 cache implementation."""
    return tuple(binding[k] for k in ('model', 'epsilon', 'point', 'derivative_basis', 'kernel'))


def check_binding(expected, record):
    assert bound_key(expected) == bound_key(record)


def main():
    start = time.monotonic()
    reg = read(HERE/'registration.json')
    for path, wanted in reg['sources'].items():
        assert sha(path) == wanted, path
    obstruction = read(WORK/'fixed_H64_obstruction_batch0009.json')
    state = read(WORK/'batch0009/state.json')
    summary = read(WORK/'batch0009/summary.json')
    validation = read(WORK/'batch0009/validation.json')
    source_reg = read(WORK/'registration.json')
    assert obstruction['epsilon'] == source_reg['epsilon'] == '1/64'
    assert obstruction['target'] == source_reg['target'] == reg['target']
    assert state['source_sha256'] == reg['inherited_source086_hash_not_recomputed_in_099']
    assert summary['state_sha256'] == validation['state_sha256'] == obstruction['overlay_sha256']
    assert not state['faults']
    assert summary['passed'] == 100804 and len(state['new_passed']) == 4700
    queue = [key(n) for n in state['queue']]
    inflight = {key(e['node']): e for e in state['inflight'].values()}
    assert len(queue) == len(set(queue)) == 190 and len(inflight) == 4
    assert not set(queue) & inflight.keys()
    pending = set(queue) | inflight.keys()
    assert len(pending) == 194
    assert len(pending) == obstruction['pending_total'] == summary['queue'] + summary['inflight']
    assert (obstruction['retained_queue'], obstruction['retained_inflight']) == (190, 4)
    for a, b in zip(sorted(pending), sorted(pending)[1:]):
        assert a[0] != b[0] or not b[1].startswith(a[1])
    old_new_passed = [key(row['node']) for row in state['new_passed']]
    assert len(old_new_passed) == len(set(old_new_passed)) == 4700
    for root, path in pending:
        assert all(root != r or (not path.startswith(p) and not p.startswith(path))
                   for r, p in old_new_passed)
    assert validation['tree'] == {'antichain': True, 'internal_count': 48638,
        'leaf_count': 100998, 'roots_with_Kraft_sum_one': 52360}
    assert summary['passed'] + len(pending) == validation['tree']['leaf_count']

    rows = []
    counter_nodes = set()
    target = F(reg['target'])
    for c in obstruction['counterexamples']:
        assert sha(c['path']) == c['sha256']
        record = read(c['path'])
        assert record['record'] == c['record']
        matches = [(node, e) for node, e in inflight.items()
                   if e['point_checkpoint'] == c['path']]
        assert len(matches) == 1
        node, e = matches[0]
        assert e['returned']['kind'] == c['outcome'] == 'counterexample'
        assert e['returned']['record'] == c['record']
        assert e['point'] == c['point'] == c['record']['point']
        assert c['H_interval'] == c['record']['values'][0]
        point = list(map(F, c['point']))
        assert point[0] == 1 and all(0 <= x <= 1 for x in point[1:])
        lower, upper = map(F, c['H_interval'])
        assert 0 < lower <= upper < target
        gap = target - upper
        assert gap == F(c['strict_target_minus_H_upper'])
        counter_nodes.add(node)
        rows.append({'node': list(node), 'point': c['point'], 'H_upper': str(upper),
                     'strict_target_minus_H_upper': str(gap), 'point_sha256': c['sha256']})
    assert len(rows) == 4 and counter_nodes == inflight.keys()
    assert key(state['counterexample']['node']) in counter_nodes
    migration = [{'node': list(node), 'origin': 'inflight_counter_H64' if node in inflight else 'queue',
                  'new_epsilon_intended': '1/128', 'certified_in_099': False}
                 for node in sorted(pending)]
    (HERE/'pending194.json').write_text(json.dumps({
        'scope': 'Frozen expected migration set only; no production state changed',
        'source_state_sha256': obstruction['overlay_sha256'],
        'nodes': migration}, sort_keys=True, indent=2)+'\n')

    constants = []
    for eps in (F(1,64), F(1,128)):
        constants.append({'epsilon': str(eps), 'hessian_magnitude': str(1/(8*eps*eps)),
                          'Taylor_penalty': str(1/(16*eps*eps)), 'log_score_bound': str(1/eps)})
    assert F(constants[1]['Taylor_penalty']) == 4*F(constants[0]['Taylor_penalty'])

    # Synthetic geometry test. delta=1024 and gradients=0 are ASSUMED abstract
    # certificate inputs, never evaluated or asserted for the actual H_epsilon.
    left = [(F(0), F(1,2)), *[(F(0), F(1))]*3]
    right = [(F(1,2), F(1)), *[(F(0), F(1))]*3]
    assert left[0][1] == right[0][0] and left[0][0] == 0 and right[0][1] == 1
    abstract = []
    for eps, box in zip((F(1,64),F(1,128)),(left,right)):
        anchor = [(a+b)/2 for a,b in box]
        margin = corners_min(eps, F(1024), [(0,0)]*4, box, anchor)
        assert margin >= 0
        abstract.append({'epsilon':str(eps),'assumed_Delta':'1024','assumed_gradient':'0',
                         'minimum_log_margin':str(margin),'actual_H_certificate':False})
    assert [r['minimum_log_margin'] for r in abstract] == ['816','192']
    b = {'model':'synthetic_model_identity','epsilon':'1/64','point':'1,1/4,1/2,1/2,1/2',
         'derivative_basis':'physical_dx2_dx3_dx4_dx5','kernel':'synthetic_source_hash'}
    check_binding(b,dict(b))
    rejected=[]
    for field in b:
        wrong=dict(b);wrong[field]+='_wrong'
        try:
            check_binding(b,wrong)
        except AssertionError:
            rejected.append(field)
        else:
            raise AssertionError('missed binding mismatch: '+field)
    eps128=dict(b,epsilon='1/128')
    assert bound_key(b) != bound_key(eps128)
    result = {'status':'passed', 'pid':os.getpid(), 'python':platform.python_version(),
        'registration_sha256':sha(HERE/'registration.json'),'code_sha256':sha(__file__),
        'new_H_calls':0,'new_gradient_calls':0,'old_baseline_recomputations':0,
        'new_log_exp_Arb_calls':0, 'saved_counterexamples':rows,
        'source_passed_count_read':summary['passed'],'source_pending_count_checked':len(pending),
        'source_full_cover_geometry_read_only':validation['tree'],
        'pending194_sha256':sha(HERE/'pending194.json'),
        'epsilon_constants':constants,'abstract_two_leaf_cases':abstract,
        'binding_mismatches_rejected':rejected,
        'g105_closed_cover_certified':False,'production098_audited':False,
        'wall_seconds':time.monotonic()-start}
    (HERE/'results.json').write_text(json.dumps(result,sort_keys=True,indent=2)+'\n')
    print(json.dumps({k:result[k] for k in ('status','pid','wall_seconds','new_H_calls',
          'source_pending_count_checked','g105_closed_cover_certified')},sort_keys=True))


if __name__ == '__main__':
    main()
