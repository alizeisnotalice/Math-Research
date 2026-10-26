#!/usr/bin/env python3
"""Finite frozen batch9 check: counterpoints first; durable, unrepeated prefix."""
import sys
sys.dont_write_bytecode = True
import gzip, json, os, time
from pathlib import Path
from fractions import Fraction as F
from interfaces import REG, HERE, common as C, arith as A


def atomic(path, value):
    path = Path(path); temporary = Path(str(path) + '.tmp')
    with temporary.open('wb') as stream:
        stream.write(json.dumps(value, sort_keys=True, indent=2).encode() + b'\n')
        stream.flush(); os.fsync(stream.fileno())
    temporary.replace(path)


def decimal_outward(value, digits=18):
    scale = 10**digits
    return [str(F(value[0].numerator*scale//value[0].denominator, scale)),
            str(F(-((-value[1].numerator*scale)//value[1].denominator), scale))]


def main():
    began = time.monotonic()
    registration_path = HERE/'obstruction9_values_registration.json'
    registration = json.loads(registration_path.read_text())
    for path, digest in registration['input_locks'].items():
        assert C.file_sha(path) == digest, ('input changed', path)
    for path, digest in REG['protected_old_locks'].items():
        assert C.file_sha(path) == digest, ('old lock changed', path)
    manifest = json.loads(Path(registration['point_manifest']).read_text())
    assert manifest['target'] == REG['target'] == '21/200000'
    assert manifest['epsilon'] == REG['epsilon'] == '1/64'
    points = manifest['points']; assert len(points) == manifest['new_points'] == 12
    assert sum(p['outcome'] == 'counterexample' for p in points) == 4
    point_keys = {','.join(p['point']) for p in points}
    assert len(point_keys) == 12
    reader = C.Stream(REG['source086'])
    while True:
        key = reader.key(); assert key is not None
        if key == 'cache':
            old_keys = {key for key, _ in reader.items(True)}
            reader.f.close(); break
        reader.value()
    assert len(old_keys) == 180040 and not point_keys.intersection(old_keys)
    prior_keys = set()
    for path in registration['prior_numeric_registrations']:
        previous = json.loads(Path(path).read_text())
        assert C.file_sha(previous['point_manifest']) == previous['point_manifest_sha256']
        previous_manifest = json.loads(Path(previous['point_manifest']).read_text())
        prior_keys.update(','.join(p['point']) for p in previous_manifest['points'])
    assert len(prior_keys) == 4744 and not point_keys.intersection(prior_keys)
    fan = json.loads(gzip.decompress(Path(REG['fan']).read_bytes()))
    order = sorted(range(12), key=lambda i: (points[i]['outcome'] != 'counterexample', i))
    assert registration['evaluation_order'] == order
    result_path = HERE/'obstruction9_values_result.json'
    if result_path.exists():
        result = json.loads(result_path.read_text())
        assert result['registration_sha256'] == C.file_sha(registration_path)
    else:
        result = {'task': 'audit_g105090_obstruction9', 'registration_sha256': C.file_sha(registration_path),
                  'status': 'pending_values', 'points': [], 'old180040_recomputed': False,
                  'prior4744_recomputed': False, 'all5_components': True,
                  'method': '081 readonly direct rational per-cone partial fractions, integer768 log and Machin pi; saved512 outward',
                  'scope': 'fixed H_(1/64) values only; no full G counterexample or Cn conclusion',
                  'runs': []}
    this_run = {'pid': os.getpid(), 'initial_verified': len(result['points'])}
    result['runs'].append(this_run); atomic(result_path, result)
    target = F(REG['target'])
    for position, index in enumerate(order):
        job = points[index]
        assert C.file_sha(job['path']) == job['sha256']
        raw = json.loads(Path(job['path']).read_text())
        assert raw['stage'] == job['outcome'] and raw['stage'] in ('complete', 'counterexample')
        assert raw['partial'] is None and raw['record'] == job['record'] and raw['point'] == job['point']
        assert len(job['record']['values']) == 5
        point = list(map(F, job['point']))
        assert point[0] == 1 and all(0 <= x <= 1 for x in point)
        assert all(point[j] >= point[j+1] for j in range(4)), 'not closed ordered face'
        if position < len(result['points']):
            proof = result['points'][position]
            assert proof['point_manifest_index'] == index and proof['point'] == job['point']
            assert proof['raw_checkpoint_sha256'] == job['sha256']
            for ball, production in zip(proof['independent_values_512_dyadic'], job['record']['values']):
                lo, hi = map(F, ball); L, U = map(F, production); assert L < lo <= hi < U
            continue
        start = time.monotonic()
        values, geometry = A.direct_values(job['point'], fan)
        saved = []
        for interval, production in zip(values, job['record']['values']):
            L, U = map(F, production); assert L < interval[0] <= interval[1] < U
            output = A.dyadic_outward(interval); lo, hi = map(F, output)
            assert L < lo <= hi < U
            saved.append(list(output))
        hlo, hhi = map(F, saved[0]); assert hlo > 0
        production_hi = F(job['record']['values'][0][1])
        counter = job['outcome'] == 'counterexample'
        assert counter == (production_hi < target)
        if counter: assert hhi < target
        proof = {'point_manifest_index': index, 'point': job['point'], 'raw_checkpoint_sha256': job['sha256'],
                 'raw_outcome': job['outcome'], 'independent_values_512_dyadic': saved,
                 'all5_strict_containment': True, 'closed_ordered_face': True,
                 'H_upper_below_target': counter, 'target_minus_independent_H_upper': str(target-hhi),
                 'H_micro_decimal18_outward_rational': decimal_outward((hlo*10**6, hhi*10**6)),
                 'geometry': geometry, 'seconds': time.monotonic()-start, 'pid': os.getpid()}
        result['points'].append(proof)
        result['verified_points'] = len(result['points'])
        result['verified_counters'] = sum(p['H_upper_below_target'] for p in result['points'])
        result['status'] = 'four_fixed_H64_obstructions_independently_confirmed' if result['verified_counters'] == 4 else 'pending_values'
        atomic(result_path, result)
        print(json.dumps({'pid': os.getpid(), 'manifest_index': index, 'verified': len(result['points']),
                          'verified_counters': result['verified_counters'], 'H_micro': proof['H_micro_decimal18_outward_rational']}), flush=True)
    assert len(result['points']) == 12 and result['verified_counters'] == 4
    for path, digest in registration['input_locks'].items(): assert C.file_sha(path) == digest
    for path, digest in REG['protected_old_locks'].items(): assert C.file_sha(path) == digest
    this_run.update({'final_verified': 12, 'wall_seconds': time.monotonic()-began, 'exit_intent': 0})
    result['status'] = 'all12_H_plus4grad_verified_four_fixed_H64_obstructions'
    atomic(result_path, result)
    print(json.dumps({'pid': os.getpid(), 'status': result['status'], 'wall_seconds': this_run['wall_seconds'],
                      'result_sha256': C.file_sha(result_path)}), flush=True)


if __name__ == '__main__': main()
