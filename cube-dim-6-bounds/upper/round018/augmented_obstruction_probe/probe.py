#!/usr/bin/env python3
"""B018: prepare exact new-frame tuple, then at most one authorized native call."""
import sys
sys.dont_write_bytecode = True
import argparse
import hashlib
import itertools
import json
import math
import os
import subprocess
import time
from fractions import Fraction as Q
from pathlib import Path

HERE = Path(__file__).resolve().parent
PROJECT = HERE.parents[2]
ENGINE = Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube6-audit/cube-dim-6-bounds/upper/round014/generic_short_engine')
AUG = Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube6-curvature-astra/cube-dim-6-bounds/upper/round012/augmented_frame')
LIT = Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube6-literature/cube-dim-6-bounds')
BASE_FAN = PROJECT/'upper/round003/direct_kernel/fan_integer_input.txt'
FAN = AUG/'fan_beta1.txt'
POINT = PROJECT/'upper/round015/spatial_minimum_search/checkpoint.json'
SOURCE_SHA = '247a1f6c6967ad52f34320b999818214488700651b46cb8e5feb3eac7064953b'
BINARY_SHA = '67119aa799041ed75b66f043444640a2a62556f58517cbdf3bd6f7cccd564f04'
FAN_SHA = '450caf04d124334474b91f3d9fdabce4eadf3259ea47def9094ae5eb056b6c4e'
BASE_SHA = '7a57caad9b7713d9ce0896e51ab1dd402566feb1cba50829104c593a6d12a812'

def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def dump(path, obj):
    Path(path).write_text(json.dumps(obj, indent=2, ensure_ascii=False)+'\n')

def read_fan(path):
    it = iter(map(int, Path(path).read_text().split()))
    n = next(it)
    assert n == 21
    cones = []
    for _ in range(n):
        determinant = next(it)
        rows = [(next(it), tuple(next(it) for _ in range(6))) for _ in range(6)]
        cones.append((determinant, rows))
    assert next(it, None) is None
    return cones

def frame():
    return [(1, *s) for s in itertools.product((-1, 1), repeat=5)] + [tuple(int(i == j) for i in range(6)) for j in range(6)]

def zero_forms(rays):
    forms = set()
    for ray in rays:
        for perm in set(itertools.permutations(ray)):
            nz = [i for i, v in enumerate(perm) if v]
            for signs in itertools.product((-1, 1), repeat=len(nz)):
                a = list(perm)
                for i, s in zip(nz, signs):
                    a[i] *= s
                g = math.gcd(*a)
                a = tuple(v//g for v in a)
                if next(v for v in a if v) < 0:
                    a = tuple(-v for v in a)
                forms.add(a)
    return sorted(forms)

def prepare():
    import sympy
    import flint
    assert not (HERE/'input_manifest.json').exists(), 'Do not overwrite a registered batch.'
    tic = time.perf_counter()
    started = time.time()
    assert sha(BASE_FAN) == BASE_SHA and sha(FAN) == FAN_SHA
    assert sha(ENGINE/'generic_short.cpp') == SOURCE_SHA
    assert sha(ENGINE/'scratch/generic_short') == BINARY_SHA
    base, aug = read_fan(BASE_FAN), read_fan(FAN)
    vs = frame()
    audits = []
    for index, ((det0, rows0), (det1, rows1)) in enumerate(zip(base, aug)):
        assert det0 == det1 > 0
        assert [r for _, r in rows0] == [r for _, r in rows1]
        actual_det = abs(int(sympy.Matrix.hstack(*[sympy.Matrix(r) for _, r in rows1]).det()))
        assert actual_det == det1
        for (q0, r), (q1, _) in zip(rows0, rows1):
            assert all(r[i] >= r[i+1] >= 0 for i in range(5))
            explicit_old = sum(abs(sum(a*b for a, b in zip(v, r))) for v in vs[:32])
            explicit_new = sum(abs(sum(a*b for a, b in zip(v, r))) for v in vs)
            assert q0 == explicit_old and q1 == q0+sum(map(abs, r)) == explicit_new
            assert 0 < q1 <= 1000000 and all(abs(x) <= 1000000 for x in r)
        for v in vs:
            products = [sum(a*b for a, b in zip(v, r)) for _, r in rows1]
            assert min(products) >= 0 or max(products) <= 0
        audits.append({'cone_index': index, 'determinant': det1,
                       'old_q': [q for q, _ in rows0], 'new_q': [q for q, _ in rows1],
                       'rays_unchanged': True, 'all38_cone_forms_constant_weak_sign': True})
    prior = json.loads(POINT.read_text())
    assert prior['counterexample_native_index'] == 2
    point = prior['native_results'][1]
    nums, den = point['numerators'], point['denominator']
    assert nums == [100000000000000,99999998999999,99999997999988,90268521812170,53714088833376,157107163435]
    assert den == 10**14 and nums[0] == den and all(nums[i] > nums[i+1] > 0 for i in range(5))
    unique_rays = sorted(set(r for _, rows in aug for _, r in rows))
    assert len(unique_rays) == 18
    forms = zero_forms(unique_rays)
    dots = [sum(a*n for a, n in zip(f, nums)) for f in forms]
    assert all(dots)
    max_num = max(sum(abs(a)*max(nums) for a in r) for r in unique_rays)
    max_den = den*max(q for _, rows in aug for q, _ in rows)
    assert max_num < 2**63-1 and max_den <= 2**63-1 and max_num*max_den < 2**127-1
    old_manifest = json.loads((AUG/'input_manifest.json').read_text())
    assert old_manifest['beta'] == '1' and old_manifest['T'] == '2/45'
    assert all(tuple(Q(n, den) for n in nums) != tuple(Q(n, p['denominator']) for n in p['numerators']) for p in old_manifest['points'])
    model = {'dimension': 6, 'beta': '1', 'directions': vs, 'weights': [1]*38,
             'symbol': 'sum_body32 abs(v.xi)+sum_i abs(xi_i)', 'measure': 'Lebesgue',
             'fan_sha256': FAN_SHA, 'Fourier_inverse_prefactor': '(2*pi)^-6'}
    frame_sha = hashlib.sha256(json.dumps(model, sort_keys=True, separators=(',', ':')).encode()).hexdigest()
    build = json.loads((ENGINE/'build_record.json').read_text())
    files = [BASE_FAN, FAN, POINT, PROJECT/'upper/round015/spatial_minimum_search/report.md',
             ENGINE/'generic_short.cpp', ENGINE/'scratch/generic_short', ENGINE/'build_record.json',
             ENGINE/'report.md', AUG/'report.md', AUG/'theorem.md', AUG/'input_manifest.json',
             LIT/'audit/round012/augmented_frame_theory/report.md',
             LIT/'audit/round014/generic_short_engine_review/report.md', HERE/'registration.md', Path(__file__).resolve()]
    for path, expected in build['runtime_dylib_sha256'].items():
        assert sha(path) == expected
        files.append(Path(path))
    output = {'status': 'REGISTERED_AND_RELEASED_BEFORE_SOLE_NATIVE_PROBE',
              'started_epoch': started, 'beta': '1', 'T': '2/45', 'T_numerator': 2, 'T_denominator': 45,
              'numerators': nums, 'denominator': den, 'max_native_attempts': 1,
              'max_native_seconds': 120, 'max_total_seconds': 1200,
              'model': model, 'frame_sha256': frame_sha,
              'canonical_point_key': [frame_sha, '2/45', nums, den],
              'same_frame_previous_two_points_are_distinct': True,
              'inherits_base_kernel_value': False,
              'preflight': {'cone_count': 21, 'ray_occurrences': 126, 'unique_rays': 18,
                            'zero_form_count': len(forms), 'zero_nodes': 0, 'minimum_abs_primitive_form': min(map(abs, dots)),
                            'max_raw_node_numerator_bound': max_num, 'max_raw_node_denominator': max_den,
                            'comparator_int128_product_bound': max_num*max_den, 'cone_audits': audits},
              'command': [str(ENGINE/'scratch/generic_short'), 'full', str(FAN), '2', '45', str(den), *map(str, nums)],
              'inputs': [{'path': str(p), 'sha256': sha(p)} for p in files],
              'dependencies': {'python': sys.version, 'sympy': sympy.__version__, 'python_flint': flint.__version__},
              'preparation_elapsed_seconds': time.perf_counter()-tic, 'preparation_pid': os.getpid(),
              'native_calls_so_far': 0}
    dump(HERE/'input_manifest.json', output)
    print(json.dumps({k:output[k] for k in ('status','frame_sha256','numerators','denominator','preparation_elapsed_seconds','preparation_pid')}, indent=2))
    print(json.dumps({k:v for k,v in output['preflight'].items() if k != 'cone_audits'}, indent=2))

def run():
    assert not (HERE/'run_receipt.json').exists(), 'One attempt only; no replay.'
    manifest = json.loads((HERE/'input_manifest.json').read_text())
    assert manifest['status'] == 'REGISTERED_AND_RELEASED_BEFORE_SOLE_NATIVE_PROBE'
    for source in manifest['inputs']:
        assert sha(source['path']) == source['sha256']
    remaining = manifest['max_total_seconds']-(time.time()-manifest['started_epoch'])
    assert remaining > 0, 'registered total time expired before execution'
    record = {'status': 'NATIVE_ATTEMPT_REGISTERED', 'pid': os.getpid(), 'native_attempts': 1,
              'input_manifest_sha256': sha(HERE/'input_manifest.json'), 'started_epoch': time.time(),
              'command': manifest['command'], 'additional_native_calls_authorized': 0}
    dump(HERE/'run_receipt.json', record)
    tic = time.perf_counter()
    stdout, stderr = HERE/'native.stdout.json', HERE/'native.stderr.txt'
    process = None
    try:
        with stdout.open('w') as out, stderr.open('w') as err:
            process = subprocess.Popen(manifest['command'], stdout=out, stderr=err, cwd=HERE)
            record['native_pid'] = process.pid
            record['status'] = 'NATIVE_RUNNING'
            dump(HERE/'run_receipt.json', record)
            try:
                code = process.wait(timeout=min(remaining, manifest['max_native_seconds']))
            except subprocess.TimeoutExpired:
                process.kill()
                code = process.wait()
                record['timed_out'] = True
            except BaseException:
                if process.poll() is None:
                    process.kill()
                record['returncode'] = process.wait()
                raise
        record['returncode'] = code
        record['status'] = 'SINGLE_NATIVE_EXIT0_STOPPED' if code == 0 else 'SINGLE_NATIVE_FAILED_STOPPED'
    finally:
        if process is not None and process.poll() is None:
            process.kill()
            process.wait()
        record['elapsed_seconds'] = time.perf_counter()-tic
        record['all_native_children_waited'] = True
        record['stdout_sha256'] = sha(stdout) if stdout.exists() else None
        record['stderr_sha256'] = sha(stderr) if stderr.exists() else None
        dump(HERE/'run_receipt.json', record)
    print(json.dumps(record, indent=2))
    if record['status'] != 'SINGLE_NATIVE_EXIT0_STOPPED':
        raise SystemExit(2)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('mode', choices=('prepare', 'run'))
    args = parser.parse_args()
    (prepare if args.mode == 'prepare' else run)()
