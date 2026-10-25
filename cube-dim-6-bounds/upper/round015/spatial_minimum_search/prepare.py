#!/usr/bin/env python3
"""B015 preparation only: exact inputs, Fourier-node preflight, no native call."""
import sys
sys.dont_write_bytecode = True
import hashlib
import itertools
import json
import math
import platform
from fractions import Fraction as Q
from pathlib import Path
from time import time

HERE = Path(__file__).resolve().parent
PROJECT = HERE.parents[2]
ENGINE = Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube6-audit/cube-dim-6-bounds/upper/round014/generic_short_engine')
FAN = PROJECT / 'upper/round003/direct_kernel/fan_integer_input.txt'
POINT2 = PROJECT / 'upper/round010/short_minimum_diagnostic/point2.json'
OLD = PROJECT / 'upper/round011/cutoff48/input_manifest.json'
DEN = 10**14
GAP = 10**6  # 1e-8 on the DEN lattice
INT64_MAX = 2**63 - 1
GAP_WEIGHTS = tuple(11**i for i in range(6))
PERTURBATION_K = (1, 2, 3, 4, 5, 6)

def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def dump(path, obj):
    Path(path).write_text(json.dumps(obj, indent=2, ensure_ascii=False) + '\n')

def nearest(q):
    """Nearest integer, exact; ties to even."""
    a, b = divmod(q.numerator, q.denominator)
    if 2*b < q.denominator:
        return a
    if 2*b > q.denominator:
        return a+1
    return a + (a % 2)

def fan_data():
    it = iter(map(int, FAN.read_text().split()))
    n = next(it)
    assert n == 21
    rays = set()
    cones = []
    for _ in range(n):
        det = next(it)
        rows = [(next(it), tuple(next(it) for _ in range(6))) for _ in range(6)]
        rays.update(rows)
        cones.append((det, rows))
    assert next(it, None) is None
    assert len(rays) == 18
    return cones, sorted(rays)

def all_zero_forms(rays):
    """Exact full hyperoctahedral orbit, primitive and modulo global sign.

    Every native node numerator is +/- a positive integer times one form.
    A denominator q cannot affect its zero locus. Native uses 32 sign pairs;
    including all 64 signs here and quotienting global sign is identical for
    zero detection, with no floating tolerance.
    """
    forms = set()
    for _, ray in rays:
        for perm in set(itertools.permutations(ray)):
            nz = [i for i, x in enumerate(perm) if x]
            for bits in itertools.product((-1, 1), repeat=len(nz)):
                a = list(perm)
                for i, s in zip(nz, bits):
                    a[i] *= s
                g = math.gcd(*a)
                a = tuple(v//g for v in a)
                if next(v for v in a if v) < 0:
                    a = tuple(-v for v in a)
                forms.add(a)
    return sorted(forms)

def feasible(nums, den):
    return (len(nums) == 6 and nums[0] == den and
            Q(nums[-1], den) >= Q(1, 10**8) and
            all(Q(nums[i]-nums[i+1], den) >= Q(1, 10**8)
                for i in range(5)))

def preflight(nums, den, forms, rays):
    assert all(type(v) is int for v in nums) and type(den) is int
    dots = [sum(a*b for a, b in zip(form, nums)) for form in forms]
    zeros = [list(form) for form, v in zip(forms, dots) if v == 0]
    max_n = max(sum(abs(a)*max(map(abs, nums)) for a in ray) for _, ray in rays)
    max_d = den * max(q for q, _ in rays)
    return {'feasible': feasible(nums, den), 'zero_form_count': len(zeros),
            'zero_forms': zeros, 'min_abs_primitive_form': min(map(abs, dots)),
            'int64_safe': (0 < den <= INT64_MAX and max_n < INT64_MAX and
                           max_d <= INT64_MAX and max(map(abs, nums)) < INT64_MAX),
            'conservative_max_abs_raw_numerator': max_n,
            'max_raw_node_denominator': max_d}

def gap_repairs(nums, den):
    """Six prescribed exact NEW points; no kernel calls and no auto-acceptance.

    Gap i=0 is 1-x2 and gap i=5 is x6. Put the mass subtraction in
    the exact largest gap (first index on ties) and add to every other gap.
    Every k is based on the same rounded input, never on the previous repair.
    """
    base = [nearest(Q(n, den)*DEN) for n in nums]
    gaps = [base[i]-base[i+1] for i in range(5)] + [base[-1]]
    j = max(range(6), key=lambda i: gaps[i])
    for k in PERTURBATION_K:
        delta = [k*w if i != j else 0 for i, w in enumerate(GAP_WEIGHTS)]
        delta[j] = -sum(delta)
        moved = [g+d for g, d in zip(gaps, delta)]
        candidate = [DEN] + [DEN-sum(moved[:i]) for i in range(1, 6)]
        yield {'k': k, 'largest_gap_index': j,
               'quantized_base_numerators': base, 'quantized_base_denominator': DEN,
               'base_rounding_errors': [str(Q(a, DEN)-Q(n, den)) for a, n in zip(base, nums)],
               'base_is_same_exact_point': all(Q(a, DEN) == Q(n, den) for a, n in zip(base, nums)),
               'base_gap_numerators': gaps, 'gap_delta_numerators': delta,
               'candidate_gap_numerators': moved,
               'numerators': candidate, 'denominator': DEN,
               'is_new_exact_input': True}

def make_preparation():
    import scipy
    import numpy
    import flint
    tic = time()
    cones, rays = fan_data()
    forms = all_zero_forms(rays)
    p2 = json.loads(POINT2.read_text())
    old = json.loads(OLD.read_text())['point_order'][1]
    seed1 = [nearest(Q(n, p2['denominator'])*DEN) for n in p2['candidate_numerators']]
    seeds = [
        {'label': 'B010_point2_nearest_den_1e14_NEW_INPUT', 'numerators': seed1,
         'denominator': DEN, 'source_path': str(POINT2), 'source_sha256': sha(POINT2),
         'coordinate_rounding_errors': [str(Q(n, DEN)-Q(o, p2['denominator']))
            for n, o in zip(seed1, p2['candidate_numerators'])],
         'inherits_prior_point_certificate': False},
        {'label': 'A007_reference_same_point_NEW_T48_GRADIENT',
         'numerators': old['numerators'], 'denominator': old['denominator'],
         'source_path': str(OLD), 'source_sha256': sha(OLD),
         'inherits_prior_gradient_certificate': False}
    ]
    for seed in seeds:
        seed['preflight'] = preflight(seed['numerators'], seed['denominator'], forms, rays)
        assert seed['preflight']['feasible'] and seed['preflight']['int64_safe']
        seed['canonical_rational_point'] = [str(Q(n, seed['denominator'])) for n in seed['numerators']]
    frame = {'dimension': 6, 'spectral_directions': [[1, *s] for s in itertools.product((-1, 1), repeat=5)],
             'direction_weights': [1]*32, 'fourier_exponent': 'sum_v abs(dot(v,xi))',
             'coordinate_measure': 'Lebesgue', 'fan_sha256': sha(FAN)}
    frame_hash = hashlib.sha256(json.dumps(frame, sort_keys=True, separators=(',', ':')).encode()).hexdigest()
    forms_bytes = json.dumps(forms, separators=(',', ':')).encode()
    sources = [FAN, POINT2, OLD, ENGINE/'generic_short.cpp', ENGINE/'build_record.json',
               ENGINE/'native_result.json', ENGINE/'validation_results.json',
               ENGINE/'scratch/generic_short', Path(__file__).resolve()]
    build = json.loads((ENGINE/'build_record.json').read_text())
    sources += [Path(p) for p in build['runtime_dylib_sha256']]
    inputs = [{'path': str(p), 'sha256': sha(p)} for p in sources]
    for path, digest in build['runtime_dylib_sha256'].items():
        assert sha(path) == digest
    assert sha(ENGINE/'generic_short.cpp') == build['source_sha256']
    assert sha(ENGINE/'scratch/generic_short') == build['binary_sha256']
    reg = {
        'status': 'PREPARED_NOT_RELEASED_NO_NATIVE_EVALUATIONS',
        'T_numerator': 6, 'T_denominator': 125,
        'threshold64': '3/256000', 'threshold625': '3/250000',
        'max_native_attempts': 16, 'max_seconds': 1200,
        'per_seed_new_native_cap': 8,
        'max_optimizer_callback_requests': 80,
        'max_new_nonzero_repair_attempts_per_proposal': 6,
        'repair_k': list(PERTURBATION_K), 'repair_gap_weights': list(GAP_WEIGHTS),
        'repair_rule': 'largest exact gap j (first index on ties); add k*11^i/D to other gaps, subtract their sum from j',
        'proposal_denominator': DEN, 'ordered_min_gap': '1/100000000',
        'frame': frame, 'frame_sha256': frame_hash,
        'native_source_path': str(ENGINE/'generic_short.cpp'),
        'native_binary_path': str(ENGINE/'scratch/generic_short'),
        'native_report_expected_path': str(ENGINE/'report.md'),
        'fan_path': str(FAN), 'seeds': seeds,
        'zero_form_count': len(forms), 'zero_forms_sha256_canonical_json': hashlib.sha256(forms_bytes).hexdigest(),
        'base_rays': [{'q': q, 'ray': list(r)} for q, r in rays],
        'precision_bits': 512,
        'inputs': inputs,
        'dependencies': {'python': sys.version, 'scipy': scipy.__version__, 'numpy': numpy.__version__,
                         'python_flint': flint.__version__, 'platform': platform.platform()},
        'native_evaluations_launched': 0,
        'preparation_elapsed_seconds': time()-tic,
        'gate': 'requires main explicit release, frozen C014 report SHA, and matching source/binary/runtime hashes',
        'scientific_scope': 'bounded candidate search only, never a global minimum certificate'
    }
    dump(HERE/'preparation.json', reg)
    print(json.dumps({'status': reg['status'], 'seeds': seeds,
                      'zero_form_count': len(forms), 'frame_sha256': frame_hash,
                      'preparation_seconds': reg['preparation_elapsed_seconds']}, indent=2))

if __name__ == '__main__':
    make_preparation()
