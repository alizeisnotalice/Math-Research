#!/usr/bin/env python3
"""Three-orbit rank-two diagnostic; no random samples or spatial grid."""
import hashlib
import itertools
import json
import math
import os
import sys
import time
from pathlib import Path
import numpy as np
import scipy
from scipy.integrate import quad
from scipy.special import gamma
import sympy as sp

ROOT = Path(__file__).resolve().parent
START = time.monotonic()
V = sp.Matrix([(1,) + signs for signs in itertools.product((-1, 1), repeat=4)])
TRIPLES = ((0, 1, 2), (0, 1, 6), (0, 3, 5))
WEIGHTS = (3/16, 3/8, 7/16)
BETA = 3/4
B = 3/20
RHO = 247/256

def digest(path):
    path = Path(path)
    return {'path': str(path), 'sha256': hashlib.sha256(path.read_bytes()).hexdigest()}

def dump(name, obj):
    tmp = ROOT / (name + '.tmp')
    tmp.write_text(json.dumps(obj, ensure_ascii=False, indent=2) + '\n')
    tmp.replace(ROOT / name)

def geometry(triple):
    selected = V[list(triple), :]
    ns = selected.nullspace()
    q1 = ns[0] / sp.sqrt(ns[0].dot(ns[0]))
    q2raw = ns[1] - q1.dot(ns[1]) * q1
    q2 = sp.simplify(q2raw / sp.sqrt(q2raw.dot(q2raw)))
    Q = q1.row_join(q2)
    W = sp.simplify(V * Q)
    assert sp.simplify(Q.T * Q) == sp.eye(2)
    assert sp.simplify(W.T * W) == 16 * sp.eye(2)
    return W, int((selected * selected.T).det())

def one(triple):
    exact, gram = geometry(triple)
    w = np.array(exact, dtype=float)
    roots = [0., 2 * math.pi]
    for x, y in w:
        if abs(x) + abs(y) > 1e-14:
            root = math.atan2(-x, y) % math.pi
            roots.extend((root, root + math.pi))
    roots.sort()
    breaks = []
    for root in roots:
        if not breaks or root - breaks[-1] > 1e-12:
            breaks.append(root)
    def integrand(theta):
        e = np.array([math.cos(theta), math.sin(theta)])
        a = (w @ e) ** 2 / 16
        return float(np.sum(a ** BETA) ** (-1 / BETA))
    contributions = [quad(integrand, a, z, epsabs=1e-12, epsrel=1e-12, limit=200)
                     for a, z in zip(breaks[:-1], breaks[1:])]
    angular = sum(x[0] for x in contributions) / (2 * math.pi)
    error = sum(x[1] for x in contributions) / (2 * math.pi)
    k = gamma(1 + 1 / BETA) * angular
    return {'triple': triple, 'gram_determinant': gram,
            'projected_rows_exact': [[str(x) for x in exact.row(i)] for i in range(16)],
            'zero_rows': [i for i in range(16) if all(x == 0 for x in exact.row(i))],
            'breakpoints_diagnostic': breaks, 'angular_mean': angular,
            'quad_error_estimate_not_certified': error,
            'K_E_det_minus_half': k, 'K_pow_2b': k ** (2 * B),
            'elapsed_seconds': time.monotonic() - START}

def main():
    registration = json.loads((ROOT / 'registration.json').read_text())
    assert registration['status'] == 'registered_before_computation'
    script_hash = digest(__file__)
    cp = ROOT / 'checkpoint.json'
    rows = []
    if cp.exists():
        prior = json.loads(cp.read_text())
        if prior.get('script_sha256') == script_hash['sha256']:
            rows = prior['representatives']
    for triple in TRIPLES[len(rows):]:
        if time.monotonic() - START > 120:
            return
        rows.append(one(triple))
        dump('checkpoint.json', {'status': 'representative_completed',
              'pid': os.getpid(), 'script_sha256': script_hash['sha256'],
              'representatives': rows, 'next_representative': len(rows)})
        print(json.dumps(rows[-1], ensure_ascii=False), flush=True)
    m = lambda s: gamma(1 + s / BETA) / gamma(1 + s)
    geometric = math.exp(sum(t * (2*B) * math.log(r['K_E_det_minus_half'])
                             for t, r in zip(WEIGHTS, rows)))
    j_plain = 2 * gamma(7/4) * m(B)**3 * geometric
    j_fractional = 2 * gamma(7/4) * m(B*RHO)**(3/RHO) * geometric
    c = 3 / (2 ** (15/2) * math.pi ** (5/2))
    i_from_j = lambda j: c**(10/7) * math.pi**2 / 60 * j**(10/7)
    result = {'status': 'diagnostic_complete_not_certificate', 'script': script_hash,
              'registration': digest(ROOT / 'registration.json'),
              'source_inputs': [digest('/Users/zhengzhihao/Desktop/T/Math-Research-cube5-upper-audit-cont/cube-dim-5-bounds/' + p)
                                for p in ('upper_I021/report.md','upper_I021/results.json','upper_I022/report.md','upper_I022/results.json')],
              'representatives': rows, 'geometric_K_factor': geometric,
              'J_plain': j_plain, 'J_fractional': j_fractional,
              'I_plain': i_from_j(j_plain), 'I_fractional': i_from_j(j_fractional),
              'runtime': {'python': sys.executable, 'scipy': scipy.__version__, 'sympy': sp.__version__,
                          'pid': os.getpid(), 'seconds': time.monotonic() - START},
              'new_C5_upper_bound': False}
    dump('diagnostic.json', result)
    print(json.dumps({key: value for key, value in result.items()
                      if key not in ('representatives', 'source_inputs')}, indent=2), flush=True)

if __name__ == '__main__':
    main()
