#!/usr/bin/env python3
"""Primary exact interface checks. No new union or Green evaluations."""
import hashlib
import json
import sys
import time
from itertools import product
from pathlib import Path

import sympy as s


def main():
    started = time.perf_counter()
    root = Path(__file__).resolve().parents[1]
    directions = [s.Matrix([1, *tail]) for tail in product((-1, 1), repeat=6)]
    gram = sum((v*v.T for v in directions), s.zeros(7))
    assert gram == 64*s.eye(7)
    hadamard = s.Matrix([[1]])
    for _ in range(3):
        hadamard = hadamard.row_join(hadamard).col_join(
            hadamard.row_join(-hadamard))
    basis = hadamard[1:, 1:]
    for j in range(7):
        if basis[0, j] == -1:
            basis[:, j] = -basis[:, j]
    assert abs(basis.det()) == 512
    assert all(basis[:, j] in directions for j in range(7))
    a, c, g, J, I = s.symbols('a c g J I', positive=True)
    q = s.Rational(7, 6)
    radius = (a/c)**s.Rational(1, 6)
    radial = a*radius-c*radius**7/7
    assert s.simplify(radial-s.Rational(6, 7)*a**q/c**s.Rational(1, 6)) == 0
    mass = s.Rational(6, 7)*J*c**(-s.Rational(1, 6))/(128*(g-c))
    assert s.simplify(s.diff(mass, c)/mass-(7*c-g)/(6*c*(g-c))) == 0
    optimum = s.simplify(mass.subs({c: g/7, J: 896*I}))
    assert s.simplify(optimum-7**q*I/g**q) == 0
    assert s.simplify(optimum**6-7**7*I**6/g**7) == 0
    C = s.gamma(3)/(4*s.pi**s.Rational(7, 2))
    area = 2*s.pi**s.Rational(7, 2)/s.gamma(s.Rational(7, 2))
    gaussian_coefficient = C**q*area/896
    expected = 1/(105*2**s.Rational(25, 6)*s.pi**s.Rational(13, 12))
    assert s.simplify(gaussian_coefficient/expected) == 1
    short_coefficient = 2**(-s.Rational(36, 7))*s.Rational(7, 8)*(512/s.pi**7)**s.Rational(1, 7)*(64/s.pi)**s.Rational(6, 7)
    assert s.simplify(short_coefficient/(s.Rational(7, 8)*2**s.Rational(9, 7)*s.pi**(-s.Rational(13, 7)))) == 1
    bindings = {}
    for rel in ('lower/round001', 'upper/round001/theory', 'baseline/round001'):
        directory = root/rel
        manifest = json.loads((directory/'manifest.json').read_text())
        if 'files' in manifest:
            entries = [(name, value if isinstance(value, str) else value['sha256'])
                       for name, value in manifest['files'].items()]
        else:
            entries = [(item['path'], item['sha256']) for item in manifest['artifacts']]
        assert entries, directory
        for name, digest in entries:
            path = directory/name
            assert hashlib.sha256(path.read_bytes()).hexdigest() == digest, path
            bindings[str(path.relative_to(root))] = digest
    print(json.dumps({
        'status': 'PASS exact primary interfaces and integrated file bindings',
        'python': sys.version, 'sympy': s.__version__,
        'directions': 64, 'gram': '64*Id7', 'basis_abs_det': 512,
        'exact_identity_checks': 7,
        'cost_sixth_power': '7^7 I^6 / g^7',
        'gaussian_I_coefficient': str(expected),
        'short_time_norm_coefficient': str(s.simplify(short_coefficient)),
        'file_bindings': bindings,
        'new_union_runs': 0, 'new_Green_evaluations': 0,
        'seconds': time.perf_counter()-started,
    }, ensure_ascii=False, indent=2))


if __name__ == '__main__':
    main()
