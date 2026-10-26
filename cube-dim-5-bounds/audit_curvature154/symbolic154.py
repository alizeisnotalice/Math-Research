"""Prepared new symbolic identities only; no physical/kernel/corner cases."""
import argparse
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import time

HERE = Path(__file__).resolve().parent


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--reviewed-registration-sha', required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    registration = HERE/'symbolic_registration.json'
    assert sha(registration) == args.reviewed_registration_sha
    reg = json.loads(registration.read_bytes())
    assert sha(__file__) == reg['code_sha256']
    for path, h in reg['source_locks'].items():
        assert sha(path) == h
    assert args.output.is_absolute() and args.output.is_relative_to(HERE) and not args.output.exists()
    args.output.mkdir()
    with (args.output/'registration.json').open('x') as handle:
        json.dump({'pid': os.getpid(), 'utc': datetime.now(timezone.utc).isoformat(),
                   'reviewed_registration_sha256': args.reviewed_registration_sha,
                   'symbolic_groups': 5, 'physical_cases': 0, 'new_H_log_corners': 0}, handle)
    start = time.monotonic()
    import sympy as S
    t, eps = S.symbols('t eps', positive=True)
    r = S.symbols('r', real=True)
    H0, H1, L, U = S.symbols('H0 H1 L U', positive=True)
    integral = lambda power, a, b: S.integrate(t**power, (t, a, b))
    global_mean = S.cancel(integral(-7, eps, S.oo)/integral(-5, eps, S.oo))
    early_mean = S.cancel(integral(-7, eps, 2*eps)/integral(-5, eps, 2*eps))
    tail_mean = S.cancel(integral(-7, 2*eps, S.oo)/integral(-5, 2*eps, S.oo))
    assert global_mean == S.Rational(2, 3)/eps**2
    assert early_mean == S.Rational(7, 10)/eps**2
    assert tail_mean == S.Rational(1, 6)/eps**2
    split = S.expand(eps**2*((1-r)*early_mean+r*tail_mean))
    assert split == S.Rational(7, 10)-S.Rational(8, 15)*r
    residual = S.Rational(8, 15)*(H1/H0-L/U)
    decomposition = S.Rational(8, 15)*(U*(H1-L)+L*(U-H0))/(H0*U)
    assert S.cancel(residual-decomposition) == 0
    assert S.solve(split-S.Rational(2, 3), r) == [S.Rational(1, 16)]
    assert split.subs(r, 1) == S.Rational(1, 6)
    assert S.expand(1-S.Rational(3, 4)*r-split) == (18-13*r)/60
    assert S.solve(1024*split-256, r) == [S.Rational(27, 32)]
    assert S.Rational(1, 16)*64**2*S.Rational(2, 3) == S.Rational(512, 3)
    assert S.cancel(integral(-5, 2*eps, S.oo)/integral(-5, eps, S.oo)) == S.Rational(1, 16)
    for path, h in reg['source_locks'].items():
        assert sha(path) == h
    result = {'status': 'new_symbolic_identities_passed', 'pid': os.getpid(),
              'seconds': time.monotonic()-start, 'sympy_version': S.__version__,
              'symbolic_groups': 5, 'physical_cases': 0, 'new_H_log_corners': 0,
              'global_weighted_mean': str(global_mean), 'early_weighted_mean': str(early_mean),
              'tail_weighted_mean': str(tail_mean), 'split_normalized': str(split),
              'ratio_endpoint_residual': str(decomposition), 'switch_ratio': '1/16',
              'EPS128_q256_threshold_ratio': '27/32', 'EPS64_cap': '512/3',
              'all_sources_unchanged': True}
    with (args.output/'results.json').open('x') as handle:
        json.dump(result, handle, sort_keys=True, indent=2)
        handle.write('\n')


if __name__ == '__main__':
    main()
