#!/usr/bin/env python3
"""New136 quartic cancellation only; no physical endpoint calculations."""
from datetime import datetime, timezone
from pathlib import Path
import hashlib
import json
import os
import time
import sympy as s

HERE = Path(__file__).resolve().parent

def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def write(name, data):
    with (HERE / name).open('x') as stream:
        json.dump(data, stream, ensure_ascii=False, indent=2, sort_keys=True)
        stream.write('\n')

def main():
    assert not (HERE / 'symbolic_quartic_results.json').exists()
    registration = {
        'task': 'audit_upper136_new_quartic_identity',
        'registered_utc': datetime.now(timezone.utc).isoformat(),
        'pid': os.getpid(), 'code_sha256': sha(__file__),
        'finite_objects': [
            'Jensen cap fourth-power cancellation to an affine Cbar,M expression',
            'positive Cbar and M coefficients of the fourth-power ceiling',
        ],
        'physical_endpoint_cases': 0, 'old_toys_repeated': 0,
        'new_H_R_gradient_log_exp_Gamma_calls': 0,
        'old_baseline_recomputations': 0,
    }
    write('symbolic_quartic_registration.json', registration)
    begin = time.monotonic()
    C, M, eps, P = s.symbols('C M eps P', positive=True)
    D = C / (2 * eps**2)
    from_cap = (5 * eps**2 / (6 * P))**4 * (4 * D + 3 * M) / 7
    affine = 625 * eps**6 * (2 * C + 3 * M * eps**2) / (9072 * P**4)
    assert s.simplify(from_cap - affine) == 0
    assert s.diff(affine, C) == 625 * eps**6 / (4536 * P**4)
    assert s.diff(affine, M) == 625 * eps**8 / (3024 * P**4)
    write('symbolic_quartic_results.json', {
        'task': registration['task'], 'status': 'PASS_two_new_symbolic_objects',
        'pid': os.getpid(), 'sympy_version': s.__version__,
        'seconds': time.monotonic() - begin,
        'registration_sha256': sha(HERE / 'symbolic_quartic_registration.json'),
        'cap_fourth_power': str(affine),
        'derivative_Cbar': str(s.diff(affine, C)),
        'derivative_M': str(s.diff(affine, M)),
        'physical_endpoint_cases': 0, 'old_toys_repeated': 0,
        'new_H_R_gradient_log_exp_Gamma_calls': 0,
        'old_baseline_recomputations': 0,
    })
    print('PASS: two new quartic symbolic objects; no physical endpoints')

if __name__ == '__main__':
    main()
