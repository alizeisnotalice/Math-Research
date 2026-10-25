#!/usr/bin/env python3
"""One beta=1 exact frame and one dependency example; no heat kernel."""
import sys
sys.dont_write_bytecode = True
from pathlib import Path
import hashlib, itertools, json, os, time
import sympy as sp

start = time.monotonic()
here = Path(__file__).resolve().parent
body = [sp.Matrix((1,) + s) for s in itertools.product((-1, 1), repeat=5)]
axes = [sp.eye(6)[:, i] for i in range(6)]
A0 = sp.Matrix.hstack(*body)
A1 = sp.Matrix.hstack(*body, *axes)
assert A0.shape == (6, 32) and A1.shape == (6, 38)
assert A0*A0.T == 32*sp.eye(6)
assert A1*A1.T == 33*sp.eye(6)
assert A1.rank() == 6
assert sp.Matrix.hstack(*axes).det() == 1
v = sp.Matrix([1]*6)
w = sp.Matrix([1]*5 + [-1])
e = axes[5]
assert v in body and w in body
assert v-w == 2*e
assert sp.Matrix.hstack(v, w, e).rank() == 2

out = {
    'status': 'PASS_BETA1_FRAME_AND_DEPENDENT_TRIPLE',
    'beta_exact': '1',
    'dimension': 6, 'body_columns': 32, 'augmented_columns': 38,
    'rank': 6,
    'body_Gram_diagonal': [32]*6,
    'augmented_Gram_diagonal': [33]*6,
    'all_off_diagonal_Gram_entries_zero': True,
    'coordinate_basis_determinant': 1,
    'dependent_triple': {'v': list(map(int,v)), 'w': list(map(int,w)),
                         'e6': list(map(int,e)),
                         'identity': 'v-w=2e6', 'rank': 2},
    'independence_of_arbitrary_three_columns_is_false': True,
    'symbol_definition': 'sum_body |v dot xi| + sum_i |xi_i|',
    'general_beta_Gram_proved_in_report': '(32+beta^2)I6',
    'new_kernel_evaluations': 0, 'MC_samples': 0, 'matlab_calls': 0,
    'fan_replayed': False, 'minimum_certified': False, 'new_C6_bound': False,
    'pid': os.getpid(), 'python': sys.version, 'sympy': sp.__version__,
    'elapsed_seconds': time.monotonic()-start,
    'script_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    'exit_code': 0
}
(here/'structure_results.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
