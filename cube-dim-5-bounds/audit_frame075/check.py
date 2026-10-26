#!/usr/bin/env python3
"""One fixed exact matrix qualification. No kernel, gradients or special functions."""
import time
START = time.monotonic()
from fractions import Fraction as F
from itertools import product
from pathlib import Path
import hashlib
import json
import os
import signal
import sys

HERE = Path(__file__).resolve().parent


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def transpose(A):
    return [list(x) for x in zip(*A)]


def mm(A, B):
    BT = transpose(B)
    return [[sum((x*y for x, y in zip(row, col)), F(0)) for col in BT] for row in A]


def eye(n):
    return [[F(i == j) for j in range(n)] for i in range(n)]


def determinant(A):
    A = [list(map(F, row)) for row in A]
    out = F(1)
    for j in range(len(A)):
        k = next(k for k in range(j, len(A)) if A[k][j])
        if k != j:
            A[j], A[k] = A[k], A[j]
            out = -out
        out *= A[j][j]
        for i in range(j + 1, len(A)):
            factor = A[i][j] / A[j][j]
            A[i] = [a - factor*b for a, b in zip(A[i], A[j])]
    return out


def inverse(A):
    n = len(A)
    R = [list(map(F, row)) + e for row, e in zip(A, eye(n))]
    for j in range(n):
        k = next(k for k in range(j, n) if R[k][j])
        R[j], R[k] = R[k], R[j]
        d = R[j][j]
        R[j] = [x/d for x in R[j]]
        for i in range(n):
            if i != j:
                a = R[i][j]
                R[i] = [x-a*y for x, y in zip(R[i], R[j])]
    return [row[n:] for row in R]


def encode(A):
    return [[str(x) for x in row] for row in A]


def main():
    def timeout(*_):
        raise TimeoutError('160s internal stop, below 180s total wall budget')
    signal.signal(signal.SIGALRM, timeout)
    signal.alarm(160)
    reg = json.loads((HERE/'registration.json').read_text())
    assert reg['status'] == 'registered_before_single_fixed_matrix_batch'
    assert reg['source_sha256'] == sha(__file__)
    output = HERE/'results.json'
    assert not output.exists(), 'Do not repeat completed or failed finite batch'
    for item in reg['readonly_inputs']:
        assert sha(item['path']) == item['sha256']
    out = {'pid': os.getpid(), 'python': sys.executable, 'status': 'running',
           'source_sha256': sha(__file__), 'registration_sha256': sha(HERE/'registration.json')}
    print(json.dumps(out), flush=True)
    try:
        V = [(1,) + s for s in product((-1, 1), repeat=4)]
        B = transpose(V)
        A = [[F(x, 16) for x in row] for row in V]
        BBt = mm(B, transpose(B))
        assert BBt == [[16*x for x in row] for row in eye(5)]
        assert mm(B, A) == eye(5)
        AtA = mm(transpose(A), A)
        assert AtA == [[x/16 for x in row] for row in eye(5)]
        P = mm(A, B)
        assert transpose(P) == P and mm(P, P) == P
        R = [[eye(16)[i][j] - P[i][j] for j in range(16)] for i in range(16)]
        assert mm(B, R) == [[F(0) for _ in range(16)] for _ in range(5)]
        # Each R column is e_j - B^dagger v_j, proving equality of fiber translations.
        assert set(V + [tuple(-x for x in v) for v in V]) == set(product((-1, 1), repeat=5))
        Ccolumns = [(1, 1, 1, 1, 1)]
        for j in range(1, 5):
            v = [1]*5
            v[j] = -1
            Ccolumns.append(tuple(v))
        C = transpose(Ccolumns)
        detC = determinant(C)
        assert detC == 16
        D = transpose([v for v in V if v not in Ccolumns])
        minusCiD = [[-x for x in row] for row in mm(inverse(C), D)]
        tangent = minusCiD + eye(11)
        gram = mm(transpose(tangent), tangent)
        gramdet = determinant(gram)
        detBBt = determinant(BBt)
        assert detBBt == 1024**2
        assert gramdet == detBBt / detC**2 == 64**2
        assert F(64, 1024) == 1/abs(detC)
        # Verify the sixteen former basis matrices exactly, then compare to the full-frame Q.
        Ncolumns = [(1,1,1,1,1),(1,-1,-1,-1,1),(1,-1,-1,1,-1),
                    (1,-1,1,-1,-1),(1,1,-1,-1,-1)]
        basis_rows = []
        for s in product((-1, 1), repeat=4):
            diagonal = (1,) + s
            Ns = transpose([tuple(diagonal[i]*v[i] for i in range(5)) for v in Ncolumns])
            assert all(tuple(v) in V for v in transpose(Ns))
            assert abs(determinant(Ns)) == 48
            Ni = inverse(Ns)
            Q = [row[1:] for row in mm(transpose(Ni), Ni)[1:]]
            expected = [[F(i == j, 4)-F(s[i]*s[j], 36) for j in range(4)] for i in range(4)]
            assert Q == expected
            # Q-I/16 = (11/144)I + (4I-ss^T)/36; last term PSD by Cauchy.
            for i in range(4):
                for j in range(4):
                    assert Q[i][j]-F(i == j, 16) == F(11*(i == j), 144)+F(4*(i == j)-s[i]*s[j], 36)
            basis_rows.append({'s': list(s), 'Q': encode(Q), 'minimum_gap_eigenvalue': '11/144'})
        for item in reg['readonly_inputs']:
            assert sha(item['path']) == item['sha256']
        out.update({'status': 'all_fixed_matrix_and_fiber_identities_passed',
                    'directions': V, 'B_Btranspose': encode(BBt), 'right_inverse': encode(A),
                    'right_inverse_Gram': encode(AtA), 'coarea_J_B': '1024',
                    'coarea_density_prefactor': '1/1024', 'coordinate_fiber_basis': Ccolumns,
                    'coordinate_basis_determinant': str(detC), 'fiber_parameter_Gram_determinant': str(gramdet),
                    'fiber_parameter_volume_factor': '64', 'effective_coordinate_prefactor': '1/16',
                    'all_16_coordinate_lifts_differ_by_kernel_vectors': True,
                    'signed_body_directions_are_all_32_cube_vertices': True,
                    'tail_Q_full': encode([[F(i == j, 16) for j in range(4)] for i in range(4)]),
                    'old_basis_comparisons': basis_rows, 'new_H_gradient_log_exp_Gamma_calls': 0,
                    'old_162_box_replay': False, 'new_global_certificate': False})
    except Exception as error:
        out.update({'status': 'stopped_failure_no_rerun', 'error': repr(error)})
        raise
    finally:
        out['seconds_from_process_initialization_before_imports'] = time.monotonic()-START
        output.write_text(json.dumps(out, indent=2)+'\n')
        signal.alarm(0)
    print(json.dumps({k:out[k] for k in ('status','pid','seconds_from_process_initialization_before_imports')}), flush=True)


if __name__ == '__main__':
    main()
