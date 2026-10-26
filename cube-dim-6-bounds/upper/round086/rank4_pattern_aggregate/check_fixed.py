"""Exactly four new registered symbolic configurations; no kernel evaluation."""
import signal
import time
WALL_START = time.monotonic()
CPU_START = time.process_time()
signal.signal(signal.SIGALRM, signal.SIG_DFL)
signal.alarm(170)
import resource
cpu_hard = resource.getrlimit(resource.RLIMIT_CPU)[1]
CPU_LIMIT = min(25, cpu_hard) if cpu_hard != resource.RLIM_INFINITY else 25
resource.setrlimit(resource.RLIMIT_CPU, (CPU_LIMIT, CPU_LIMIT))
import hashlib
import json
import os
from pathlib import Path
import sys
import traceback
from itertools import combinations, product
sys.dont_write_bytecode = True

P = Path(__file__).resolve().parent
REG = json.loads((P / 'registration.json').read_text())
INPUT = json.loads((P / 'input_manifest.json').read_text())
ROWS = []
CURRENT = None


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def serial(value):
    if isinstance(value, dict):
        return {str(k): serial(v) for k, v in value.items()}
    if isinstance(value, (tuple, list)):
        return [serial(v) for v in value]
    if isinstance(value, (str, int, float, bool)) or value is None:
        return value
    return str(value)


def write_new(path, data):
    with path.open('x') as out:
        json.dump(serial(data), out, indent=2, sort_keys=True)
        out.write('\n')


def pins():
    assert sha(P / 'registration.json') == INPUT['registration_sha256']
    for path, digest in INPUT['readonly_sha256'].items():
        assert sha(Path(path)) == digest, path
    for name, digest in json.loads((P / 'pre_run_manifest.json').read_text())['sha256'].items():
        assert sha(P / name) == digest, name


def main():
    global CURRENT
    pins()
    write_new(P / 'run_started.json', dict(pid=os.getpid(), python=sys.executable,
        registration_sha256=sha(P / 'registration.json'), pre_run_manifest_sha256=sha(P / 'pre_run_manifest.json'),
        checker_CPU_hard_limit=CPU_LIMIT, default_alarm_seconds=170, automatic_retries=0))
    code, error = 0, None
    try:
        import sympy as s
        import aggregate as g

        CURRENT = 'C1_assignment_and_fixed_completion'
        a = g.matrix(REG['A_rows'])
        records = g.independent_four_bases(a)
        groups = {record['J']: [] for record in records}
        selected_count = 0
        for bits in product((0, 1), repeat=8):
            large = tuple(i for i, bit in enumerate(bits) if bit)
            candidates = [r['J'] for r in records if set(r['J']).issubset(large)]
            high_rank = g.rank(g.selected(a, large)) >= 4
            assert bool(candidates) == high_rank
            if candidates:
                first = min(candidates)
                groups[first].append(large)
                selected_count += 1
        for r in records:
            n, inv = s.Matrix(r['N']), s.Matrix(r['inverse'])
            assert n.det() == r['determinant'] != 0 and n * inv == s.eye(6)
            assert tuple(r['J']) == tuple(sorted(r['J'])) and not set(r['J']) & set(r['K'])
        assert len(set(large for vals in groups.values() for large in vals)) == selected_count
        ROWS.append(dict(id=CURRENT, status='PASS', toy_pattern_count=256,
            independent_four_count=len(records), rank_at_least_four_patterns=selected_count,
            assigned_group_sizes={','.join(map(str, j)): len(vals) for j, vals in groups.items()},
            zero_column_index=7, parallel_column_indices=[0, 6], determinant_values=sorted(set(r['determinant'] for r in records))))
        print(CURRENT, 'PASS', flush=True)

        CURRENT = 'C2_positive_factorization'
        # Abstract nonnegative derivative envelopes, not derivatives of events.
        small = {(i, r): s.Symbol(f's{i}_{r}', nonnegative=True) for i in range(4, 8) for r in range(3)}
        large = {(i, r): s.Symbol(f'l{i}_{r}', nonnegative=True) for i in range(8) for r in range(3)}
        allocations = [alpha for alpha in product(range(3), repeat=6) if sum(alpha) <= 2]
        assert len(allocations) == 28
        for alpha in allocations:
            fixed = s.prod(large[i, alpha[i]] for i in range(4))
            summed = 0
            for bits in product((0, 1), repeat=4):
                term = fixed
                for offset, i in enumerate(range(4, 8)):
                    r = alpha[i] if i < 6 else 0
                    term *= large[i, r] if bits[offset] else small[i, r]
                summed += term
            factored = fixed * s.prod(small[i, alpha[i]] + large[i, alpha[i]] for i in (4, 5))
            factored *= s.prod(small[i, 0] + large[i, 0] for i in (6, 7))
            assert s.expand(summed - factored) == 0
        ROWS.append(dict(id=CURRENT, status='PASS', derivative_allocations=28,
            enlarged_patterns_per_allocation=16, selected_tail_factors=4,
            interpretation='Order zero sums become c_t exactly; positive absolute derivative sums remain on completion columns.'))
        print(CURRENT, 'PASS', flush=True)

        CURRENT = 'C3_fixed_basis_chain_rule'
        n = s.Matrix(REG['A_rows'])[:, :6]
        inv = n.inv()
        u, v = map(s.Matrix, REG['directions'])
        z, w = inv * u, inv * v
        e, f = s.symbols('epsilon delta')
        ell = s.symbols('ell0:6')
        q = s.symbols('q0:6')
        jets = [1 + ell[i]*(e*z[i] + f*w[i]) + q[i]*(e*z[i]+f*w[i])**2/2 for i in range(6)]
        jet_product = s.prod(jets)
        first = s.diff(jet_product, e).subs({e: 0, f: 0})
        mixed = s.diff(jet_product, e, f).subs({e: 0, f: 0})
        expected_first = sum(ell[i]*z[i] for i in range(6))
        expected_mixed = sum(q[i]*z[i]*w[i] for i in range(6))
        expected_mixed += sum(ell[i]*ell[j]*z[i]*w[j] for i in range(6) for j in range(6) if i != j)
        assert s.expand(first - expected_first) == 0
        assert s.expand(mixed - expected_mixed) == 0
        ROWS.append(dict(id=CURRENT, status='PASS', determinant=n.det(), inverse=inv.tolist(),
            physical_directions=REG['directions'], solved_directions=[list(z), list(w)],
            proof_scope='Fixed affine solved variables only; no selector differentiation.'))
        print(CURRENT, 'PASS', flush=True)

        CURRENT = 'C4_time_powers_and_PSD_majorant'
        c = g.scalar_constants(REG['R'], REG['T'])
        a1, a2, b1, b2, t = [s.Rational(c[name]) for name in ('a1', 'a2', 'b1', 'b2', 'T')]
        dt, dc = s.Rational(c['d_tail']), s.Rational(c['d_completion'])
        x = s.symbols('x0:6', nonnegative=True)
        tail_pairs = list(combinations(range(4), 2))
        cross_pairs = [(i, j) for i in range(4) for j in (4, 5)]
        exact_integrated = t**3/3*(a2*sum(x[i]**2 for i in range(4)) + 2*a1**2*sum(x[i]*x[j] for i, j in tail_pairs))
        exact_integrated += a1*b1*t**2*sum(x[i]*x[j] for i, j in cross_pairs)
        exact_integrated += t*(b2*(x[4]**2+x[5]**2)+2*b1**2*x[4]*x[5])
        diagonal = dt*sum(x[i]**2 for i in range(4)) + dc*(x[4]**2+x[5]**2)
        squares = a1**2*t**3/3*sum((x[i]-x[j])**2 for i, j in tail_pairs)
        squares += a1*b1*t**2/2*sum((x[i]-x[j])**2 for i, j in cross_pairs)
        squares += b1**2*t*(x[4]-x[5])**2
        assert s.expand(diagonal - exact_integrated - squares) == 0
        assert min(a1, a2, b1, b2, dt, dc, t) > 0
        for alpha in allocations:
            k, qtail = sum(alpha), sum(alpha[:4])
            exponent = 4 - 2 - sum(alpha[4:])
            assert exponent == 2 - k + qtail and exponent >= 0
        weighted = inv.T * s.diag(*([dt]*4 + [dc]*2)) * inv
        assert s.simplify(n.T * weighted * n - s.diag(*([dt]*4 + [dc]*2))) == s.zeros(6)
        certificate = g.majorant_from_bases(records, REG['R'], REG['T'])
        bound = s.Matrix(certificate['pi6_M4'])
        rebuilt = s.zeros(6)
        for r in records:
            b = s.Matrix(r['inverse'])
            rebuilt += b.T * s.diag(*([dt]*4 + [dc]*2)) * b / (s.Rational(REG['R'])**8 * abs(r['determinant']))
        assert bound == rebuilt and bound == bound.T
        ROWS.append(dict(id=CURRENT, status='PASS', constants=c,
            positive_square_terms=len(tail_pairs)+len(cross_pairs)+1,
            time_power_allocations_checked=28, toy_pi6_trace_M4=certificate['pi6_trace_M4'],
            toy_pi6_K012=certificate['pi6_K012'], toy_pi6_M4=certificate['pi6_M4'],
            actual_body32_constant_instantiations=0))
        print(CURRENT, 'PASS', flush=True)
        pins()
    except BaseException as exc:
        error = dict(configuration=CURRENT, error=repr(exc), traceback=traceback.format_exc())
        code = 1
    result = dict(schema='A086-rank4-aggregate-results-v1',
        status='PASS_FOUR_FIXED_SYMBOLIC_CONFIGURATIONS' if code == 0 else 'FAIL_STOP_NO_RETRY',
        results=ROWS, actual_started_configurations=len(ROWS)+(1 if error else 0),
        passed_configurations=len(ROWS), failure=error,
        CPU_seconds=time.process_time()-CPU_START, wall_seconds=time.monotonic()-WALL_START,
        body_kernel_calls=0, full6D_density_evaluations=0, actual_body32_constant_instantiations=0,
        old076_configuration_reruns=0, automatic_retries=0, MATLAB_calls=0, recursive_agents=0)
    write_new(P / 'results.json', result)
    print(json.dumps(serial(result), sort_keys=True), flush=True)
    return code


if __name__ == '__main__':
    sys.exit(main())
