#!/usr/bin/env python3
"""Read-only exact arithmetic: saved N/D/P, complete roots, dyadics, payments.

No Arb evaluation, PF integral, full tuple traversal, native call or root search.
"""
import sys
sys.dont_write_bytecode = True
sys.set_int_max_str_digits(0)
import hashlib, json, math, os, time
from fractions import Fraction as Q
from pathlib import Path
from flint import fmpq, fmpq_poly, fmpz_poly

HERE = Path(__file__).resolve().parent

def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def read(path):
    return json.loads(Path(path).read_text())

def rational(value):
    value = Q(value)
    return fmpq(value.numerator, value.denominator)

def bounds(row):
    lo = Q(int(row['lower_mantissa'])) * Q(2) ** int(row['lower_binary_exponent'])
    hi = Q(int(row['upper_mantissa'])) * Q(2) ** int(row['upper_binary_exponent'])
    assert lo <= hi
    return lo, hi

def signs(values):
    result = []
    for value in values:
        if value:
            result.append(1 if value > 0 else -1)
    return sum(a != b for a, b in zip(result, result[1:]))

def verify():
    tic = time.perf_counter()
    reg = read(HERE / 'registration.json')
    for path, digest in reg['source_sha256'].items():
        assert sha(path) == digest, ('changed frozen input', path)
    for name, digest in reg['local_code_sha256'].items():
        assert sha(HERE / name) == digest
    manifest = HERE / 'manifest.json'
    if manifest.exists():
        for name, digest in read(manifest)['whitelist_sha256'].items():
            assert sha(HERE / name) == digest, ('changed final output', name)
    constants = read(reg['constants_path'])
    roots = read(reg['roots_path'])
    base = Path(reg['constants_path']).parent
    kernel = read(base / 'kernel_polynomials.json')
    critical = read(base / 'critical_polynomial.json')
    table = read(reg['partial_table_path'])
    old = read(Path(reg['partial_table_path']).parent / 'results.json')
    assert sha(base / 'kernel_polynomials.json') == critical['kernel_polynomials_sha256']
    assert sha(base / 'critical_polynomial.json') == roots['critical_polynomial_sha256']
    assert sha(reg['partial_table_path']) == kernel['partial_table_sha256']
    ell, upper_pi, c = map(Q, (constants['ell'], constants['upper_pi'], constants['c']))
    assert (ell, upper_pi, c) == (Q(3141592653,10**9), Q(3141592654,10**9), Q(1,4100))
    pl, pu = bounds(constants['pi'])
    assert ell < pl <= pu < upper_pi
    gamma = Q(old['Gamma_G'])
    safe_G = [gamma/(128*pu**5), gamma/(128*pl**5)]
    assert safe_G[0] > 0 and safe_G == list(map(Q, constants['safe_G_rational_bounds_from_pi']))
    # Rebuild the positive denominator and recover every pole coefficient from N/D.
    factors = [(rational(row['s']), row['exponent']) for row in kernel['D_factors']]
    assert dict(factors) == {rational(row['b_squared']):row['exponent'] for row in table['denominator_factors']}
    assert len(factors) == 359 and sum(m for s,m in factors) == 417
    D = fmpq_poly([1])
    for s,m in factors:
        assert s > 0 and 1 <= m <= 3
        D *= fmpq_poly([s,1]) ** m
    ndata = kernel['N']
    denominator = int(ndata['positive_common_denominator'])
    assert denominator > 0
    N = fmpq_poly([fmpq(int(v),denominator) for v in ndata['integer_coefficients_ascending']])
    assert N.degree() == 414 and D.degree() == 417 and N[0] == 0
    assert N.gcd(D).degree() == 0 and D.leading_coefficient() == 1
    assert N.leading_coefficient() == rational(old['full_weight_checksum'])
    expected = {(rational(r['b_squared']),r['order']):rational(r['coefficient']) for r in table['real_z_principal']}
    assert len(expected) == 417
    derivN = [N, N.derivative(), N.derivative().derivative()]
    recovered = 0
    for s,m in factors:
        E, remainder = divmod(D, fmpq_poly([s,1])**m)
        assert not remainder
        ed = [E]
        for k in range(1,m):
            ed.append(ed[-1].derivative())
        ec = [ed[k](-s)/math.factorial(k) for k in range(m)]
        nc = [derivN[k](-s)/math.factorial(k) for k in range(m)]
        assert ec[0] != 0
        series = []
        for k in range(m):
            value = (nc[k] - sum((ec[h]*series[k-h] for h in range(1,k+1)),fmpq(0)))/ec[0]
            series.append(value)
            assert value == expected.get((s,m-k),fmpq(0))
            recovered += 1
    assert recovered == 417
    P = N - 64*rational(ell)**6*rational(c)*D
    Z = fmpz_poly([int(v) for v in critical['primitive_integer_P_ascending']])
    scale = rational(critical['positive_scale_P_equals_scale_times_primitive'])
    assert scale > 0 and fmpq_poly(Z)*scale == P and int(Z.content()) == 1
    assert Z.degree() == 417 and Z[0] < 0 and Z.leading_coefficient() < 0
    S = fmpz_poly([int(v) for v in critical['squarefree_part_ascending']])
    assert critical['squarefree_content'] == '-1' and S == -Z
    sf = critical['squarefree_factors']
    assert len(sf) == 1 and sf[0]['multiplicity'] == 1
    assert S == fmpz_poly([int(v) for v in sf[0]['integer_coefficients_ascending']])
    assert S.gcd(S.derivative()).degree() == 0 and S[0] != 0
    assert signs(reversed(S.coeffs())) == roots['Descartes_variations'] == 2
    assert roots['completeness_proof'] == 'Descartes_upper_equals_exact_disjoint_bracket_lower'
    assert roots['positive_root_count'] == roots['positive_roots_with_multiplicity'] == len(roots['roots']) == 2
    lastz = lastt = Q(0)
    sq = fmpq_poly(S)
    for index, row in enumerate(roots['roots']):
        za,zb,ta,tb = map(Q,(row['z_lower'],row['z_upper'],row['T_lower'],row['T_upper']))
        assert 0 < lastz < za if index else 0 < za
        assert za < zb and sq(rational(za))*sq(rational(zb)) < 0
        assert ta > lastt and ta*ta <= za < zb <= tb*tb
        assert row['multiplicity'] == 1 and row['squarefree_factor_index'] == 0
        lastz,lastt = zb,tb
    # All zero even channels disappear only after the frozen full symbolic sum.
    for row in table['zero_pole_channels']:
        if row['order'] % 2 == 0:
            assert Q(row['A']) == 0
    bvalues = set()
    for s,order in expected:
        sqrat = Q(str(s))
        n,d = math.isqrt(sqrat.numerator),math.isqrt(sqrat.denominator)
        assert n*n == sqrat.numerator and d*d == sqrat.denominator
        assert Q(n,d) > 0
        bvalues.add(Q(n,d))
    assert len(bvalues) == 359
    # Exact arithmetic on complete saved dyadics; this does not reevaluate any K.
    out, receipt = read(HERE/'results.json'), read(HERE/'run_receipt.json')
    assert receipt['status'] == 'COMPLETED' and receipt['exit_code'] == 0
    assert receipt['positive_time_PF_values'] == 2 and receipt['elapsed_seconds'] < 1200
    assert receipt['registration_sha256'] == sha(HERE/'registration.json')
    assert receipt['results_sha256'] == sha(HERE/'results.json')
    assert all(receipt[key] == 0 for key in ('new_native','new_H','new_G_full','extra_positive_times','child_processes'))
    assert out['status'] == 'PASS_TWO_LOCKED_CRITICAL_PAYMENTS' and out['precision_bits'] == 640
    assert len(out['evaluations']) == 2 and out['point'] == reg['point']
    rho = (upper_pi/ell)**6
    assert out['rho_upper'] == str(rho)
    payments = []
    for index,(record,root) in enumerate(zip(out['evaluations'],roots['roots']),1):
        assert record['root_index'] == index
        assert record['T_lower'] == root['T_lower'] and record['T_upper'] == root['T_upper']
        lo,hi = bounds(record['K_at_T_upper'])
        assert 0 < lo <= hi < safe_G[1]
        payment = c*Q(root['T_lower']) - rho*hi
        assert payment == Q(record['payment_lower_exact']) > 0
        assert 0 < Q(record['payment_simple_lower']) <= payment
        assert record['payment_strictly_positive']
        payments.append(record['payment_simple_lower'])
    assert c < Q(1,4096) and Q(1,64)/c == Q(1025,16)
    return dict(status='PASS_COMPLETE_SAVED_ALL_TIME_CERTIFICATE',pid=os.getpid(),
                elapsed_seconds=time.perf_counter()-tic,principal_slots_recovered=417,
                N_degree=414,D_degree=417,complete_positive_root_count=2,simple_positive_roots=True,
                P_sign_sequence=['negative','positive','negative'],
                F_derivative_sign_sequence=['positive','negative','positive'],
                exact_positive_payment_lower_bounds=payments,Arb_or_kernel_evaluations=0,
                root_isolation_or_search=0,new_C6_bound=False)

if __name__ == '__main__':
    assert sys.argv[1:] in ([],['--write-receipt'])
    result = verify()
    if sys.argv[1:] == ['--write-receipt']:
        assert not (HERE/'checked_results.json').exists()
        (HERE/'checked_results.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result),flush=True)
