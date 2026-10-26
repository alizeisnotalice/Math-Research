#!/usr/bin/env python3
"""Exact-root Lobatto upper rule and one fixed 1920-element Arb orbit."""
from fractions import Fraction as F
from pathlib import Path
from itertools import product, permutations
import argparse
import hashlib
import json
import os
import sys
import time
import sympy as sp
import flint
from flint import arb, arb_mat, fmpq, fmpq_mat, ctx

BASE = Path(__file__).resolve().parent
START = time.monotonic()
parser = argparse.ArgumentParser()
parser.add_argument('--output-dir', type=Path, default=BASE)
args = parser.parse_args()
OUT = args.output_dir
OUT.mkdir(parents=True, exist_ok=True)
ctx.prec = 512


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_json(name, data):
    path = OUT/name
    tmp = path.with_suffix(path.suffix+'.tmp')
    tmp.write_text(json.dumps(data, indent=2, sort_keys=True)+'\n')
    tmp.replace(path)


def ar(x):
    x = F(x)
    return arb(x.numerator)/x.denominator


def ends(x):
    return F(str(x.lower().fmpq())), F(str(x.upper().fmpq()))


def from_ends(lo, hi):
    return ar(F(lo)).union(ar(F(hi)))


def outward_decimal(x, up=False, digits=65):
    scale = 10**digits
    n = (x*scale).__ceil__() if up else (x*scale).__floor__()
    sign = '-' if n < 0 else ''
    a, b = divmod(abs(n), scale)
    return f'{sign}{a}.{b:0{digits}d}'


def encode_ball(x):
    lo, hi = ends(x)
    return {'lower_fraction': str(lo), 'upper_fraction': str(hi),
            'decimal_outward': [outward_decimal(lo), outward_decimal(hi, True)]}


def decode_ball(row):
    return from_ends(row['lower_fraction'], row['upper_fraction'])


def sf(x):
    return F(int(x.p), int(x.q))


def fq(x):
    x = F(x)
    return fmpq(x.numerator, x.denominator)


def beta_moment(a, n):
    out = F(1)
    for k in range(n):
        out *= (a+k)/(2*a+k)
    return out


def integrate_poly(poly, a):
    return sum((sf(poly.nth(k))*beta_moment(a, k) for k in range(poly.degree()+1)), F(0))


def polynomial_ball(coefficients_descending, x):
    value = arb(0)
    for coefficient in coefficients_descending:
        value = value*x+ar(coefficient)
    return value


registration_hash = sha(BASE/'registration.json')
script_hash = sha(Path(__file__))
reg = json.loads((BASE/'registration.json').read_text())
if (OUT/'results.json').exists():
    saved = json.loads((OUT/'results.json').read_text())
    assert saved['registration_sha256'] == registration_hash and saved['script_sha256'] == script_hash
    print(json.dumps({'status': 'existing_complete_certificate_read_only', 'results_sha256': sha(OUT/'results.json')}))
    sys.exit(0)
if (OUT/'failure.json').exists():
    print('Recorded certificate failure exists; no retry performed.')
    sys.exit(2)
resume_state = json.loads((OUT/'checkpoint.json').read_text()) if (OUT/'checkpoint.json').exists() else None
if resume_state is not None:
    assert resume_state['registration_sha256'] == registration_hash and resume_state['script_sha256'] == script_hash

P049, P050 = BASE.parent/'upper_orbit049', BASE.parent/'upper_orbit050'
expected = reg['single_fixed_shape']
assert sha(P049/'manifest.json') == expected['source049_manifest_sha256']
assert sha(P050/'manifest.json') == expected['source050_manifest_sha256']
assert sha(P049/'pending.json') == expected['source_pending_sha256']
assert sha(P050/'repair99.json') == expected['source050_repair_sha256']
assert sha(P050/'results.json') == expected['source050_results_sha256']
for directory in (P049, P050):
    for item in json.loads((directory/'manifest.json').read_text())['white_list']:
        assert sha(directory/item['path']) == item['sha256'], item['path']


def checkpoint(stage, completed=0, partial=None, **extra):
    data = {'pid': os.getpid(), 'stage': stage, 'completed_groups': completed,
            'total_groups': 1920, 'elapsed_seconds_this_run': time.monotonic()-START,
            'precision_bits': ctx.prec, 'registration_sha256': registration_hash,
            'script_sha256': script_hash, 'resume': 'python -B certify.py', **extra}
    if partial is not None:
        data['partial_sum_lobatto_orbit_values'] = encode_ball(partial)
    write_json('checkpoint.json', data)


if resume_state is None or resume_state.get('completed_groups', 0) == 0:
    checkpoint('registered_start')
print(json.dumps({'pid': os.getpid(), 'precision_bits': 512,
                  'registration_sha256': registration_hash, 'script_sha256': script_hash}), flush=True)
try:
    # Exact rational orthogonal polynomials for normalized modified Beta(9/4,9/4).
    t = sp.symbols('t')
    a, q = F(9, 4), F(5, 4)
    p = [sp.Poly(1, t, domain=sp.QQ), sp.Poly(t-sp.Rational(1, 2), t, domain=sp.QQ)]
    b = {}
    h = [F(1)]
    for k in range(1, 6):
        b[k] = F(k)*(k+2*a-2)/(4*(2*k+2*a-1)*(2*k+2*a-3))
        assert b[k] > 0
        bk = sp.Rational(b[k].numerator, b[k].denominator)
        p.append(sp.Poly((t-sp.Rational(1, 2))*p[k].as_expr()-bk*p[k-1].as_expr(), t, domain=sp.QQ))
        h.append(h[-1]*b[k])
    assert len(p) == 7 and p[6].degree() == 6 and p[6].LC() == 1
    assert sp.expand(p[6].as_expr().subs(t, 1-t)-p[6].as_expr()) == 0
    orthogonality_checks = 0
    for k in range(7):
        for j in range(k):
            assert integrate_poly(sp.Poly(p[k].as_expr()*t**j, t), a) == 0
            orthogonality_checks += 1
    for k in range(6):
        assert integrate_poly(sp.Poly(p[k].as_expr()**2, t), a) == h[k]
    assert beta_moment(q, 1)-beta_moment(q, 2) == F(5, 28)
    denom, integer_poly = p[6].clear_denoms()
    _, integer_poly = integer_poly.primitive()
    isolation_eps = sp.Rational(1, 2**560)
    isolated = integer_poly.intervals(eps=isolation_eps)
    assert len(isolated) == 6
    root_intervals = []
    for bounds, multiplicity in isolated:
        lo, hi = map(sf, bounds)
        assert multiplicity == 1 and F(0) < lo < hi < F(1)
        assert hi-lo <= F(1, 2**560)
        count = integer_poly.count_roots(sp.Rational(lo.numerator, lo.denominator), sp.Rational(hi.numerator, hi.denominator))
        assert count == 1
        assert integer_poly.eval(bounds[0])*integer_poly.eval(bounds[1]) < 0
        root_intervals.append((lo, hi))
    assert integer_poly.count_roots(0, 1) == 6
    for (lo, hi), (nlo, nhi) in zip(root_intervals, root_intervals[1:]):
        assert hi < nlo
    coefficients = [[sf(x) for x in pk.all_coeffs()] for pk in p]
    nodes_inside = [from_ends(lo, hi) for lo, hi in root_intervals]
    gauss_weights, interior_weights = [], []
    for node in nodes_inside:
        residual = polynomial_ball(coefficients[6], node)
        assert ends(residual)[0] <= 0 <= ends(residual)[1]
        christoffel = sum((polynomial_ball(coefficients[k], node)**2/ar(h[k]) for k in range(6)), arb(0))
        assert christoffel > 0
        g = 1/christoffel
        w = ar(F(5, 28))*g/(node*(1-node))
        assert g > 0 and w > 0
        gauss_weights.append(g)
        interior_weights.append(w)
    p0, p1 = sf(p[6].eval(0)), sf(p[6].eval(1))
    assert p0 == p1 and p0 != 0
    end0_poly = sp.Poly((1-t)*p[6].as_expr()**2/sp.Rational(p0.numerator, p0.denominator)**2, t)
    end1_poly = sp.Poly(t*p[6].as_expr()**2/sp.Rational(p1.numerator, p1.denominator)**2, t)
    w0, w1 = integrate_poly(end0_poly, q), integrate_poly(end1_poly, q)
    assert w0 == w1 and w0 > 0
    nodes = [arb(0)]+nodes_inside+[arb(1)]
    weights = [ar(w0)]+interior_weights+[ar(w1)]
    moment_checks = []
    for k in range(14):
        moment_ball = sum((w*x**k for w, x in zip(weights, nodes)), arb(0))
        exact_moment = beta_moment(q, k)
        lo, hi = ends(moment_ball)
        assert lo <= exact_moment <= hi
        moment_checks.append({'degree': k, 'exact_beta_moment': str(exact_moment),
                              'outward_quadrature_moment': encode_ball(moment_ball),
                              'contains_exact_moment': True})
    gl, gh = ends(sum(gauss_weights, arb(0)))
    assert gl <= 1 <= gh
    roots_certificate = {
        'degree': 6, 'monic_coefficients_descending': [str(x) for x in coefficients[6]],
        'primitive_integer_coefficients_descending': [str(x) for x in integer_poly.all_coeffs()],
        'b_k': {str(k): str(v) for k, v in b.items()}, 'h_0_to_5': [str(x) for x in h],
        'exact_orthogonality_checks': orthogonality_checks, 'exact_squared_norm_checks': 6,
        'exact_root_count_in_0_1': 6, 'root_isolation_width_at_most': str(F(1, 2**560)),
        'roots': [{'lower_rational': str(lo), 'upper_rational': str(hi), 'exact_root_count': 1,
                   'arb_node': encode_ball(node), 'gauss_weight': encode_ball(g), 'lobatto_weight': encode_ball(w)}
                  for (lo, hi), node, g, w in zip(root_intervals, nodes_inside, gauss_weights, interior_weights)],
        'endpoint_weight_0_exact': str(w0), 'endpoint_weight_1_exact': str(w1),
        'normalizing_modified_mass': '5/28', 'all_weights_positive': True,
        'degree_13_exactness': 'Proved by subtracting endpoint interpolation and exact six-node Gauss integration of degree<=11 under the modified measure.',
        'directed_moment_checks_0_to_13': moment_checks, 'precision_bits': 512,
        'no_float_root_solver_used': True
    }
    write_json('roots_weights.json', roots_certificate)
    if not resume_state or resume_state.get('completed_groups', 0) == 0:
        checkpoint('exact_roots_and_positive_weights_complete')

    # One and only one exact rational shape, from the fixed original binary floats.
    pending = json.loads((P049/'pending.json').read_text())
    assert pending['sample_id'] == 99
    yy = [F(*float(y).as_integer_ratio()) for y in pending['gamma_draws']]
    assert all(y > 0 for y in yy)
    inverse_y_sum = sum((1/y for y in yy), F(0))
    pp = [(1/y)/inverse_y_sum for y in yy]
    assert sum(pp, F(0)) == 1 and all(x > 0 for x in pp)
    vv = [(1,)+e for e in product((-1, 1), repeat=4)]
    M = [[sum((weight*v[i]*v[j] for weight, v in zip(pp, vv)), F(0)) for j in range(5)] for i in range(5)]
    assert all(M[i][i] == 1 for i in range(5))
    MQ = fmpq_mat([[fq(x) for x in row] for row in M])
    d = F(str(MQ.det()))
    assert d > 0
    AQ = MQ.inv()
    AA = [[F(str(AQ[i, j])) for j in range(5)] for i in range(5)]
    assert MQ*AQ == fmpq_mat([[int(i == j) for j in range(5)] for i in range(5)])
    minors = [F(str(fmpq_mat([[fq(M[i][j]) for j in range(k)] for i in range(k)]).det())) for k in range(1, 6)]
    assert all(x > 0 for x in minors)
    repaired = json.loads((P050/'repair99.json').read_text())
    assert [str(y) for y in yy] == [str(F(s)) for s in repaired['input_exact_binary_fractions']]
    fixed_shape = {'sample_id': 99, 'exact_gamma_binary_rationals': [str(y) for y in yy],
                   'exact_P': [str(x) for x in pp], 'exact_M': [[str(x) for x in row] for row in M],
                   'exact_determinant': str(d), 'positive_leading_principal_minors': [str(x) for x in minors],
                   'exact_inverse_identity_verified': True, 'exact_diagonal_one': True,
                   'source_pending_sha256': sha(P049/'pending.json'),
                   'exact_inverse_serialization_sha256': hashlib.sha256(json.dumps([[str(x) for x in row] for row in AA]).encode()).hexdigest()}
    write_json('fixed_shape.json', fixed_shape)
    A_arb = [[ar(x) for x in row] for row in AA]
    d_arb = ar(d)
    group = [(permutation, (1,)+sign) for permutation in permutations(range(5)) for sign in product((-1, 1), repeat=4)]
    assert len(group) == len(set(group)) == 1920
    if not resume_state or resume_state.get('completed_groups', 0) == 0:
        checkpoint('fixed_exact_shape_ready')

    start_group = resume_state.get('completed_groups', 0) if resume_state else 0
    if start_group:
        assert start_group % 120 == 0 and start_group < 1920
        partial = decode_ball(resume_state['partial_sum_lobatto_orbit_values'])
        min_group_lo = F(resume_state['minimum_group_lower'])
        max_group_hi = F(resume_state['maximum_group_upper'])
        max_endpoint_error = F(resume_state['maximum_endpoint_error'])
        assert sha(OUT/'roots_weights.json') == resume_state['roots_weights_sha256']
        assert sha(OUT/'fixed_shape.json') == resume_state['fixed_shape_sha256']
        block_records = json.loads((OUT/'blocks.json').read_text())['blocks']
        assert block_records[-1]['completed_groups'] == start_group
    else:
        partial = arb(0)
        min_group_lo, max_group_hi = None, None
        max_endpoint_error = F(0)
        block_records = []
    determinant_count = 8*start_group
    for group_index in range(start_group, 1920):
        permutation, signs = group[group_index]
        B_arb = [[signs[i]*signs[j]*A_arb[permutation[i]][permutation[j]] for j in range(5)] for i in range(5)]
        value = arb(0)
        for node_index, (node, weight) in enumerate(zip(nodes, weights)):
            mixed = arb_mat([[node*A_arb[i][j]+(1-node)*B_arb[i][j] for j in range(5)] for i in range(5)])
            det_ball = mixed.det()
            assert det_ball > 0, ('determinant_enclosure_not_positive', group_index, node_index, str(det_ball))
            norm_det = d_arb*det_ball
            assert norm_det > 0
            kernel = 1/norm_det.sqrt()
            klo, khi = ends(kernel)
            assert klo > 0 and klo <= 1
            if node_index in (0, 7):
                assert klo <= 1 <= khi
                max_endpoint_error = max(max_endpoint_error, 1-klo, khi-1)
                kernel = arb(1)  # proved exact endpoint identity, not numerical clipping.
            value += weight*kernel
            determinant_count += 1
        vlo, vhi = ends(value)
        assert vlo > 0 and vlo <= 1
        min_group_lo = vlo if min_group_lo is None else min(min_group_lo, vlo)
        max_group_hi = vhi if max_group_hi is None else max(max_group_hi, vhi)
        partial += value
        if (group_index+1) % 120 == 0:
            record = {'completed_groups': group_index+1, 'determinants': determinant_count,
                      'partial_sum': encode_ball(partial), 'elapsed_seconds': time.monotonic()-START}
            block_records.append(record)
            write_json('blocks.json', {'blocks': block_records})
            checkpoint('group_blocks', group_index+1, partial,
                       roots_weights_sha256=sha(OUT/'roots_weights.json'), fixed_shape_sha256=sha(OUT/'fixed_shape.json'),
                       minimum_group_lower=str(min_group_lo), maximum_group_upper=str(max_group_hi),
                       maximum_endpoint_error=str(max_endpoint_error))
            print(json.dumps({'groups': group_index+1, 'elapsed_seconds': time.monotonic()-START}), flush=True)
            if time.monotonic()-START >= 110 and group_index+1 < 1920:
                checkpoint('time_limit_checkpoint_preserved', group_index+1, partial,
                           roots_weights_sha256=sha(OUT/'roots_weights.json'), fixed_shape_sha256=sha(OUT/'fixed_shape.json'),
                           minimum_group_lower=str(min_group_lo), maximum_group_upper=str(max_group_hi),
                           maximum_endpoint_error=str(max_endpoint_error))
                sys.exit(75)
    assert determinant_count == 15360
    r_lobatto = partial/1920
    rlo, rhi = ends(r_lobatto)
    assert 0 < rlo <= rhi < 1
    assert rlo >= w0+w1
    norm_factor = (r_lobatto.log()/5).exp()
    previous_float = F.from_float(float(repaired['orbit_r16']))
    comparison_difference = from_ends(rlo-previous_float, rhi-previous_float)
    results = {
        'status': 'one_fixed_rational_shape_full_group_lobatto_upper_expression_certified',
        'shape_count': 1, 'group_count': 1920, 'lobatto_nodes': 8, 'arb_determinants': determinant_count,
        'precision_bits': 512, 'r_lobatto_directed_enclosure': encode_ball(r_lobatto),
        'strict_mathematical_direction': '0 < true_r_Gamma(M) <= r_Lobatto(M) <= 1 for every fixed SPD M',
        'uniform_endpoint_floor_exact': str(w0+w1),
        'norm_factor_r_lobatto_to_one_fifth': encode_ball(norm_factor),
        'smallest_individual_group_lower': str(min_group_lo),
        'largest_individual_group_upper': str(max_group_hi),
        'largest_endpoint_computed_enclosure_radius_about_one': str(max_endpoint_error),
        'saved050_r16_float_comparison': {'value': float(previous_float), 'exact_binary_rational': str(previous_float),
                                         'r_lobatto_minus_saved_float': encode_ball(comparison_difference),
                                         'saved_value_is_certified_lower_bound': False, 'old_rule_rerun': False},
        'root_weights_sha256': sha(OUT/'roots_weights.json'), 'fixed_shape_sha256': sha(OUT/'fixed_shape.json'),
        'registration_sha256': registration_hash, 'script_sha256': script_hash,
        'R_expectation_certified': False, 'new_full_I_bound': False, 'new_C5_bound': False,
        'all_near_degenerate_shapes_domination': 'Universal Lobatto value<=1; prior integrable determinant moment remains a dominating function.'
    }
    write_json('results.json', results)
    # Verify source bytes once more after all fixed-shape computations.
    for directory in (P049, P050):
        for item in json.loads((directory/'manifest.json').read_text())['white_list']:
            assert sha(directory/item['path']) == item['sha256'], item['path']
    checkpoint('complete', 1920, partial)
    receipt = {'pid': os.getpid(), 'elapsed_seconds': time.monotonic()-START, 'exit_status': 0,
               'python_flint': flint.__version__, 'sympy': sp.__version__, 'precision_bits': 512,
               'results_sha256': sha(OUT/'results.json'), 'script_sha256': script_hash,
               'registration_sha256': registration_hash, 'all_source_manifest_files_unchanged': True}
    write_json('run_receipt.json', receipt)
    print(json.dumps(receipt, sort_keys=True), flush=True)
    print(json.dumps({'r_lobatto': results['r_lobatto_directed_enclosure']['decimal_outward'],
                      'endpoint_floor': str(w0+w1), 'saved_noncertified_r16': float(previous_float),
                      'certifies_R_expectation': False}), flush=True)
except Exception as exc:
    failure = {'pid': os.getpid(), 'stage': 'certificate_failure', 'exception': type(exc).__name__,
               'detail': str(exc), 'elapsed_seconds': time.monotonic()-START,
               'registration_sha256': registration_hash, 'script_sha256': script_hash}
    write_json('failure.json', failure)
    print(json.dumps(failure), flush=True)
    sys.exit(2)
