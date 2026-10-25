#!/usr/bin/env python3
"""One cutoff Cauchy-potential box, rational/log Taylor jets and full tail.

Uses the locked017 table only. No production code from017/018/035 is run.
The three center candidates, degree, cutoff, width and precision are fixed.
"""
from __future__ import annotations
import argparse
from fractions import Fraction as F
from itertools import product
from math import comb, factorial, lcm, prod
from pathlib import Path
import gzip
import hashlib
import json
import os
import resource
import sys
import time

import flint
from flint import arb, ctx, fmpq

if hasattr(sys, 'set_int_max_str_digits'):
    sys.set_int_max_str_digits(0)

HERE = Path(__file__).resolve().parent


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def rq(value):
    value = F(str(value))
    return fmpq(value.numerator, value.denominator)


def endpoint(value, upper=False):
    return str((value.upper() if upper else value.lower()).fmpq())


def atomic(path, value):
    temporary = path.with_suffix(path.suffix + '.tmp')
    with temporary.open('w') as output:
        json.dump(value, output, sort_keys=True, separators=(',', ':'))
        output.write('\n')
        output.flush()
        os.fsync(output.fileno())
    os.replace(temporary, path)


def simplex_jet(ids, vertices, nodes, D, a, alphas, parents):
    """Exact fmpq coefficients of the fourth DD of h_a; six symbols.

    Symbol0 is rational1, others are log(A^2+(32*k)^2) at the local
    nodes. Constant -log(32) drops from every fourth divided difference.
    """
    ids = sorted(ids, key=lambda i: nodes[i])
    k = [nodes[i] for i in ids]
    assert len(set(k)) == 5
    slopes = [tuple(D*z for z in vertices[i][1:]) for i in ids]
    a_den = int(a.denominator)
    a_num = int(a.numerator)
    log_keys = [a_num*a_num + (a_den*n)**2 for n in k]
    table = []
    for right, node in enumerate(k):
        denominator = a*a + node*node
        derivatives = {
            1: node/denominator,
            2: (a*a-node*node)/(denominator**2),
            3: 2*node*(node*node-3*a*a)/(denominator**3),
        }
        jet = [[fmpq(0) for _ in range(6)] for _ in alphas]
        jet[0][right+1] = fmpq(1, 2)
        for ai, alpha in enumerate(alphas[1:], 1):
            coefficient = derivatives[sum(alpha)] / prod(factorial(v) for v in alpha)
            for slope, power in zip(slopes[right], alpha):
                if power:
                    coefficient *= slope**power
            jet[ai][0] = coefficient
        row = [jet]
        for order in range(1, right+1):
            left = right-order
            d = node-k[left]
            dl = [u-v for u, v in zip(slopes[right], slopes[left])]
            P, Q = row[-1], table[right-1][order-1]
            quotient = []
            for ai in range(len(alphas)):
                coefficient = []
                for symbol in range(6):
                    numerator = P[ai][symbol]-Q[ai][symbol]
                    for j, previous in parents[ai]:
                        numerator -= dl[j]*quotient[previous][symbol]
                    coefficient.append(numerator/d)
                quotient.append(coefficient)
            row.append(quotient)
        table.append(row)
    return table[-1][-1], log_keys


def time_absolute_upper(a, breakpoints):
    """Upper dyadic for integral_a^infinity prod min(1/b_j,1/s) ds.

    Each rational b_j is a directed lower bound for sqrt(a^2+k_j^2),
    and b_j>=a. Thus this bounds the actual resolvent absolute integral.
    """
    b = sorted(breakpoints)
    assert a > 0 and b[0] >= a
    answer = arb((b[0]-a)/prod(b))
    if b[1] != b[0]:
        answer += arb(b[1]/b[0]).log()/arb(prod(b[1:]))
    for m in range(2, 5):
        answer += arb((b[m-1]**(1-m)-b[m]**(1-m))/
                      ((m-1)*prod(b[m:])))
    answer += arb(1/(4*b[4]**4))
    assert answer.is_finite() and answer > 0
    return rq(endpoint(answer, True))


def qualification(reg, vertices, source_hash, registration_hash):
    h = F(reg['four_tail_halfwidth'])
    rejected = []
    for ci, row in enumerate(reg['center_candidates_in_order']):
        point = tuple(map(F, row))
        D = lcm(*(x.denominator for x in point))
        center_n = [int(x*D) for x in point]
        nodes = [sum(v*w for v, w in zip(center_n, z)) for z in vertices]
        if len(set(nodes)) != 280:
            rejected.append({'candidate_index': ci, 'point': row,
                             'reason': 'center integer projections not all distinct',
                             'distinct_count': len(set(nodes))})
            continue
        assert point[0] == 1 and D == 104729
        assert all(0 < x-h < x+h < 1 for x in point[1:])
        crossings = []
        for vi, (z, node) in enumerate(zip(vertices, nodes)):
            norm = D*sum(abs(t) for t in z[1:])
            if norm and abs(F(node)) < h*norm:
                u = [F(-node*(1 if t > 0 else -1 if t < 0 else 0), norm)
                     for t in z[1:]]
                witness = [F(1)] + [x+t for x, t in zip(point[1:], u)]
                assert all(abs(t) < h for t in u)
                assert sum(x*t for x, t in zip(witness, z)) == 0
                crossings.append({'vertex_index': vi, 'z': z, 'center_k': node,
                                  'projection_interval_in_integer_scale':
                                  [str(F(node)-h*norm), str(F(node)+h*norm)],
                                  'zero_point': list(map(str, witness))})
        assert crossings, 'fixed box does not strictly cross any true-ray zero plane'
        # Prefer the requested x2=x3 plane; an antipodal ray gives repeated
        # zero nodes at the same witness, which the resolvent proof covers.
        preferred = [r for r in crossings if r['z'][0] == r['z'][3] == r['z'][4] == 0
                     and r['z'][1] == -r['z'][2] != 0]
        witness = preferred[0] if preferred else crossings[0]
        negative = [-v for v in witness['z']]
        opposite = vertices.index(negative)
        witness['opposite_vertex_index'] = opposite
        assert sum(F(x)*t for x, t in zip(witness['zero_point'], vertices[opposite])) == 0
        data = {'status': 'accepted_one_preregistered_center', 'candidate_index': ci,
                'rejected_candidates': rejected, 'point': list(map(str, point)),
                'D': D, 'all_280_center_projections_distinct': True,
                'minimum_abs_center_integer_projection': min(map(abs, nodes)),
                'closed_box_strictly_inside_unitface': True,
                'strict_zero_crossing_vertex_count': len(crossings),
                'zero_and_repeated_node_witness': witness,
                'source_sha256': source_hash, 'registration_sha256': registration_hash}
        atomic(HERE/'qualification.json', data)
        return point, D, nodes, data
    atomic(HERE/'qualification.json', {'status': 'all_three_registered_centers_rejected',
                                    'rejected_candidates': rejected})
    raise ValueError('all registered centers nongeneric; no new candidates allowed')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--work-dir', type=Path, required=True)
    parser.add_argument('--resume', action='store_true')
    args = parser.parse_args()
    began = time.monotonic()
    regpath = HERE/'registration.json'
    reg = json.loads(regpath.read_text())
    ctx.prec = reg['precision_bits']
    assert reg['Taylor_total_degree'] == 3 and ctx.prec == 384
    for row in reg['sources']:
        assert sha(row['source']) == row['sha256'], row['source']
    source_hash, reg_hash = sha(Path(__file__)), sha(regpath)
    fingerprint = hashlib.sha256(json.dumps({'source': source_hash, 'registration': reg_hash,
                                             'inputs': reg['sources']}, sort_keys=True).encode()).hexdigest()
    work = args.work_dir.resolve()
    assert not any((p/'.git').exists() for p in (work, *work.parents)), 'cache outside Git required'
    work.mkdir(parents=True, exist_ok=True)
    checkpoint = work/'checkpoint.json'
    if checkpoint.exists() and not args.resume:
        raise ValueError('existing checkpoint requires --resume')
    tablepath = Path(reg['sources'][0]['source'])
    table = json.loads(gzip.decompress(tablepath.read_bytes()))
    assert table['denominator'] == 144 and len(table['vertices']) == 280
    assert len(table['simplices']) == 3504
    vertices, simplices = table['vertices'], table['simplices']
    point, D, nodes, qual = qualification(reg, vertices, source_hash, reg_hash)
    h, eps = rq(reg['four_tail_halfwidth']), rq(reg['eps'])
    a = 144*D*eps
    b_lower = [rq(endpoint(arb(a*a+k*k).sqrt())) for k in nodes]
    assert min(b_lower) >= a > 0
    q_vertex = [h*D*sum(abs(v) for v in z[1:])/b
                for z, b in zip(vertices, b_lower)]
    assert max(q_vertex) < 1, 'fixed box resolvent condition failed; no width adjustment'
    alphas = sorted((x for x in product(range(4), repeat=4) if sum(x) <= 3),
                    key=lambda x: (sum(x), x))
    assert len(alphas) == 35
    index = {v: i for i, v in enumerate(alphas)}
    parents = [[(j, index[tuple(v-int(j == t) for t, v in enumerate(alpha))])
                for j in range(4) if alpha[j]] for alpha in alphas]
    stencil = [{} for _ in alphas]
    remainder_sum = arb(0)
    next_cone, segments = 0, []
    if args.resume:
        saved = json.loads(checkpoint.read_text())
        assert saved['fingerprint'] == fingerprint
        next_cone, segments = saved['next_simplex'], saved['segments']
        stencil = [{int(key): rq(value) for key, value in row.items()}
                   for row in saved['stencil']]
        remainder_sum = arb(rq(saved['remainder_sum_upper']))

    def save(done, status):
        atomic(checkpoint, {'fingerprint': fingerprint, 'status': status,
                           'next_simplex': done,
                           'stencil': [{str(k): str(v) for k, v in row.items()} for row in stencil],
                           'remainder_sum_upper': endpoint(remainder_sum, True),
                           'segments': segments+[{'pid': os.getpid(),
                                                  'seconds': time.monotonic()-began,
                                                  'end_simplex': done}]})

    time_cache = {}
    for ci in range(next_cone, 3504):
        simplex = simplices[ci]
        det, ids = int(simplex['det_numerator']), simplex['vertices']
        jet, log_keys = simplex_jet(ids, vertices, nodes, D, a, alphas, parents)
        for ai, row in enumerate(jet):
            for symbol, value in enumerate(row):
                if value == 0:
                    continue
                key = 0 if symbol == 0 else log_keys[symbol-1]
                combined = stencil[ai].get(key, fmpq(0))-det*value
                if combined:
                    stencil[ai][key] = combined
                else:
                    stencil[ai].pop(key, None)
        qc = max(q_vertex[i] for i in ids)
        tail = (1-qc)**(-5)-sum((fmpq(comb(n+4, 4))*qc**n for n in range(4)), fmpq(0))
        assert tail >= 0
        key = tuple(sorted(abs(nodes[i]) for i in ids))
        if key not in time_cache:
            time_cache[key] = time_absolute_upper(a, [b_lower[i] for i in ids])
        remainder_sum += arb(det*time_cache[key]*tail)
        done = ci+1
        if done % reg['checkpoint_stride_full_simplices'] == 0 or done == 3504:
            save(done, 'complete_stencils' if done == 3504 else 'checkpoint')
            print(f'completed {done}/3504 cutoff cones; seconds {time.monotonic()-began:.3f}', flush=True)
        if time.monotonic()-began >= reg['checkpoint_seconds'] and done < 3504:
            save(done, 'checkpoint')
            raise SystemExit(75)

    logs = {key: arb(key).log() for row in stencil for key in row if key}
    factor = arb(D)**4/(144*(2*arb.pi())**5)
    coefficients = [factor*sum((arb(c)*(logs[k] if k else 1) for k, c in row.items()), arb(0))
                    for row in stencil]
    variation = sum((abs(c)*arb(h**sum(alpha))
                     for alpha, c in zip(alphas[1:], coefficients[1:])), arb(0))
    full_remainder = factor*arb(rq(endpoint(remainder_sum, True)))
    lower_expression = coefficients[0]-variation-full_remainder
    upper_expression = coefficients[0]+variation+full_remainder
    L = max(fmpq(0), rq(endpoint(lower_expression)))
    U = rq(endpoint(upper_expression, True))
    assert U > L >= 0
    mean_polynomial = arb(0)
    for alpha, coefficient in zip(alphas, coefficients):
        if all(n % 2 == 0 for n in alpha):
            mean_polynomial += coefficient*arb(h**sum(alpha)/prod(n+1 for n in alpha))
    mean_lower = max(L, rq(endpoint(mean_polynomial-full_remainder)))
    mean_upper = min(U, rq(endpoint(mean_polynomial+full_remainder, True)))
    assert L <= mean_lower <= mean_upper <= U
    q = arb(5)/4
    chord_slope = (arb(U)**q-arb(L)**q)/arb(U-L)
    average_power_upper = arb(L)**q+chord_slope*arb(mean_upper-L)
    volume = (2*h)**4
    integral_upper = arb(volume)*average_power_upper
    constant_power_upper = arb(volume)*arb(U)**q
    assert integral_upper.is_finite() and integral_upper > 0
    final_bits = max(abs(int(s)).bit_length() for row in stencil for v in row.values()
                     for s in str(v).split('/'))
    rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    result = {
        'status': 'certified_one_closed_cutoff_potential_box',
        'object': reg['object'], 'model': reg['model'],
        'eps': str(eps), 'time_interval': '[1/512,infinity)',
        'box_center': list(map(str, point)), 'four_tail_halfwidth': str(h),
        'box_volume': str(volume), 'Taylor_total_degree': 3, 'precision_bits': ctx.prec,
        'all_simplices_processed': 3504, 'vertices': 280,
        'integer_scaling': {'D': D, 'a_144D_eps': str(a),
                            'factor': 'D^4/(144*(2*pi)^5)'},
        'maximum_resolvent_q_rational_upper': str(max(q_vertex)),
        'maximum_resolvent_q_display': str(arb(max(q_vertex))),
        'center_has_distinct_projections': True,
        'box_crosses_zero_projection': True,
        'zero_and_repeated_node_witness': qual['zero_and_repeated_node_witness'],
        'strict_zero_crossing_vertex_count': qual['strict_zero_crossing_vertex_count'],
        'Taylor_coefficients': [{'alpha': list(alpha), 'ball': str(c)}
                                for alpha, c in zip(alphas, coefficients)],
        'center_H_eps_ball': str(coefficients[0]),
        'Taylor_variation_upper_ball': str(variation),
        'cutoff_time_resolvent_remainder_upper_ball': str(full_remainder),
        'uniform_lower_expression_ball': str(lower_expression),
        'uniform_upper_expression_ball': str(upper_expression),
        'H_eps_uniform_lower_dyadic': str(L),
        'H_eps_uniform_upper_dyadic': str(U),
        'mean_Taylor_polynomial_ball': str(mean_polynomial),
        'mean_H_eps_lower_dyadic': str(mean_lower),
        'mean_H_eps_upper_dyadic': str(mean_upper),
        'chord_integral_upper_expression_ball': str(integral_upper),
        'integral_H_eps_5over4_upper_dyadic': endpoint(integral_upper, True),
        'constant_range_integral_upper_ball': str(constant_power_upper),
        'exact_stencil_entries': sum(len(row) for row in stencil),
        'maximum_stencil_integer_bits': final_bits,
        'checkpoint': {'path': str(checkpoint), 'sha256': sha(checkpoint),
                       'bytes': checkpoint.stat().st_size, 'fingerprint': fingerprint},
        'runtime': {'pid': os.getpid(), 'seconds': time.monotonic()-began,
                    'peak_rss_bytes': rss if sys.platform == 'darwin' else rss*1024,
                    'python': sys.executable, 'flint_version': flint.__version__,
                    'boxes': 1, 'mesh_cells': 0},
        'source_sha256': source_hash, 'registration_sha256': reg_hash,
        'external_sources': reg['sources'],
        'new_global_I_or_g_or_C5_bound': False,
        'small_time_0_to_eps_not_in_this_object': True,
    }
    atomic(HERE/'results.json', result)
    print(json.dumps({k: result[k] for k in ['status', 'center_H_eps_ball',
                      'maximum_resolvent_q_display', 'Taylor_variation_upper_ball',
                      'cutoff_time_resolvent_remainder_upper_ball',
                      'uniform_lower_expression_ball', 'uniform_upper_expression_ball',
                      'chord_integral_upper_expression_ball', 'runtime']}, indent=2))


if __name__ == '__main__':
    main()
