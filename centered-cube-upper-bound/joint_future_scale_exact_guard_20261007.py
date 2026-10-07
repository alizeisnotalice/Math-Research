#!/usr/bin/env python3
"""Exact finite selector/marginal and closed-face guards; not an actual FIRST sample.
No randomization; no old Beta D/M/K recomputation or cube phi numerics.
"""
from fractions import Fraction as F
from math import isqrt
from pathlib import Path
import hashlib
import json

HERE = Path(__file__).resolve().parent
CH = F(16, 3)


def frac_record(v):
    v = F(v)
    if v.numerator.bit_length()+v.denominator.bit_length() < 500:
        return f'{v.numerator}/{v.denominator}'
    num = v.numerator.to_bytes((v.numerator.bit_length()+7)//8, 'big')
    den = v.denominator.to_bytes((v.denominator.bit_length()+7)//8, 'big')
    raw = len(num).to_bytes(8, 'big')+num+len(den).to_bytes(8, 'big')+den
    return {'sha256_big_endian_length_framed': hashlib.sha256(raw).hexdigest(),
            'numerator_bits': v.numerator.bit_length(),
            'denominator_bits': v.denominator.bit_length()}


def block(N, width):
    return [F(1, width) if d < width else F(0) for d in range(N)]


def case(n, N, persistent):
    alpha = isqrt(n) + (isqrt(n)**2 != n)
    h = [block(N, k) for k in (1, 2, 4)]
    widths = (N, N, N) if persistent else (2, 4, N)
    phi = [[F(1, 8*N) + F(7, 8)*z for z in block(N, k)] for k in widths]
    avg0 = [[(h[j][d] + alpha*phi[j][d])/(alpha+1)
             for d in range(N)] for j in range(3)]
    # In this finite *linear* family avg(s) interpolates avg0 and phi.
    S = [max(avg0[j][d] for j in range(3)) for d in range(N)]
    S = [max(S[d], *(phi[j][d] for j in range(3))) for d in range(N)]
    C = sum(S)
    assert all(sum(row) == 1 for row in h + phi + avg0)
    assert all(sum(S[(x-y) % N] for x in range(N)) == C for y in range(N))
    receivers = (0, N//4, N//2, 3*N//4)
    q = F(1, N)  # Whole uniform source mu has mass one.
    rows = []
    for ix, x in enumerate(receivers):
        j = ix % 3
        s = F(2, n) if ix % 2 == 0 else F(1, 2)
        p = [(1-s)*h[j][(x-y) % N] + s*phi[j][(x-y) % N]
             for y in range(N)]
        av = [(1-s)*avg0[j][(x-y) % N] + s*phi[j][(x-y) % N]
              for y in range(N)]
        assert sum(p)/N == q and all(z > 0 for z in p)
        # Two history weights, a normalized hard source, and a correlated gate.
        r_weights = [F(1+(x+z) % 3) for z in range(N)]
        r = [z/sum(r_weights) for z in r_weights]
        w = [[F(1+(x+y) % 3, 8), F(1+(2*x+y) % 3, 8)] for y in range(N)]
        assert all(sum(z) <= 1 for z in w)
        accepted = []
        for y in range(N):
            accepted.append(sum(w[y][history] * sum(
                r[z] * F((x+3*y+5*z+history) % 7, 6) for z in range(N))
                for history in range(2)))
        assert all(F(0) <= a <= 1 for a in accepted)
        rows.append((x, p, av, accepted))
    reports = []
    # Include an exact ratio threshold: equality must be assigned to paid.
    equality_tau = S[0]/rows[0][1][0]
    for tau in (F(1, 4), F(1, 8), equality_tau):
        actual = F(0)
        marginal = F(0)
        kernel_bound = F(0)
        old_paid = new_paid = rescued = hot = equality_count = 0
        for x, p, av, accepted in rows:
            paid = [S[(x-y) % N] >= tau*p[y] for y in range(N)]
            old = [av[y] >= tau*p[y] for y in range(N)]
            assert all(not old[y] or paid[y] for y in range(N))
            actual += CH*q*sum(p[y]*accepted[y] for y in range(N) if paid[y])
            marginal += CH*sum(p[y]/N for y in range(N) if paid[y])
            kernel_bound += CH/tau*sum(S[(x-y) % N]/N for y in range(N))
            for y in range(N):
                old_paid += old[y]
                new_paid += paid[y]
                rescued += paid[y] and not old[y]
                hot += not paid[y]
                equality_count += S[(x-y) % N] == tau*p[y]
            hot_response = sum(p[y]/N for y in range(N) if not paid[y])
            for jj in range(3):
                # Endpoints certify every s in this affine surrogate family.
                for kernel in (avg0[jj], phi[jj]):
                    test_response = sum(kernel[(x-y) % N]/N
                                        for y in range(N) if not paid[y])
                    assert test_response <= tau*hot_response
        global_bound = CH*C/tau  # W=1; receivers form a subset of full group.
        assert actual <= marginal <= kernel_bound <= global_bound
        reports.append({'tau': frac_record(tau), 'actual': frac_record(actual),
                        'original_y_marginal_bound': frac_record(marginal),
                        'selected_receiver_S_bound': frac_record(kernel_bound),
                        'one_W_bound': frac_record(global_bound),
                        'old_paid_count': old_paid, 'new_paid_count': new_paid,
                        'rescued_old_low_count': rescued, 'retained_hot_count': hot,
                        'paid_equality_count': equality_count, 'exact_PASS': True})
    return {'n': n, 'alpha': alpha, 'N': N,
            'family': 'hard-persistent' if persistent else 'sibling-rescue',
            'whole_source_mass_W': '1/1', 'q': frac_record(q),
            'fixed_S_column_mass_C': frac_record(C), 'reports': reports,
            'scope': 'Finite cyclic linear kernel family with correlated positive gates; not original cubes/FIRST/future-cap samples.'}


def closed_faces(n):
    a, b = F(1), F(2)
    records = []
    for threshold in (a, F(3, 2), b):
        exact_sup = threshold**(-n)
        if threshold < b:
            previous = F(0)
            approximants = []
            for k in (1, 2, 4, 8):
                L = threshold + (b-threshold)/2**k
                value = L**(-n)
                assert previous <= value <= exact_sup
                assert value/exact_sup == (threshold/L)**n
                approximants.append({'right_L': frac_record(L),
                                     'relative_value': frac_record(value/exact_sup)})
                previous = value
            records.append({'face_threshold': frac_record(threshold),
                            'right_approximation_PASS': True,
                            'approximants': approximants})
        else:
            # A dense interior alone misses the closed face exactly at b.
            interior_value = F(0)  # Every L<b misses |v|_infty=b/2.
            assert interior_value == 0 and exact_sup > 0
            records.append({'face_threshold': '2/1', 'b_endpoint_required': True,
                            'interior_sup': '0/1',
                            'endpoint_closed_value': frac_record(exact_sup)})
    return {'n': n, 'exact_PASS': True, 'records': records,
            'scope': 'Exact hard centered-cube kernel closed-face/dense-right checks; no soft phi or cap certification.'}


def main():
    rounds = []
    for n, N in ((512, 8), (1024, 16), (4096, 32)):
        rounds.append({'n': n, 'selector_cases': [case(n, N, False), case(n, N, True)],
                       'closed_faces': closed_faces(n), 'exact_PASS': True})
    data = {'seed': None, 'randomization': 'none',
            'scope': 'Only new sup-selector algebra, marginal domination/source-once chain, inclusion/equality and closed-face measurability guards. No D/M/K recheck, actual FIRST simulation, or space-budget closure.',
            'all_exact_PASS': True, 'rounds': rounds}
    dest = HERE/'joint_future_scale_exact_guard_20261007_results.json'
    dest.write_text(json.dumps(data, indent=2, ensure_ascii=False)+'\n')
    print(json.dumps({'all_exact_PASS': True, 'rounds': 3, 'selector_cases': 6,
                      'threshold_cases': 18, 'closed_face_cases': 9,
                      'rescued': sum(r['rescued_old_low_count'] for a in rounds
                                     for c in a['selector_cases'] for r in c['reports']),
                      'retained': sum(r['retained_hot_count'] for a in rounds
                                      for c in a['selector_cases'] for r in c['reports']),
                      'results': str(dest)}, ensure_ascii=False))


if __name__ == '__main__':
    main()
