#!/usr/bin/env python3
"""Finite fine-ball cap sums. Decimal crosschecks are not interval certificates.

For odd n, F_n is the exact polynomial CDF of one coordinate in B(0,1).
Let s=(rho^2+t^2-1)/(2rho). The outside probability is
F_n(s/t)-t^(-n) F_n(s-rho), by the two intersection caps.
All transcendental evaluations use Python's standard-library Decimal.
"""
from decimal import Decimal as D, localcontext, ROUND_CEILING, ROUND_FLOOR
from fractions import Fraction as F
from math import comb, factorial
from pathlib import Path
import argparse, hashlib, json, time

BATCHES = ((1, 3, 9), (17, 33, 65), (97, 129, 161))


def ceil_log2(k):
    return (k - 1).bit_length()


def coordinate_cdf(n):
    assert n >= 1 and n % 2
    m = (n - 1) // 2
    normal = F(factorial(2*m+1), 2**(2*m+1)*factorial(m)**2)
    rational = [normal*F((-1)**k*comb(m, k)*2**(m-k), m+k+1)
                for k in range(m+1)]
    coefficients = [D(c.numerator)/D(c.denominator) for c in rational]

    def cdf(x):
        if x <= -1:
            return D(0)
        if x >= 1:
            return D(1)
        if x > 0:
            return 1-cdf(-x)
        u = 1+x
        p = D(0)
        for c in reversed(coefficients):
            p = p*u+c
        return u**(m+1)*p
    return cdf


def geometry_case(n, kind, parameter, precision):
    with localcontext() as ctx:
        ctx.prec = precision
        log2 = D(2).ln()
        if kind == 'closed_form':
            b = parameter
            rho = 1-D(2)**(-b)
        else:
            rho = (-D(parameter)*log2/n).exp()
        delta = 1-rho*rho
        gap = 1-rho
        K = (1+D(n)/2*(D(2)/delta).ln()/log2
             +D(n)/2*(1+16*rho*rho/(n*delta)).ln()/log2)
        triangle_depth = D(n)*(1/gap).ln()/log2
        N = ((parameter-1) if kind == 'closed_form' else
             int(triangle_depth.to_integral_value(rounding=ROUND_CEILING))-1)
        cdf = coordinate_cdf(n)
        values, masses = [], []
        # Direct powers avoid accumulated drift in recursively generated radii.
        for d in range(1, N+1):
            t = (-D(d)*log2/n).exp() if n != 1 else D(2)**(-d)
            assert t > gap
            s = (rho*rho+t*t-1)/(2*rho)
            a = cdf(s/t)-cdf(s-rho)/(t**n)
            assert 0 <= a <= 1, (n, kind, parameter, d, a)
            values.append(a)
            masses.append(D(2)**(-d)*a)
        following = (-D(N+1)*log2/n).exp() if n != 1 else D(2)**(-N-1)
        # The n=1 equality cases are evaluated exactly as powers of two.
        assert following <= gap
        A = sum(values, D(0))
        a1 = values[0] if values else D(0)
        G = (A+a1)/2
        # Disjoint annuli C_d\C_(d+1) attain the integrated supremum.
        annuli = [q-(masses[i+1] if i+1 < len(masses) else D(0))
                  for i, q in enumerate(masses)]
        assert all(q >= 0 for q in annuli)
        annulus_integral = sum((D(2)**(d+1)*q for d, q in enumerate(annuli)), D(0))
        tol = D(10)**(-(precision//2))
        assert abs(annulus_integral-G) <= tol
        assert A <= K and A <= N and G <= A
        assert 2*rho**n*G <= rho**n*(K+1)
        row = dict(n=n, kind=kind, parameter=parameter, rho=str(rho),
                   delta=str(delta), terms=N, A=str(A), a1=str(a1), G=str(G),
                   K=str(K), triangle_N=N, K_minus_A=str(K-A),
                   triangle_minus_A=str(D(N)-A), A_minus_G=str(A-G),
                   annulus_G_difference=str(abs(annulus_integral-G)),
                   annuli_disjoint_by_nested_ball_definition=True)
        if kind == 'closed_form':
            exact_A = D(parameter-2)/2+D(2)**(-parameter)
            exact_G = D(parameter-1)/4
            assert abs(A-exact_A) <= tol and abs(G-exact_G) <= tol
            row.update(exact_A=str(exact_A), exact_G=str(exact_G))
        elif kind == 'lower_bound':
            floor_value = int((D(n)/2*(1/delta).ln()/log2-1)
                              .to_integral_value(rounding=ROUND_FLOOR))
            lower = D(max(0, floor_value))/4
            assert A >= lower
            row.update(lower_bound=str(lower), A_minus_lower_bound=str(A-lower))
        elif kind == 'cutoff':
            B = n*ceil_log2(2*n)
            # Exact integer characterization of ceil(log2(B/2+1)).
            L = 0
            while 2**(L+1) < B+2:
                L += 1
            assert L == parameter
            envelope = rho**n*(K+1)
            assert K <= D(B)/2 and envelope <= 1
            row.update(B=B, L=L, B_half_minus_K=str(D(B)/2-K),
                       rho_n_K_plus_one=str(envelope), one_minus_envelope=str(1-envelope))
        return row


def rational_kernel_batches():
    rows = []
    # Annulus masses are exact rational numbers; each q_d is nested ball mass.
    # The density at shell d is max(2^j:j<=d)=2^d.
    for batch, lengths in enumerate(((1, 2, 4), (8, 16, 32), (64, 96, 128)), 1):
        count = 0
        for length in lengths:
            for profile in range(3):
                shell = [F((d+profile)%7+1, (d+2)**2*(profile+2))
                         for d in range(1, length+1)]
                q = [sum(shell[d:], F(0)) for d in range(length)]
                a = [2**(d+1)*mass for d, mass in enumerate(q)]
                left = sum((2**(d+1)*mass for d, mass in enumerate(shell)), F(0))
                right = (sum(a, F(0))+a[0])/2
                assert left == right
                assert all(q[d]-q[d+1] == shell[d] for d in range(length-1))
                assert q[-1] == shell[-1]
                # Volume coordinate v=|w|^n/R_j^n makes all endpoints rational.
                # The final truncated shell includes the remaining inner ball.
                shells = [(F(1, 2**(d+1)) if d < length else F(0),
                           F(1, 2**d)) for d in range(1, length+1)]
                assert shells[-1][0] == 0
                for i, (lo, hi) in enumerate(shells):
                    assert 0 <= lo < hi
                    if i+1 < len(shells):
                        assert shells[i+1][1] == lo
                    v = (lo+hi)/2
                    assert sum(low <= v < high for low, high in shells) == 1
                    active = [j for j in range(1, length+1) if v < F(1, 2**j)]
                    assert active == list(range(1, i+2))
                for boundary in [F(0)]+[lo for lo, hi in shells]:
                    assert sum(lo <= boundary < hi for lo, hi in shells) == 1
                for d in range(1, length+1):
                    assert max(2**j for j in range(1, d+1)) == 2**d
                    Vj = F(length+profile+1, profile+2)
                    kernels = [F(2**j)/Vj for j in range(1, d+1)]
                    assert max(kernels) == (sum(kernels, F(0))+kernels[0])/2
                    count += 1
        rows.append(dict(batch=batch, lengths=list(lengths), profiles=3,
                         exact_shell_checks=count, identity_checks=3*len(lengths),
                         arithmetic='Fraction', all_passed=True))
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    start = time.monotonic()
    rows = []
    for batch, dimensions in enumerate(BATCHES, 1):
        batch_rows = []
        for n in dimensions:
            if n == 1:
                cases = [('closed_form', b) for b in range(1, 5)]
            else:
                cases = [('lower_bound', a) for a in sorted({1, ceil_log2(n)})]
                if n >= 32:
                    B = n*ceil_log2(2*n)
                    L = 0
                    while 2**(L+1) < B+2:
                        L += 1
                    cases.append(('cutoff', L))
            for kind, parameter in cases:
                low = geometry_case(n, kind, parameter, 120)
                high = geometry_case(n, kind, parameter, 160)
                with localcontext() as ctx:
                    ctx.prec = 180
                    difference = max(abs(D(low[key])-D(high[key]))
                                     for key in ('A', 'a1', 'G', 'K'))
                assert difference < D('1e-70'), (n, kind, parameter, difference)
                high.update(batch=batch, precisions=[120, 160],
                            cross_precision_max_absolute_difference=str(difference),
                            interval_certificate=False)
                rows.append(high)
                batch_rows.append(high)
        print('batch', batch, 'dimensions', ','.join(map(str, dimensions)),
              'cases', len(batch_rows), 'max_terms', max(r['terms'] for r in batch_rows),
              'all_passed', flush=True)
    out = dict(status='passed', geometry=rows, dimensions=[list(b) for b in BATCHES],
               geometry_case_count=len(rows), rational_kernel=rational_kernel_batches(),
               arithmetic='standard-library Decimal at 120/160 digits; exact Fraction checks',
               cross_precision_difference_is_certified_error_bound=False,
               actual_E_budget_counterexample=False, mathematical_main_goal_solved=False,
               scope='fine-ball geometric cap bound and integrated disjoint-mask supremum only',
               elapsed_seconds=time.monotonic()-start,
               script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    args.output.write_text(json.dumps(out, indent=2)+'\n')
    print('PASSED cases', len(rows), 'seconds', round(out['elapsed_seconds'], 2), flush=True)


if __name__ == '__main__':
    main()
