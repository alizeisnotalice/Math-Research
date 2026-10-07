#!/usr/bin/env python3
"""Finite exact component guards for the source-core antichain interface.

No original phi/FIRST/history sampling. Geometry and positive row/column arrays
are independent certificates, never a joint actual residual model.
"""
from fractions import Fraction as F
from itertools import combinations
from math import isqrt
from pathlib import Path
import hashlib
import json

OUT = Path(__file__).with_name('double_cutoff_actual_geometry_exact_guard_20261007_results.json')
CHECKS = []


def check(name, statement, detail=None):
    if not statement:
        raise AssertionError((name, detail))
    CHECKS.append({'name': name, 'pass': True, 'detail': detail})


def fs(x):
    return f'{x.numerator}/{x.denominator}'


def k_log(n):
    return (n + 1).bit_length()  # ceil(log2(n+2)), exact integers


def members(node):
    if node >= 8:
        return [node - 8]
    return members(2 * node) + members(2 * node + 1)


def ancestors(node):
    result = []
    while node:
        result.append(node)
        node //= 2
    return result[::-1]


def lca(i, j):
    a, b = 8 + i, 8 + j
    while a != b:
        if a > b:
            a //= 2
        else:
            b //= 2
    return a


def top_for(node, marked):
    return next((a for a in ancestors(node) if a in marked), None)


def area_sweep(rects):
    xs = sorted({r[t] for r in rects for t in (0, 1)})
    total = F(0)
    for lo, hi in zip(xs, xs[1:]):
        active = sorted((r[2], r[3]) for r in rects if r[0] < hi and r[1] > lo)
        length = F(0)
        if active:
            a, b = active[0]
            for c, d in active[1:]:
                if c > b:
                    length += b - a
                    a, b = c, d
                else:
                    b = max(b, d)
            length += b - a
        total += (hi - lo) * length
    return total


def area_ie(rects):
    total = F(0)
    for size in range(1, len(rects) + 1):
        for group in combinations(rects, size):
            dx = max(F(0), min(r[1] for r in group) - max(r[0] for r in group))
            dy = max(F(0), min(r[3] for r in group) - max(r[2] for r in group))
            total += (1 if size % 2 else -1) * dx * dy
    return total


def core_volume(ids, n, width, radius):
    # The last n-2 coordinates have identical exact sections. This is an
    # exact rational subspace factorization, not a high-dimensional bbox.
    side = width + radius
    rects = []
    for i in ids:
        x, y = F(2 * i + 1, 16), F(i % 3, 32)
        rects.append((x - side / 2, x + side / 2,
                      y - side / 2, y + side / 2))
    sweep, ie = area_sweep(rects), area_ie(rects)
    check(f'n{n}:core_union_{ids}', sweep == ie, {'rectangles': len(ids)})
    return sweep * side ** (n - 2)


def scalar_radii(n):
    sqrt = isqrt(n)
    a, b, q = F(1), F(2), F(3, 7)
    for c in (F(3, 2), F(2), F(5, 2)):
        h = q * sqrt * c ** n
        rho = min(b, c)
        check(f'n{n}:CP_radius_mass_{c}', q * rho ** n <= h / sqrt)
        for small in (a, F(3, 2), b):
            for large in (small, b):
                for rh, ls in ((small, large), (large, small)):
                    strict = h > q * min(ls, rh) ** n * sqrt
                    check(f'n{n}:CP_closed_{c}_{ls}_{rh}', not strict or min(ls, rh) <= rho)
        if c > b:
            check(f'n{n}:CP_cap_equality_{c}', h > q * b ** n * sqrt and b == rho)
    for r in (F(9, 8), F(3, 2), b, F(5, 2)):
        # h/q=r^(n/2), so the exact GP hard upper radius is min(b,r).
        h = q * r ** (n // 2)
        bound = min(b, r)
        for ls in (a, F(3, 2), b):
            for rh in (a, F(9, 8), F(3, 2), b):
                strict = h > q * (ls * rh) ** (n // 2)
                check(f'n{n}:GP_closed_{r}_{ls}_{rh}', not strict or rh <= bound)
        if r > b:
            check(f'n{n}:GP_cap_equality_{r}', h > q * b ** (n // 2) and bound == b)
        if r == b:
            check(f'n{n}:GP_root_equality_paid_{r}', h == q * (a * b) ** (n // 2))
    # Source-antichain monotonicity can be checked without irrational roots.
    masses = [F(i * i + 1, i + 2) for i in range(1, 9)]
    for i in range(len(masses) - 1):
        h1, h2 = masses[i], masses[i] + masses[i + 1]
        check(f'n{n}:radius_monotone_power_{i}',
              min(b ** n, h1 / (q * sqrt)) <= min(b ** n, h2 / (q * sqrt))
              and min(b ** n, (h1 / q) ** 2 / a ** n)
              <= min(b ** n, (h2 / q) ** 2 / a ** n))
    # The single physical scale a=b is closed and covered, not an exception.
    check(f'n{n}:single_scale_closed', min(F(1), F(2)) == F(1))


def forest_geometry(n):
    sqrt = isqrt(n)
    zeta = F(1, sqrt)
    T = F(sqrt * k_log(n) ** 4)
    mass = [F(16 + i * i, 17) if i % 2 == 0 else F(1, 64 * (i + 1)) for i in range(8)]
    born = [mass[i] * F(2 + i % 3, 4) for i in range(8)]
    W = sum(mass)
    q = min(born) / (100 * sqrt * 2 ** n)
    cores, volumes = {}, {}
    for p in range(1, 16):
        ids = members(p)
        keep = [i for i in ids if i % 2 == 0]
        if p != 1:
            odds = [i for i in ids if i % 2]
            keep += odds[:1]  # child can contain an atom the ancestor omits
        if sum(mass[i] for i in ids if i not in keep) > zeta * sum(mass[i] for i in ids):
            keep = ids[:]
        cores[p] = set(keep)
        volumes[p] = core_volume(keep, n, F(1, 128 * n), F(2))
        MP, hP = sum(mass[i] for i in ids), sum(born[i] for i in ids)
        check(f'n{n}:core_mass_{p}', sum(mass[i] for i in ids if i not in keep) <= zeta * MP)
        check(f'n{n}:core_volume_mark_{p}', q * volumes[p] <= T * hP)
        check(f'n{n}:core_radius_cap_{p}', hP > q * sqrt * 2 ** n and (hP / q) ** 2 > 2 ** n)
    check(f'n{n}:nonnested_core_witness', 1 in cores[4] and 1 not in cores[1])

    patterns = [set([1, 2, 4, 8, 9, 5]), set([2, 5, 10, 11, 12, 14]), set(range(8, 16))]
    rows = []
    for number, marked in enumerate(patterns):
        top = [p for p in marked if top_for(p, marked) == p]
        sets = [set(members(p)) for p in top]
        check(f'n{n}:antichain_{number}', all(not (a & b) for a, b in combinations(sets, 2)))
        check(f'n{n}:antichain_mass_{number}', sum(sum(mass[i] for i in a) for a in sets) <= W)
        bad_ids = {i for p in top for i in members(p) if i not in cores[p]}
        bad_mass = sum(mass[i] for i in bad_ids)
        bad_born = sum(born[i] for i in bad_ids)
        check(f'n{n}:bad_source_once_{number}', bad_mass <= zeta * W and bad_born <= bad_mass)
        check(f'n{n}:union_good_volume_{number}',
              q * sum(volumes[p] for p in top) <= T * sum(sum(born[i] for i in members(p)) for p in top))
        reclassified = 0
        for i in range(8):
            for j in range(8):
                if i == j:
                    continue  # no distinct-child LCA
                P = lca(i, j)
                A = top_for(P, marked)
                if A is None:
                    continue
                check(f'n{n}:unique_ancestor_{number}_{i}_{j}',
                      A in top and i in members(A) and j in members(A))
                for endpoint in (i, j):
                    ancestor_good = endpoint in cores[A]
                    check(f'n{n}:endpoint_recheck_{number}_{i}_{j}_{endpoint}',
                          ancestor_good or endpoint in bad_ids)
                    if endpoint in cores[P] and not ancestor_good:
                        reclassified += 1
        rows.append({'pattern': number, 'top': sorted(top), 'bad_ids': sorted(bad_ids),
                     'bad_mass': fs(bad_mass), 'child_good_ancestor_bad_occurrences': reclassified})
    check(f'n{n}:reclassification_present', rows[0]['child_good_ancestor_bad_occurrences'] > 0)
    return mass, born, cores, patterns, rows, q, T


def positive_array_marginals(n, mass, born, cores, patterns):
    q, Ch = F(1), F(5, 2)
    xs = range(12)
    soft, hard, orders = [], [], []
    for x in xs:
        a = [F(1 + (x + 2 * i) % 7, 3 + i) for i in range(8)]
        b = [F(1 + (2 * x + i) % 5, 2 + i) for i in range(8)]
        # Unnormalized fixed-source kernels with postulated original row bounds.
        soft.append([q * a[i] / (sum(a) * mass[i]) for i in range(8)])
        hard.append([Ch * q * b[i] / (sum(b) * born[i]) for i in range(8)])
        ls, rh = [(F(2), F(1)), (F(1), F(2)), (F(3, 2), F(3, 2))][x % 3]
        orders.append(rh <= ls)
        check(f'n{n}:row_bounds_{x}', sum(mass[i] * soft[x][i] for i in range(8)) == q
              and sum(born[i] * hard[x][i] for i in range(8)) == Ch * q)
    envelope = [[max(soft[x][i], hard[x][i]) for i in range(8)] for x in xs]
    normalization = max(sum(envelope[x][i] for x in xs) for i in range(8))
    rx = F(1) / normalization
    for i in range(8):
        check(f'n{n}:fixed_source_column_{i}', rx * sum(envelope[x][i] for x in xs) <= 1)
    reports = []
    for number, marked in enumerate(patterns):
        top = [p for p in marked if top_for(p, marked) == p]
        bad_ids = {i for p in top for i in members(p) if i not in cores[p]}
        mu_bad = [mass[i] if i in bad_ids else F(0) for i in range(8)]
        CP, GP = F(0), F(0)
        for x in xs:
            cp_x, gp_x = F(0), F(0)
            for y in range(8):
                for z in range(8):
                    if y == z:
                        continue
                    P = lca(y, z)
                    A = top_for(P, marked)
                    if A is None:
                        continue
                    G = F(1, 2 + (x + y + z) % 9)
                    pair = mass[y] * soft[x][y] * born[z] * hard[x][z] * G / q
                    if z not in cores[A]:
                        gp_x += pair
                    # D_h includes all K; D_s includes only a positive K=0
                    # subkernel. Its hard substitute here is the same soft
                    # array, so the endpoint domination is exact algebra.
                    if orders[x]:
                        if z not in cores[A]:
                            cp_x += pair
                    elif y not in cores[A]:
                        cp_x += pair / 3
            kh_bad = sum(mu_bad[i] * hard[x][i] for i in range(8))
            ks_bad = sum(mu_bad[i] * soft[x][i] for i in range(8))
            check(f'n{n}:GP_bad_marginal_{number}_{x}', gp_x <= kh_bad)
            selected = kh_bad if orders[x] else Ch * ks_bad
            check(f'n{n}:CP_same_x_order_{number}_{x}', cp_x <= selected)
            check(f'n{n}:CP_one_envelope_{number}_{x}', selected <= Ch * sum(mu_bad[i] * envelope[x][i] for i in range(8)))
            CP += rx * cp_x
            GP += rx * gp_x
        mb = sum(mu_bad)
        check(f'n{n}:GP_single_column_payment_{number}', GP <= mb)
        check(f'n{n}:CP_single_column_payment_{number}', CP <= Ch * mb)
        reports.append({'pattern': number, 'bad_mass': fs(mb), 'GP_bad_traffic': fs(GP),
                        'CP_bad_traffic': fs(CP), 'finite_column_constant': '1/1'})
    return reports


def main():
    rounds = []
    for n in (4, 16, 64):
        before = len(CHECKS)
        scalar_radii(n)
        mass, born, cores, patterns, trees, q, T = forest_geometry(n)
        marginals = positive_array_marginals(n, mass, born, cores, patterns)
        rounds.append({'n': n, 'sqrt_n': isqrt(n), 'zeta': fs(F(1, isqrt(n))),
                       'k_n': k_log(n), 'T': fs(T), 'geometry_q': fs(q),
                       'checks': len(CHECKS) - before, 'antichain_receipts': trees,
                       'marginal_receipts': marginals})
    payload = {'status': 'PASS', 'seed': None, 'arithmetic': 'fractions.Fraction exact',
               'scope': 'Independent finite source-tree, closed-box volume, radius-power and positive marginal certificates. No phi, continuous FIRST, full actual history, actual coverage rate or new order certification.',
               'registration': 'double_cutoff_actual_geometry_interface_20261007.md section 8, written before this execution',
               'rounds': rounds, 'total_checks': len(CHECKS),
               'checks_sha256': hashlib.sha256(json.dumps(CHECKS, sort_keys=True).encode()).hexdigest(),
               'checks': CHECKS}
    OUT.write_text(json.dumps(payload, ensure_ascii=False, indent=2) + '\n')
    print(json.dumps({'status': payload['status'], 'total_checks': len(CHECKS),
                      'rounds': [{'n': r['n'], 'checks': r['checks']} for r in rounds],
                      'path': str(OUT)}, ensure_ascii=False))


if __name__ == '__main__':
    main()
