#!/usr/bin/env python3
"""New exact source-contract guards, not original phi/FIRST/history data.

Finite geometry, positive arrays, and compressed large-n ratios are separate
components. No previously saved experiment is imported or rerun.
"""
from fractions import Fraction as F
from itertools import combinations
from pathlib import Path
import hashlib
import json

OUT = Path(__file__).with_name('adaptive_core_source_budget_exact_guard_20261007_results.json')
checks = []


def ck(label, value, detail=None):
    if not value:
        raise AssertionError((label, detail))
    checks.append({'label': label, 'pass': True, 'detail': detail})


def fs(x):
    return f'{x.numerator}/{x.denominator}'


NODES = {1: list(range(11)), 2: [0], 3: list(range(1, 11)),
         6: list(range(1, 6)), 7: list(range(6, 11)),
         12: [1, 2], 13: [3, 4, 5], 14: [6, 7], 15: [8, 9, 10],
         30: [8, 9], 31: [10]}
POINTS = [F(1), F(520), F(530), F(650), F(660), F(670),
          F(780), F(790), F(910), F(920), F(970)]


def path(i):
    return sorted((p for p, ids in NODES.items() if i in ids), key=lambda p: p.bit_length())


def min_core(ids, loss_limit, mass):
    total = sum(mass[i] for i in ids)
    for size in range(len(ids) + 1):
        for keep in combinations(ids, size):
            if total - sum(mass[i] for i in keep) <= loss_limit:
                return tuple(keep)
    raise AssertionError('full core must exist')


def geometry(n):
    sqrt = int(n ** 0.5)
    zeta = F(1, sqrt)
    k = (n + 1).bit_length()
    T = F(sqrt * k ** 4)
    mass = [F(1)] + [F(19 + (i - 1) % 3, 100) for i in range(1, 10)] + [F(1, 1000)]
    born = [mass[i] * F(2 + i % 3, 4) for i in range(11)]
    W = sum(mass)
    MP = {p: sum(mass[i] for i in ids) for p, ids in NODES.items()}
    D = max(len(path(i)) for i in range(11))
    ck(f'n{n}:root_leaf_count', D == 5)
    ck(f'n{n}:point_padding_disjoint',
       all(abs(a - b) > 2 for a, b in combinations(POINTS, 2)))
    # Cubes around these point cores are pairwise disjoint for r<=b=2.
    # Exact n-dimensional union volume is count*r^n, never a total bbox.
    q = F(99, 100) * T / (2 ** n)
    old = {}
    for p, ids in NODES.items():
        core = min_core(ids, zeta * MP[p], mass)
        certified = q * len(core) * 2 ** n <= T * MP[p]
        if certified:
            old[p] = core
    ck(f'n{n}:old_hard_child_only_marked', set(old) == {2}, {'marked': sorted(old)})
    old_top = [p for p in old if not any(a in old for a in path(NODES[p][0]) if a != p and a.bit_length() < p.bit_length())]
    ck(f'n{n}:old_coarsest_antichain', old_top == [2])
    S_b = set(NODES[2])
    weights = {p: F(1) if p == 2 else F(1, D) if not (set(ids) & S_b) else F(0)
               for p, ids in NODES.items()}
    for i in range(11):
        ck(f'n{n}:path_budget_{i}', sum(weights[p] for p in path(i)) <= 1)
    charge = sum(weights[p] * MP[p] for p in NODES)
    ck(f'n{n}:whole_mass_charge', charge <= W)
    ck(f'n{n}:unused_mass_charge',
       sum(weights[p] * MP[p] for p in NODES if p != 2) <= W - mass[0])
    radii = [F(1)]
    while radii[-1] < 2:
        radii.append(min(F(2), radii[-1] * F(n + 1, n)))
    ck(f'n{n}:closed_grid', radii[0] == 1 and radii[-1] == 2 and len(radii) == len(set(radii)))
    certs = {2: (F(2), old[2])}
    for p, ids in NODES.items():
        if p == 2 or weights[p] == 0:
            continue
        w = weights[p]
        core = min_core(ids, zeta * w * MP[p], mass)
        ok = [r for r in radii if q * len(core) * r ** n <= T * w * MP[p]]
        if ok:
            certs[p] = (max(ok), core)
            ck(f'n{n}:only_one_maximal_radius_{p}',
               all(q * len(core) * r ** n > T * w * MP[p] for r in radii if r > max(ok)))
        else:
            ck(f'n{n}:no_default_radius_{p}', q * len(core) > T * w * MP[p])
    bad = set()
    summed_volume = F(0)
    individual_loss = F(0)
    max_core_radius = {}
    for p, (r, core) in certs.items():
        loss = MP[p] - sum(mass[i] for i in core)
        ck(f'n{n}:saved_core_quota_{p}', loss <= zeta * weights[p] * MP[p])
        ck(f'n{n}:saved_core_volume_{p}', q * len(core) * r ** n <= T * weights[p] * MP[p])
        summed_volume += len(core) * r ** n
        individual_loss += loss
        bad.update(set(NODES[p]) - set(core))
        for i in core:
            max_core_radius[i] = max(r, max_core_radius.get(i, F(0)))
    # The exact U is a union of disjoint point-centred cubes; repeated
    # certificates at one point form nested cubes, even though cores do not.
    Uvolume = sum(r ** n for r in max_core_radius.values())
    bad_mass = sum(mass[i] for i in bad)
    ck(f'n{n}:union_not_bbox', Uvolume <= summed_volume)
    ck(f'n{n}:global_good_budget', q * Uvolume <= T * charge <= T * W)
    ck(f'n{n}:bad_set_union_budget', bad_mass <= individual_loss <= zeta * charge <= zeta * W)
    ck(f'n{n}:born_bad_subset', sum(born[i] for i in bad) <= bad_mass)
    reclassified = 0
    for i in range(11):
        for R in radii:
            eligible = [p for p in path(i) if p in certs and R <= certs[p][0]]
            if not eligible:
                continue
            P = eligible[0]
            ancestor_good = i in certs[P][1]
            ck(f'n{n}:endpoint_recheck_{i}_{R}', ancestor_good or i in bad)
            # Receiver lies on the ORIGINAL hard cube's positive closed face.
            x = POINTS[i] + R / 2
            if ancestor_good:
                ck(f'n{n}:closed_face_support_{i}_{R}',
                   abs(x - POINTS[i]) <= max_core_radius[i] / 2)
            if not ancestor_good and any(i in certs[a][1] for a in eligible[1:]):
                reclassified += 1
    if n >= 16:
        ck(f'n{n}:nonnested_certificate_reclassification', reclassified > 0)
    # A strictly new source-contract eligibility: y outside old marked A,
    # unmarked original LCA, hard z inside marked hard-child A.
    y, z, A = 1, 0, 2
    LCA = max(set(path(y)) & set(path(z)), key=lambda p: p.bit_length())
    ck(f'n{n}:hard_child_only_receipt', y not in NODES[A] and z in NODES[A]
       and LCA == 1 and LCA not in old and A in certs and certs[A][0] == 2)
    return mass, born, certs, bad, radii, {
        'n': n, 'D': D, 'zeta': fs(zeta), 'T': fs(T), 'q': fs(q),
        'source_mass': fs(W), 'charged_mass': fs(charge),
        'old_marked': sorted(old), 'new_certificate_nodes': sorted(p for p in certs if p != 2),
        'radii': {str(p): fs(r) for p, (r, _) in certs.items()},
        'bad_union_ids': sorted(bad), 'bad_union_mass': fs(bad_mass),
        'non_nested_reclassifications': reclassified,
        'hard_child_only': {'y': y, 'z': z, 'node': A, 'original_LCA': LCA,
                            'old_LCA_ancestor_qualification': False,
                            'new_hard_membership_qualification': True,
                            'scope': 'source selector certificate; all other actual gates unknown'}}


def positive_arrays(n, mass, born, certs, bad, radii):
    q, Ch = F(1), F(16, 3)
    R = [F(1), F(2), radii[len(radii) // 2]]
    soft, hard = [], []
    for x in range(9):
        u = [F(1 + (x + i) % 7, 2 + i) for i in range(11)]
        v = [F(1 + (2 * x + i) % 5, 3 + i) for i in range(11)]
        soft.append([q * u[i] / (mass[i] * sum(u)) for i in range(11)])
        hard.append([Ch * q * v[i] / (born[i] * sum(v)) for i in range(11)])
        ck(f'n{n}:whole_positive_rows_{x}', sum(mass[i] * soft[x][i] for i in range(11)) == q
           and sum(born[i] * hard[x][i] for i in range(11)) == Ch * q)
    norm = max(sum(hard[x][i] for x in range(9)) for i in range(11))
    rx = 1 / norm
    total_bad = F(0)
    child_positive = F(0)
    for x in range(9):
        rad = R[x % len(R)]
        bad_x, good_x = F(0), F(0)
        for y in range(11):
            for z in range(11):
                if y == z:
                    continue
                eligible = [p for p in path(z) if p in certs and rad <= certs[p][0]]
                if not eligible:
                    continue
                P = eligible[0]
                G = F(1, 2 + (x + y + z) % 11)
                value = mass[y] * soft[x][y] * born[z] * hard[x][z] * G / q
                if z in certs[P][1]:
                    good_x += value
                else:
                    bad_x += value
                if y == 1 and z == 0 and P == 2:
                    child_positive += rx * value
        ck(f'n{n}:source_selected_good_row_{x}', good_x <= Ch * q)
        ck(f'n{n}:one_fixed_bad_marginal_{x}',
           bad_x <= sum(mass[i] * hard[x][i] for i in bad))
        total_bad += rx * bad_x
    for i in range(11):
        ck(f'n{n}:finite_single_source_column_{i}', rx * sum(hard[x][i] for x in range(9)) <= 1)
    ck(f'n{n}:fixed_union_single_column_fee', total_bad <= sum(mass[i] for i in bad))
    ck(f'n{n}:hard_child_only_positive_traffic', child_positive > 0)
    return {'bad_traffic': fs(total_bad), 'column_constant': '1/1',
            'strict_new_hard_child_traffic': fs(child_positive),
            'scope': 'independent postulated positive rows, not original kernels'}


def compressed(e):
    n, sqrt, k, D = 1 << e, 1 << (e // 2), e + 1, e + 1
    T = sqrt * k ** 4
    # Rh^n = b^n/2^d and h=M; cancel b^n everywhere. No 2^n
    # integer, geometry or enormous-dimensional source is materialized.
    d = (2 * n * k ** 4).bit_length() + 1
    old_ratio, adaptive_ratio = F(2), F(2 * D, 1 << d)
    ck(f'n2^{e}:old_volume_rejected', old_ratio > 1)
    ck(f'n2^{e}:new_weighted_volume_accepted', adaptive_ratio < 1)
    ck(f'n2^{e}:CP_necessary_ratio', (1 << d) > 2 * T * sqrt)
    ck(f'n2^{e}:GP_necessary_power', ((2 * T) ** 2).bit_length() <= n + d)
    ck(f'n2^{e}:above_GP_easy_beta', (T ** 2).bit_length() <= n - d)
    ck(f'n2^{e}:physical_radius_range', 0 < d < n)
    # Old and new SOURCE REGIONS are disjoint: do not add their unit path
    # charges as though both occurred on the same source path.
    ck(f'n2^{e}:path_charge_compression',
       F(2, 5) + F(3, 5) * D * F(1, D) == 1)
    improves_n_half = 2 * T < n
    if e >= 64:
        ck(f'n2^{e}:large_n_sublinear_budget', improves_n_half)
    else:
        ck(f'n2^{e}:small_asymptotic_warning_preserved', not improves_n_half)
    return {'n': str(n), 'sqrt_n': str(sqrt), 'k_n': k, 'D': D,
            'T': str(T), 'T_over_n': fs(F(T, n)), 'T_less_n_over_2': improves_n_half,
            'volume_log2_gap_d': d, 'old_b_normalized_ratio': fs(old_ratio),
            'new_adaptive_normalized_ratio': fs(adaptive_ratio),
            'CP_GP_necessary_scalar_compatibility': True,
            'above_GP_easy_beta': True,
            'scope': 'compressed scalar certificates only; no actual input or full gate contract'}


def main():
    rounds = []
    for n, exponent in ((4, 40), (16, 64), (64, 80)):
        start = len(checks)
        mass, born, certs, bad, radii, record = geometry(n)
        record['positive_marginals'] = positive_arrays(n, mass, born, certs, bad, radii)
        record['compressed_large_dimension'] = compressed(exponent)
        record['checks'] = len(checks) - start
        rounds.append(record)
    payload = {'status': 'PASS', 'seed': None, 'arithmetic': 'fractions.Fraction exact',
               'registration': 'adaptive_core_source_budget_20261007.md section 8',
               'scope': 'New source-contract/tree/closed-core/positive-marginal and compressed-ratio components. No original phi, FIRST, actual coverage proportion, actual CP/GP/history model or full residual counterexample.',
               'total_checks': len(checks), 'rounds': rounds,
               'checks_sha256': hashlib.sha256(json.dumps(checks, sort_keys=True).encode()).hexdigest(),
               'checks': checks}
    OUT.write_text(json.dumps(payload, ensure_ascii=False, indent=2) + '\n')
    print(json.dumps({'status': 'PASS', 'total_checks': len(checks),
                      'rounds': [{'n': r['n'], 'checks': r['checks'],
                                  'new_nodes': r['new_certificate_nodes']} for r in rounds],
                      'path': str(OUT)}, ensure_ascii=False))


if __name__ == '__main__':
    main()
