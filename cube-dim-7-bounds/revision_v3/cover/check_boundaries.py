"""Small exact audit checks, not a parameter search or proof substitute."""
from fractions import Fraction as F
from itertools import product
from math import prod
import json


def inside(x, cube, dilation=F(1)):
    center, radius = cube
    return all(abs(v - c) < dilation * radius for v, c in zip(x, center))


def volume_union(cubes):
    if not cubes:
        return F(0)
    grids = [sorted({c[i] + sign * r for c, r in cubes for sign in (-1, 1)})
             for i in range(len(cubes[0][0]))]
    intervals = [list(zip(g, g[1:])) for g in grids]
    result = F(0)
    for box in product(*intervals):
        x = tuple((lo + hi) / 2 for lo, hi in box)
        if any(inside(x, cube) for cube in cubes):
            result += prod(hi - lo for lo, hi in box)
    return result


def finite_selection():
    # Equal radii, an exact duplicate (zero increment), and genuine rejection.
    cubes = [((F(x), F(0)), F(r)) for x, r in
             [(0, 1), (0, 1), (15, 1), (7, F(1, 2)), (-3, F(1, 2))]]
    n, a, t, tau = 2, F(9), F(2), F(1, 50)
    b = a * (1 + 1 / t)
    selected, weights, increments = [], [], []
    remaining = list(range(len(cubes)))
    lost = {}
    while True:
        eligible = [i for i in remaining
                    if sum(w for j, w in zip(selected, weights)
                           if inside(cubes[i][0], cubes[j], b)) <= tau]
        if not eligible:
            break
        k = max(eligible, key=lambda i: (cubes[i][1], -i))
        old = volume_union([cubes[j] for j in selected])
        new = volume_union([cubes[j] for j in selected] + [cubes[k]])
        v = new - old
        selected.append(k)
        increments.append(v)
        weights.append(v / (2 * a * cubes[k][1]) ** n)
        remaining.remove(k)
        for i in remaining:
            load = sum(w for j, w in zip(selected, weights)
                       if inside(cubes[i][0], cubes[j], b))
            if load > tau and i not in lost:
                lost[i] = len(selected)
    radii = [cubes[i][1] for i in selected]
    assert all(x >= y for x, y in zip(radii, radii[1:]))
    assert selected == [0, 1, 2, 4] and increments[1] == 0
    assert remaining == [3] and lost[3] == 3
    for i in remaining:
        prefix = lost[i]
        assert all(cubes[j][1] >= cubes[i][1] for j in selected[:prefix])
        paying = [j for j in selected[:prefix] if inside(cubes[i][0], cubes[j], b)]
        assert sum(weights[selected.index(j)] for j in paying) > tau
        for j in paying:
            assert all(abs(v - c) + cubes[i][1] <= (b + 1) * cubes[j][1]
                       for v, c in zip(cubes[i][0], cubes[j][0]))
    assert sum(increments) == volume_union([cubes[j] for j in selected])
    A = ((b + 1) / a) ** n
    assert volume_union(cubes) <= (1 + A / tau) * sum(increments)
    # All cells of the selected-star arrangement: tie handling and unpaid radii.
    stars = [(cubes[j][0], a * cubes[j][1]) for j in selected]
    grids = [sorted({c[i] + sign * r for c, r in stars for sign in (-1, 1)})
             for i in range(n)]
    checked, tied, unpaid = 0, 0, 0
    for box in product(*[list(zip(g, g[1:])) for g in grids]):
        x = tuple((lo + hi) / 2 for lo, hi in box)
        contributors = [k for k, (j, w) in enumerate(zip(selected, weights))
                        if w > 0 and inside(x, cubes[j], a)]
        if not contributors:
            continue
        rmin = min(cubes[selected[k]][1] for k in contributors)
        same = [k for k in contributors if cubes[selected[k]][1] == rmin]
        k = max(same)
        tied += len(same) > 1
        assert all(j < k for j in contributors if j != k)
        captured = [j for j in contributors if j != k and
                    inside(cubes[selected[k]][0], cubes[selected[j]], b)]
        assert sum(weights[j] for j in captured) <= tau
        assert weights[k] <= a ** (-n)
        for j in contributors:
            if j != k and j not in captured:
                ratio = cubes[selected[j]][1] / cubes[selected[k]][1]
                assert 1 <= ratio < t
                unpaid += 1
        checked += 1
    assert tied > 0 and unpaid > 0
    return {"selected_indices": selected, "rejected_indices": remaining,
            "increments": list(map(str, increments)),
            "star_cells_checked": checked, "tie_cells": tied,
            "unpaid_contributions_checked": unpaid}


def face_checks():
    c, rho = F(7, 5), F(3, 2)
    # At u_i=3/5 an entire negative face is on the equality boundary.
    vectors = [(F(3, 5),), (F(3, 5), F(4, 5)),
               (F(0), F(1, 5), F(3, 5), F(4, 5), F(2, 5), F(3, 5), F(9, 10))]
    result = []
    for u in vectors:
        n = len(u)
        direct = F(0)
        for i in range(n):
            other_lengths = []
            for j in range(n):
                if j != i:
                    lo, hi = max(-rho, u[j] - c * rho), min(rho, u[j] + c * rho)
                    other_lengths.append(max(F(0), hi - lo))
            for sign in (-1, 1):
                if abs(sign * rho - u[i]) < c * rho:
                    direct += prod(other_lengths)
        L = [min(2 * rho, (1 + c) * rho - v) for v in u]
        formula = sum((1 + int(v < (c - 1) * rho)) *
                      prod(L[j] for j in range(n) if j != i)
                      for i, v in enumerate(u))
        assert direct == formula
        # Increasing a coordinate must not increase captured face measure.
        u2 = [min(F(19, 20), v + F(1, 40)) for v in u]
        L2 = [min(2 * rho, (1 + c) * rho - v) for v in u2]
        formula2 = sum((1 + int(v < (c - 1) * rho)) *
                       prod(L2[j] for j in range(n) if j != i)
                       for i, v in enumerate(u2))
        assert formula2 <= formula
        result.append({"n": n, "captured_face_measure": str(formula),
                       "after_coordinate_increase": str(formula2)})
    return result


def finite_normalization_checks():
    t, a = F(7, 3), F(11)
    alpha, b = a / (a + 1), a * (1 + 1 / t)
    c = (b - 1) / (a + 1)
    assert a > 2 * t and c > 1 and alpha / (c - 1) > t
    result = []
    for rho in [F(1), F(17, 9), t]:
        ell = 1 + 1 / (2 * t) - 1 / (2 * rho)
        ell_a = (1 + c - alpha / rho) / 2
        assert ell_a == alpha * ell and 0 < ell_a <= alpha < 1
        for n in (1, 2, 7):
            L = 2 * rho * ell_a
            raw = F(n, 2 ** n) * rho ** (-n) * L ** (n - 1)
            normalized = F(n, 2) * ell_a ** (n - 1) / rho
            assert raw == normalized
        result.append({"rho": str(rho), "ell_a": str(ell_a)})
    # n=1 direct cap support: thresholds at rho=1 and rho=t, empty products.
    d, t1 = F(2, 5), F(9, 4)
    cases = []
    for u in [F(0), d, F(3, 5), d * t1, F(19, 20)]:
        inner_negative_length = max(F(0), u - d)
        cutoff = max(F(1), min(t1, u / d))
        assert inner_negative_length == max(F(0), u - (1 + d) + 1)
        # H_1 = inner_negative_length/2 + log(cutoff)/2.
        cases.append({"u": str(u), "inner_negative_length": str(inner_negative_length),
                      "shell_log_argument": str(cutoff)})
    return {"corner_samples": result, "n1_support_checks": cases}


if __name__ == "__main__":
    report = {"status": "all exact checks passed", "arithmetic": "Python fractions.Fraction",
              "finite_selection": finite_selection(), "faces": face_checks(),
              "normalization": finite_normalization_checks(),
              "scope": "fixed finite fixtures only; no optimizer, scan, scalar certificate or theorem proof"}
    print(json.dumps(report, ensure_ascii=False, indent=2))
