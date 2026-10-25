#!/usr/bin/env python3
"""One finite four-subset / fractional-Holder negative determinant moment batch."""
from collections import Counter, defaultdict
from fractions import Fraction as F
from itertools import combinations, permutations, product
from math import gcd, prod
from pathlib import Path
import hashlib
import json
import sys
import time

import flint
from flint import arb, ctx

HERE = Path(__file__).resolve().parent
PERM4 = [(p, (-1)**sum(p[i] > p[j] for i in range(4) for j in range(i + 1, 4)))
         for p in permutations(range(4))]


def det4(rows):
    return sum(sign * prod(rows[i][p[i]] for i in range(4)) for p, sign in PERM4)


def det_fraction(rows):
    a = [list(map(F, row)) for row in rows]
    answer = F(1)
    for i in range(len(a)):
        j = next((j for j in range(i, len(a)) if a[j][i]), None)
        if j is None:
            return F(0)
        if i != j:
            a[i], a[j] = a[j], a[i]
            answer = -answer
        v = a[i][i]
        answer *= v
        for j in range(i + 1, len(a)):
            factor = a[j][i] / v
            for k in range(i + 1, len(a)):
                a[j][k] -= factor * a[i][k]
            a[j][i] = 0
    return answer


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def ar(x):
    x = F(x)
    return arb(x.numerator) / x.denominator


def parse_ball(s):
    assert s.startswith("[") and s.endswith("]")
    c, r = map(F, s[1:-1].split(" +/- "))
    return c - r, c + r


def rational_cap(x, denominator=10**13):
    upper = F(str(x.upper().fmpq()))
    return F((upper * denominator).__ceil__(), denominator)


def main():
    began = time.monotonic()
    registration = json.loads((HERE / "registration.json").read_text())
    for source in registration["inputs"]:
        assert sha(Path(source["source"])) == source["sha256"], source["source"]
    source_row = next(x for x in registration["inputs"] if x["path"].endswith("round006/orbit_minkowski.json"))
    baseline = json.loads(Path(source_row["source"]).read_text())
    old_J = parse_ball(baseline["J_upper_ball"])
    old_I = parse_ball(baseline["I_upper_ball"])

    V = [(1,) + x for x in product((-1, 1), repeat=4)]
    frame = [[sum(v[i] * v[j] for v in V) for j in range(5)] for i in range(5)]
    assert frame == [[16 * int(i == j) for j in range(5)] for i in range(5)]
    W, denominator = 16**5, 5 * 16**5
    groups = {}
    occupancy = [0] * 16
    zero_rank, nonzero_rank, total_C = 0, 0, 0
    example = None
    stream_hash = hashlib.sha256()
    for J in combinations(range(16), 4):
        rows = [V[j] for j in J]
        normal = tuple((-1)**k * det4([row[:k] + row[k+1:] for row in rows]) for k in range(5))
        sq = sum(x*x for x in normal)
        if sq == 0:
            zero_rank += 1
            continue
        nonzero_rank += 1
        C = 16 * sq
        weights = [(sum(x*y for x, y in zip(normal, v)))**2 for v in V]
        assert sum(weights) == C and all(weights[j] == 0 for j in J)
        scope = set(J) | {i for i, w in enumerate(weights) if w}
        for i in scope:
            occupancy[i] += C
        total_C += C
        common = 0
        for x in normal:
            common = gcd(common, abs(x))
        primitive = tuple(x // common for x in normal)
        key = tuple(sorted(map(abs, primitive), reverse=True))
        distribution = Counter(F(w, C) for w in weights if w)
        if key not in groups:
            groups[key] = {"count": 0, "sum_C": 0, "probability_counts": distribution,
                           "representative_J": J, "primitive_normal": primitive,
                           "scope_size": len(scope)}
        group = groups[key]
        assert group["probability_counts"] == distribution
        assert group["scope_size"] == len(scope)
        group["count"] += 1
        group["sum_C"] += C
        stream_hash.update(json.dumps([J, normal, C, sorted(scope)], separators=(",", ":")).encode())
        if example is None:
            example = (J, normal, weights, C)

    assert zero_rank + nonzero_rank == 1820
    assert total_C == denominator
    assert len(set(occupancy)) == 1
    rho = F(occupancy[0], denominator)
    assert 0 < rho <= 1
    assert sum(F(g["sum_C"], denominator) for g in groups.values()) == 1
    assert sum(F(g["sum_C"] * g["scope_size"], 16 * denominator) for g in groups.values()) == rho

    # One independent completion example: exact 5x5 elimination vs cofactor dots.
    J, normal, weights, C = example
    for i in range(16):
        assert det_fraction([V[j] for j in J] + [V[i]])**2 == weights[i]
    assert sum(weights) == C
    # Equal S=1: each grouped normalized polynomial is 1; weighted AMGM is exact.
    assert total_C // 5 == W

    beta, b, q = F(3, 4), F(3, 20), F(10, 7)
    assert 5*b == beta and b/beta == F(1, 5)
    ctx.prec = 384
    product_R = arb(1)
    type_rows = []
    for key, group in sorted(groups.items()):
        theta = F(group["sum_C"], denominator)
        R = sum((count * ar(p)**ar(beta) for p, count in group["probability_counts"].items()), arb(0))
        factor = R**ar(-b * theta / beta)
        product_R *= factor
        type_rows.append({"normal_type_abs_sorted": key, "four_subset_count": group["count"],
                          "sum_completion_weight_C": group["sum_C"], "theta": str(theta),
                          "nonzero_completion_directions": sum(group["probability_counts"].values()),
                          "scope_size": group["scope_size"],
                          "probability_counts": {str(p): count for p, count in sorted(group["probability_counts"].items())},
                          "representative_J": group["representative_J"],
                          "primitive_normal": group["primitive_normal"], "R_ball": R.str(85)})
    common = 2 * ar(F(7, 4)).gamma()
    moment_b = ar(1 + b/beta).gamma() / ar(1 + b).gamma()
    moment_brho = ar(1 + b*rho/beta).gamma() / ar(1 + b*rho).gamma()
    J_ordinary = common * moment_b**5 * product_R
    J_finner = common * moment_brho**ar(F(5) / rho) * product_R
    pi = arb.pi()
    c = 3 / (ar(2)**ar(F(15, 2)) * pi**ar(F(5, 2)))
    I_ordinary = c**ar(q) * pi**2 / 60 * J_ordinary**ar(q)
    I_finner = c**ar(q) * pi**2 / 60 * J_finner**ar(q)
    cap = rational_cap(I_finner)
    assert I_finner < ar(cap)
    old_low_J, old_low_I = map(ar, [old_J[0], old_I[0]])
    checks = {"ordinary_J_improves_old": bool(J_ordinary < old_low_J),
              "finner_J_improves_old": bool(J_finner < old_low_J),
              "finner_I_improves_old": bool(I_finner < old_low_I),
              "finner_beats_ordinary": bool(J_finner < J_ordinary),
              "new_I_below_target_11e_minus6": bool(I_finner < ar(F(11, 10**6))),
              "new_I_below_1_over_60000": bool(I_finner < ar(F(1, 60000)))}
    out = {"status": "finite_grouped_negative_moment_batch_complete",
           "model": registration["model"], "new_C5_upper_bound": False,
           "source_sha256": sha(Path(__file__)), "registration_sha256": sha(HERE / "registration.json"),
           "source_inputs": registration["inputs"],
           "geometry": {"four_subsets": 1820, "rank_below_four": zero_rank,
                        "rank_four": nonzero_rank, "sum_C": total_C, "expected_5W": denominator,
                        "normal_type_count": len(groups), "normal_types": type_rows,
                        "coordinate_scope_occupancy_rho": str(rho),
                        "fractional_matching_weights": "p_J=theta_J/rho; coordinate sums exactly 1",
                        "geometry_stream_sha256": stream_hash.hexdigest(),
                        "completion_example_J": J, "completion_example_normal": normal,
                        "completion_example_weights": weights,
                        "completion_example_C": C},
           "rational_parameters": {"beta": str(beta), "b": str(b), "q": str(q), "rho": str(rho)},
           "arb_384bit": {"R_type_geometric_factor": product_R.str(85),
                          "ES_minus_b": moment_b.str(85), "ES_minus_b_rho": moment_brho.str(85),
                          "ordinary_Holder_J_upper": J_ordinary.str(85),
                          "fractional_Holder_J_upper": J_finner.str(85),
                          "ordinary_Holder_complete_I_upper": I_ordinary.str(85),
                          "fractional_Holder_complete_I_upper": I_finner.str(85),
                          "new_I_over_locked_old_I_lower": (I_finner / old_low_I).str(85)},
           "strict_rational_I_upper": str(cap), "checks": checks,
           "runtime": {"python": sys.executable, "flint": flint.__version__, "arb_precision_bits": ctx.prec,
                       "seconds": time.monotonic()-began, "basis_5subset_enumeration": False,
                       "spatial_or_probability_samples": 0},
           "scope": "New complete all-time/all-space I bound if the grouped AMGM/Holder proof is accepted; no g claim; not the excluded scalar majorant family"}
    (HERE / "results.json").write_text(json.dumps(out, indent=2) + "\n")
    print(json.dumps({"status": out["status"], "normal_type_count": len(groups),
                      "rho": str(rho), "ordinary_J": J_ordinary.str(25), "finner_J": J_finner.str(25),
                      "ordinary_I": I_ordinary.str(25), "finner_I": I_finner.str(25),
                      "rational_I_upper": str(cap), "checks": checks, "seconds": out["runtime"]["seconds"]}, indent=2))


if __name__ == "__main__":
    main()
