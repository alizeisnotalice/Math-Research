#!/usr/bin/env python3
"""Finite independent interface checks; no QMC, source edits, or parameter search.

Only A and K are enclosed with Arb.  All subsequent roots, budgets,
determinants, decimal endpoints, and comparisons use exact Fraction arithmetic.
The contour, obstacle, and Jensen arguments are audited in report.md, not proved
by this program.  The previously certified quotient polytope volume is an input.
"""
from __future__ import annotations

from fractions import Fraction as F
from itertools import permutations, product
from pathlib import Path
import hashlib
import json
import platform
import re
import sys

import flint
from flint import arb, ctx
import sympy as sp

HERE = Path(__file__).resolve().parent
DIGITS = 90


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def det(rows):
    """Exact row elimination, independent of production permutation expansion."""
    a = [list(map(F, row)) for row in rows]
    ans = F(1)
    for i in range(len(a)):
        j = next((j for j in range(i, len(a)) if a[j][i]), None)
        if j is None:
            return F(0)
        if j != i:
            a[i], a[j] = a[j], a[i]
            ans = -ans
        pivot = a[i][i]
        ans *= pivot
        for j in range(i + 1, len(a)):
            mul = a[j][i] / pivot
            for k in range(i + 1, len(a)):
                a[j][k] -= mul * a[i][k]
            a[j][i] = F(0)
    return ans


def iroot(n, k):
    assert n >= 0 and k >= 1
    lo, hi = 0, 1 << ((n.bit_length() + k - 1) // k)
    while lo + 1 < hi:
        mid = (lo + hi) // 2
        if mid**k <= n:
            lo = mid
        else:
            hi = mid
    return hi if hi**k == n else lo


def root_interval(x, k, digits=DIGITS):
    """Return rational decimal bounds; certification uses integer powers only."""
    assert x >= 0
    scale = 10**digits
    z = iroot((x.numerator * scale**k) // x.denominator, k)
    lo = F(z, scale)
    hi = lo if lo**k == x else F(z + 1, scale)
    assert lo**k <= x <= hi**k
    return lo, hi


def rounded_arb_interval(x, digits=110):
    scale = 10**digits
    lo, hi = F(str(x.lower().fmpq())), F(str(x.upper().fmpq()))
    a = F((lo * scale).__floor__(), scale)
    b = F((hi * scale).__ceil__(), scale)
    assert a <= lo <= hi <= b
    return a, b


def decimal_endpoint(x, upward=False, digits=42):
    scale = 10**digits
    k = (x * scale).__ceil__() if upward else (x * scale).__floor__()
    sign = "-" if k < 0 else ""
    whole, rem = divmod(abs(k), scale)
    return f"{sign}{whole}.{rem:0{digits}d}"


def encoded(interval):
    lo, hi = interval
    return {
        "lower_fraction": str(lo),
        "upper_fraction": str(hi),
        "decimal_outward": [decimal_endpoint(lo), decimal_endpoint(hi, True)],
    }


def printed_ball(s):
    m = re.fullmatch(r"\[([^ ]+) \+/- ([^ ]+)\]", s)
    assert m, s
    c, r = map(F, m.groups())
    return c - r, c + r


def main():
    registration = json.loads((HERE / "registration.json").read_text())
    checked_hashes = []
    resolved = {}
    for entry in registration["inputs"] + registration.get("extra_inputs", []):
        snapshot = HERE / entry["snapshot"]
        # In the primary archive avoid duplicating already imported source files.
        original = HERE.parents[1] / entry["path"]
        selected = snapshot if snapshot.exists() else original
        assert sha(selected) == entry["sha256"], entry["path"]
        resolved[entry["path"]] = selected
        checked_hashes.append({"snapshot": entry["snapshot"], "sha256": entry["sha256"]})

    def inp(rel):
        return resolved["cube-dim-5-bounds/" + rel]

    source = json.loads(inp("round005/alpha15_middle_budget.json").read_text())
    assert source["source_sha256"] == sha(inp("round005/alpha15_middle_budget.py"))
    grow = json.loads(inp("upper_g_round010/tilted_gaussian_bound.json").read_text())
    assert grow["source_sha256"] == sha(inp("upper_g_round010/tilted_gaussian_bound.py"))
    theory = json.loads(inp("upper/round010_alpha_strata_results.json").read_text())
    assert theory["files_sha256"]["round010_alpha_strata.py"] == sha(inp("upper/round010_alpha_strata.py"))
    assert theory["files_sha256"]["round002_quotient_peak_proof.md"] == sha(inp("upper/round002_quotient_peak_proof.md"))

    # E1: all 16 directions, quotient coordinates, and hyperoctahedral symmetry.
    directions = tuple((1,) + e for e in product((-1, 1), repeat=4))
    dset = set(directions)
    base = (1, 1, 1, 1, 1)
    complement = ((1, -1, -1, -1, 1), (1, -1, -1, 1, -1),
                  (1, -1, 1, -1, -1), (1, 1, -1, -1, -1))
    geometry = []
    for v in directions:
        reflected = tuple(tuple(w[j] * v[j] for j in range(5)) for w in complement)
        assert all(w in dset and w != v for w in reflected)
        five = (v,) + reflected
        quotient = tuple(tuple(w[j] - v[j] for j in range(1, 5)) for w in reflected)
        coordinate_map = ((1, 0, 0, 0, 0),) + tuple(
            tuple(-v[j] if k == 0 else int(k == j) for k in range(5))
            for j in range(1, 5))
        assert abs(det(five)) == abs(det(quotient)) == 48
        assert det(coordinate_map) == 1
        geometry.append({"selected_direction": v, "full_det": str(det(five)),
                         "quotient_det": str(det(quotient)), "coordinate_Jacobian": 1})
    group_count = 0
    for perm in permutations(range(5)):
        for signs in product((-1, 1), repeat=5):
            moved = set()
            for v in directions:
                w = tuple(signs[j] * v[perm[j]] for j in range(5))
                moved.add(tuple(w[0] * z for z in w))
            assert moved == dset
            group_count += 1
    frame = [[sum(v[i] * v[j] for v in directions) for j in range(5)] for i in range(5)]
    assert frame == [[16 * int(i == j) for j in range(5)] for i in range(5)]
    quotients = {tuple(w[j] - 1 for j in range(1, 5)) for w in directions if w != base}
    assert quotients == {tuple(-2 * e for e in mask) for mask in product((0, 1), repeat=4) if any(mask)}

    # E2: preselected old gate; exact powers avoid copying the production Arb path.
    alpha, n, q = F(3, 2), F(5), F(10, 7)
    assert q == n / (n - alpha)
    assert (1 - 4 / alpha) * (1 - 1 / q) + 1 / q == F(1, 5)
    assert 1 - n / alpha == -F(7, 3)
    assert 1 / (n / alpha - 1) == F(3, 7)
    ctx.prec = 384
    pi = arb.pi()
    A = rounded_arb_interval(2 * (arb(2) / 3).gamma() / (3 * pi))
    K = rounded_arb_interval((arb(25) / 6) * (arb(10) / (3 * pi)).sqrt())
    eps, T, I0 = F(1, 65536), F(4), F(11, 10**6)
    # 16^(5/2)=1024 and B^10=(5/6)^10(1024/3)^3 K^10 A^12.
    B10 = tuple(F(5, 6)**10 * F(1024, 3)**3 * k**10 * a**12 for a, k in zip(A, K))
    short = (root_interval(B10[0] * eps**12, 10)[0],
             root_interval(B10[1] * eps**12, 10)[1])
    # long^3=((3/7) A^5/48)^3 T^-7.
    long = (root_interval((F(3, 7) * A[0]**5 / 48)**3 / T**7, 3)[0],
            root_interval((F(3, 7) * A[1]**5 / 48)**3 / T**7, 3)[1])
    total_norm = root_interval(I0**7, 10)
    remainder = (total_norm[0] - short[1] - long[1],
                 total_norm[1] - short[0] - long[0])
    assert remainder[0] > 0
    gate = (root_interval(remainder[0]**10, 7)[0],
            root_interval(remainder[1]**10, 7)[1])
    row = next(r for r in source["rows"] if r["epsilon"] == "1/65536" and r["T"] == 4)
    comparisons = {}
    for name, val in [("short_Lq_upper_ball", short), ("long_Lq_upper_ball", long),
                      ("remaining_middle_Lq_limit_ball", remainder), ("middle_integral_limit_ball", gate)]:
        saved = printed_ball(row[name])
        assert saved[0] <= val[0] <= val[1] <= saved[1]
        comparisons[name] = {"independent_enclosure": encoded(val),
                             "contained_in_saved_printed_ball": True}
    safe_gate = F(1061435623, 10**14)
    assert 0 < safe_gate < gate[0]

    # E3: cube volume and face measure; conditional transfer arithmetic.
    face_octants, octants_per_full_face, full_faces = 5 * 2**5, 2**4, 2 * 5
    assert face_octants == full_faces * octants_per_full_face
    assert F(face_octants, 5) == 32
    c_over_g = (q - 1) / q
    assert c_over_g == F(3, 10)
    costs = {}
    for g in [F(1, 10000), F(1, 65000)]:
        c7 = F(10, 3)**10 * I0**7 / g**10
        costs[str(g)] = {"conditional_cost_interval": encoded(root_interval(c7, 7)),
                         "strictly_below_32": c7 < 32**7,
                         "strictly_above_32": c7 > 32**7}
    new_g = F(1, 65000)
    admissible_I = root_interval(32**7 * (F(3, 10) * new_g)**10, 7)

    # E4/E5: exact exponents only; the contour proof is checked in the report.
    quotient_dimension, quotient_terms = F(4), F(15)
    assert 1 + quotient_dimension / alpha == F(11, 3)
    assert quotient_dimension / (2 * alpha) == F(4, 3)
    # E6: independent symbolic full-time integration and exact saved-ball tests.
    t, a = sp.symbols("t a", positive=True)
    time_integral = sp.integrate(t**(-sp.Rational(10, 3)) * sp.exp(-a * t**(-sp.Rational(4, 3))), (t, 0, sp.oo))
    expected = sp.Rational(3, 4) * sp.gamma(sp.Rational(7, 4)) * a**(-sp.Rational(7, 4))
    assert sp.simplify(time_integral - expected) == 0
    beta = F(3, 4)
    assert 1 / beta == F(4, 3)
    assert -(1 / beta) * F(5, 2) == -F(10, 3)
    assert (F(10, 3) - 1) / F(4, 3) == F(7, 4)
    saved_g = printed_ball(grow["green_bound_arb"])
    assert saved_g[0] > new_g
    assert saved_g[1] < F(1, 10000)

    # E7: only the four displayed analytic basis classes, extended by symmetry.
    # This does not repeat the author's enumeration of 1365/4368 all-bases cases.
    basis_classes = {
        1: (((1, 0, 0, 0), (0, 1, 1, 0), (0, 1, 0, 1), (0, 0, 1, 1)), 0, 32),
        2: (((1, 1, 0, 0), (1, 0, 1, 0), (1, 0, 0, 1), (0, 1, 1, 1)), 0, 48),
        3: (((1, 0, 0, 1), (0, 1, 0, 1), (0, 0, 1, 1), (1, 1, 1, 0)), 3, 48),
        4: (((1, 1, 0, 0), (1, 0, 1, 0), (1, 0, 0, 1), (1, 1, 1, 1)), 3, 32),
    }
    basis_by_mask, D_by_mask = {}, {}
    for mask in product((0, 1), repeat=4):
        if not any(mask):
            continue
        columns, selected, expected_det = basis_classes[sum(mask)]
        permutation = next(p for p in permutations(range(4))
                           if tuple(columns[selected][p[j]] for j in range(4)) == mask)
        basis = tuple(tuple(-2 * col[permutation[j]] for j in range(4)) for col in columns)
        assert all(v in quotients for v in basis)
        assert basis[selected] == tuple(-2 * e for e in mask)
        assert abs(det(basis)) == expected_det
        basis_by_mask[str(mask)] = basis
        D_by_mask[mask] = expected_det
    reciprocal_sums = [sum((F(1, D) for mask, D in D_by_mask.items() if mask[j]), F(0))
                       for j in range(4)]
    assert reciprocal_sums == [F(3, 16)] * 4
    rank_exponents = [F(k) - F(5 - k) / alpha for k in range(6)]
    assert rank_exponents == [-F(10, 3), -F(5, 3), F(0), F(5, 3), F(10, 3), F(5)]
    assert 4 / (4 - q) == F(14, 9)
    assert (4 - q) / q == F(9, 5)
    assert 4 / q == F(14, 5)
    PP, CC, rr = sp.symbols("P C r", positive=True)
    split = (PP / CC)**sp.Rational(3, 5) * rr**sp.Rational(3, 2)
    split_integral = CC * rr**(-sp.Rational(5, 2)) * split + sp.Rational(3, 2) * PP * split**(-sp.Rational(2, 3))
    split_expected = sp.Rational(5, 2) * PP**sp.Rational(3, 5) * CC**sp.Rational(2, 5) / rr
    assert sp.simplify(split_integral - split_expected) == 0
    rational_checks = {"K_contour": F(607, 1000), "P_quotient_peak_Jensen_volume": F(1, 125000),
                       "B_short_time": F(27, 100), "D_full_G_face_envelope": F(3)}
    for key, bound in rational_checks.items():
        assert printed_ball(theory["constants"][key])[1] < bound
    table_upper_bounds = [
        ["1.230930e-5", "8.225108e-7", "2.029295e-9"],
        ["2.332176e-6", "6.769547e-8", "5.727199e-11"],
        ["4.418647e-7", "5.579983e-9", "1.619853e-12"],
    ]
    table_keys = ["short_time_Lq_norm_bound", "corner_tube_Lq_norm_bound", "corner_tube_integral_G_to_q_bound"]
    for row, displayed in zip(theory["finite_cutoff_budgets"], table_upper_bounds):
        for key, upper in zip(table_keys, displayed):
            assert printed_ball(row[key])[1] < F(upper)

    out = {
        "status": "pass_conditional_interfaces_and_global_g_benchmark_no_new_C5_bound",
        "source_sha256": sha(Path(__file__)),
        "registration_sha256": sha(HERE / "registration.json"),
        "runtime": {"python": sys.executable, "python_version": platform.python_version(),
                    "flint_version": flint.__version__, "sympy_version": sp.__version__,
                    "Arb_bits_for_A_K_only": ctx.prec, "Fraction_root_decimal_digits": DIGITS},
        "checked_inputs": checked_hashes,
        "E1_geometry": {"direction_count": len(directions), "all_16_det48_complements": geometry,
                        "signed_permutation_cases": group_count, "VVt": frame,
                        "quotient_set": "all 15 nonzero {-2,0}^4 vectors", "density_coordinate_Jacobian": 1},
        "E2_old_middle_gate": {"epsilon": str(eps), "T": str(T), "target_total_I": str(I0),
                               "A_input": encoded(A), "K_input": encoded(K),
                               "power_root_method": "integer-power comparisons with exact rational endpoints",
                               "comparisons": comparisons, "safe_rational_middle_integral_gate": str(safe_gate),
                               "middle_integral_itself_certified": False},
        "E3_transfer": {"positive_face_patches": face_octants, "full_faces": full_faces,
                        "octants_per_face": octants_per_full_face, "radial_Jacobian": "r^4",
                        "volume_V1_over_I": 32, "optimal_c_over_g": str(c_over_g),
                        "conditional_costs_at_I0": costs,
                        "I_limit_for_cost_32_at_new_g": encoded(admissible_I)},
        "E4_contour": {"analytic_review": "pass; see report.md", "theta": "pi/3",
                       "subtracted_constant_integral": "i/x, with zero real part",
                       "K_star": "sqrt(2/pi)*3^(-1/4)",
                       "numerical_K_or_new_budget_recomputed": False},
        "E5_quotient_spectral_peak": {"power_mean_terms": 15,
                                     "old_quotient_volume_input": "14455/191102976",
                                     "Gamma_parameter": "11/3", "15_exponent": "4/3",
                                     "Fourier_factor": "(2*pi)^(-4)",
                                     "polytope_volume_recomputed": False},
        "E6_global_g": {"conditional_covariance": "2*A", "Gaussian_factor": "(4*pi)^(-5/2)",
                        "time_integral_by_sympy": str(time_integral),
                        "saved_ball_lower_exceeds_1_over_65000": True,
                        "saved_ball_upper_below_1_over_10000": True,
                        "new_Green_g": str(new_g), "author_Arb_computation_rerun": False},
        "E7_rank_and_tube": {"analytic_review": "pass; see report.md",
                             "constructed_basis_by_selected_mask": basis_by_mask,
                             "constructed_basis_count": len(basis_by_mask),
                             "coordinate_reciprocal_sums": list(map(str, reciprocal_sums)),
                             "rank_exponents_k_0_to_5": list(map(str, rank_exponents)),
                             "weighted_quotient_time_split_by_sympy": str(sp.simplify(split_integral)),
                             "tube_layer_coefficient": "14/9", "tube_norm_powers": ["9/5", "14/5"],
                             "author_display_table_9_upper_roundings_checked": True,
                             "author_full_enumerations_or_constants_recomputed": False},
        "certified_C5_upper": 32,
        "missing_for_target": ["global essential g>=1/10000", "certified matching middle integral on [epsilon,T]"],
    }
    (HERE / "check_results.json").write_text(json.dumps(out, ensure_ascii=False, indent=2) + "\n")
    print(json.dumps({"status": out["status"], "16_direction_checks": len(geometry),
                      "symmetry_checks": group_count, "safe_old_middle_gate": str(safe_gate),
                      "conditional_costs": costs, "sympy_time_integral": str(time_integral)}, indent=2))


if __name__ == "__main__":
    main()
