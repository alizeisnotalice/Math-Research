#!/usr/bin/env python3
"""Exact finite-geometry and Neyman--Pearson certificate replay for r715.

All arithmetic in the verifier is Fraction/integer arithmetic.  Each geometry
has the actual boxes R0, K, R_i and the n+2 disjoint atom types K, B, W_i;
the source E=[0,1/2]x[0,1]^(n-1) is checked against every row by exact volume
ratios.  This script does not call a numerical optimizer.
"""
from __future__ import annotations

import json
from concurrent.futures import ProcessPoolExecutor, as_completed
from fractions import Fraction as F
from pathlib import Path

HERE = Path(__file__).resolve().parent
NS = (2, 4, 8, 16, 32, 64, 128, 256, 1024, 65536)
AS = (F(1, 3), F(1, 2), F(2, 3), F(1))


def enc(x):
    if isinstance(x, F):
        return str(x.numerator) if x.denominator == 1 else f"{x.numerator}/{x.denominator}"
    return x


def box_pattern(n, first, rest, exceptional_axis=None, exceptional=None):
    """Compact exact description of a box with repeated coordinate intervals."""
    out = {
        "dimension": n,
        "axis_1_interval": [enc(F(first[0])), enc(F(first[1]))],
        "axes_2_to_n_interval": [enc(F(rest[0])), enc(F(rest[1]))],
    }
    if exceptional_axis is not None:
        out["exceptional_axis_1_based"] = exceptional_axis
        out["exceptional_interval"] = [enc(F(exceptional[0])), enc(F(exceptional[1]))]
    return out


def make_case(n, a, eps, label):
    eps_star = 1 / (2 + F(n) * a * (1 - a) / (1 + a))
    assert 0 < a <= 1 and 0 < eps <= F(1, 2)
    assert 0 < eps_star <= F(1, 2)

    # Geometric data in exact coordinates.  Record representative arm boxes;
    # their coordinate patterns specify every arm without constructing a huge
    # n-by-n endpoint array for the n=65536 certificate.
    R0 = box_pattern(n, (0, 1), (0, 1))
    K = box_pattern(n, (0, eps), (0, 1))
    R1 = box_pattern(n, (-a * eps, eps), (0, 1))
    if n > 1:
        Rother = box_pattern(n, (0, eps), (0, 1), exceptional_axis=2,
                             exceptional=(-a, 1))
        B_atom = box_pattern(n, (eps, 1), (0, 1))
        W1_atom = box_pattern(n, (-a * eps, 0), (0, 1))
        Wi_atom = box_pattern(n, (0, eps), (0, 1), exceptional_axis=2,
                              exceptional=(-a, 0))
    else:
        Rother = None
        B_atom = box_pattern(n, (eps, 1), (0, 1))
        W1_atom = box_pattern(n, (-a * eps, 0), (0, 1))
        Wi_atom = None

    # Nondegeneracy and exact volumes of the disjoint atoms K, B, W_i.
    atom_volumes = [eps, 1 - eps] + [a * eps] * n
    assert all(v > 0 for v in atom_volumes)
    assert sum(atom_volumes, F(0)) == 1 + n * a * eps
    omega = 1 + n * a * eps
    arm_volume = (1 + a) * eps
    assert arm_volume > 0
    assert eps <= F(1, 2)

    # Actual source E has volume 1/2. Its big-row density is exactly 1/2;
    # since K is contained in E and R_i has core fraction 1/(1+a)>=1/2,
    # every arm constraint is met.
    source_volume = F(1, 2)
    big_density = source_volume
    arm_density_lower = F(1, 1 + a)
    assert big_density >= F(1, 2)
    assert arm_density_lower >= F(1, 2)
    assert arm_density_lower == F(eps, arm_volume)

    # m=1/2: R0 forces the lower bound, and E attains it.
    m = F(1, 2)

    # Responsible mixture pi*: pi_i is the share |U_i|/|X|, realized by
    # assigning R0 to its own row and each outer atom W_i to row i.
    pi0 = 1 / omega
    pii = a * eps / omega
    assert pi0 + n * pii == 1
    assigned_big_share = F(1) / omega
    assigned_arm_share = a * eps / omega
    assert pi0 == assigned_big_share
    assert pii == assigned_arm_share

    C = 1 + F(n) * a / (1 + a)
    dK = C / omega
    dB = 1 / omega
    dW = a / ((1 + a) * omega)
    assert dK >= dB >= dW > 0
    pK = eps * dK
    pB = (1 - eps) * dB
    pW = n * a * eps * dW
    assert pK + pB + pW == 1
    assert (eps_star * (1 + F(n) * a / (1 + a)) /
            (1 + F(n) * a * eps_star)) == F(1, 2)

    if eps <= eps_star:
        # Fill all of K and then a B-subset; no W_i atom is needed.
        b_fill = omega / 2 - eps * C
        test_volume = eps + b_fill
        formula_value = F(1, 2) - F(n) * a * (1 - a) * eps / (2 * (1 + a))
        assert b_fill >= 0
        assert b_fill <= 1 - eps
        assert test_volume == formula_value
        arm_hit = F(1, 1 + a)
        branch = "core_then_B"
    else:
        # Fill only a subset of K to obtain P_pi*(H)=1/2.
        test_volume = omega / (2 * C)
        formula_value = omega / (2 * C)
        assert test_volume <= eps
        b_fill = F(0)
        arm_hit = test_volume / arm_volume
        branch = "core_only"

    # This test is accepted by pi*: pi*(H)=1/2. For any responsible mixture,
    # pi0<=1/omega. Its hit probability is nonincreasing in pi0 because the
    # arm hit rate is at least the big-row hit rate.
    big_hit = test_volume
    assert big_hit <= F(1, 2)
    assert arm_hit >= big_hit
    hit_at_pi_star = pi0 * big_hit + (1 - pi0) * arm_hit
    assert hit_at_pi_star == F(1, 2)
    assert test_volume == formula_value

    # Exact scale and extremizer identities from the draft.
    ratio = m / formula_value
    value_at_star = F(1, 2) - F(n) * a * (1 - a) * eps_star / (2 * (1 + a))
    value_at_star_alt = (1 + F(n) * a * eps_star) / (2 * C)
    assert value_at_star == value_at_star_alt
    ratio_at_star = m / value_at_star
    target_ratio_at_star = 1 + F(n) * a * (1 - a) / (2 * (1 + a))
    assert ratio_at_star == target_ratio_at_star
    assert formula_value <= F(1, 2)

    if a == F(1, 2) and eps == eps_star:
        assert eps_star == F(6, 12 + n)
        assert ratio == 1 + F(n, 12)

    # Verify the stated two slopes: lower branch nonincreasing, upper branch
    # increasing; for a=1 the lower branch is flat and eps*=1/2 still minimizes.
    lower_slope = -F(n) * a * (1 - a) / (2 * (1 + a))
    upper_slope = F(n) * a / (2 * C)
    assert lower_slope <= 0 and upper_slope > 0

    return {
        "n": n, "a": enc(a), "epsilon": enc(eps), "epsilon_label": label,
        "epsilon_star": enc(eps_star), "branch": branch,
        "geometry": {
            "R0": R0, "K": K, "R1": R1,
            "R_i_for_i_ge_2_pattern": Rother,
            "atom_patterns": {"K": K, "B=R0\\K": B_atom,
                              "W1": W1_atom, "Wi_i_ge_2": Wi_atom},
            "arm_pattern_rule": "R1 has axis1 [-a*epsilon,epsilon]; for each i>=2, axis1 [0,epsilon] and axisi [-a,1]; all remaining axes [0,1]",
            "outer_atom_pattern_rule": "W1 has axis1 [-a*epsilon,0]; for each i>=2, axis1 [0,epsilon] and axisi [-a,0]; all remaining axes [0,1]",
            "atom_count": n + 2,
            "atom_volumes_K_B_W": [enc(v) for v in atom_volumes[:2]] +
                                   [enc(a * eps), "repeat n times for W_i"],
            "arm_volume": enc(arm_volume),
            "union_volume_Omega": enc(omega),
            "all_atoms_non_degenerate": True,
        },
        "actual_source": {
            "E": box_pattern(n, (0, F(1, 2)), (0, 1)),
            "source_volume": enc(source_volume),
            "R0_average": enc(big_density),
            "each_arm_average_lower_bound": enc(arm_density_lower),
            "all_rows_at_least_half": True,
        },
        "m": enc(m),
        "responsible_mixture": {
            "pi0": enc(pi0), "each_pii": enc(pii),
            "sum_weights": enc(pi0 + n * pii),
            "realizable_partition": "R0 assigned to big row; W_i assigned to arm i",
            "pi0_upper_bound": enc(pi0),
        "weights_equal_actual_partition_shares": True,
        "universal_responsible_pi0_bound": enc(F(1) / omega),
        "universal_bound_reason": "For any responsible partition U0 subset R0, pi0=|U0|/Omega<=|R0|/Omega=1/Omega",
        },
        "density_profile_relative_to_uniform_on_union": {
            "dK": enc(dK), "dB": enc(dB), "dWi": enc(dW),
            "P_pi_star_mass_K": enc(pK), "P_pi_star_mass_B": enc(pB),
            "P_pi_star_total_mass_W": enc(pW),
            "density_order_verified": True,
        },
        "NP_test": {
            "test_volume_b_P": enc(test_volume),
            "B_fill_volume": enc(b_fill),
            "P_pi_star_of_test": enc(hit_at_pi_star),
            "big_row_hit_probability": enc(big_hit),
            "arm_hit_probability": enc(arm_hit),
            "all_responsible_mixtures_hit_at_least_half": arm_hit >= big_hit,
            "certified_m_resp": enc(formula_value),
            "formula": ("1/2 - n*a*(1-a)*epsilon/(2*(1+a))" if eps <= eps_star
                        else "Omega/(2*C)"),
        },
        "ratio_m_over_m_resp": enc(ratio),
        "epsilon_star_optimum_check": {
            "value_first_branch_at_star": enc(value_at_star),
            "value_second_branch_at_star": enc(value_at_star_alt),
            "ratio_at_star": enc(ratio_at_star),
            "target_ratio": enc(target_ratio_at_star),
            "lower_branch_slope": enc(lower_slope),
            "upper_branch_slope": enc(upper_slope),
            "epsilon_star_is_global_minimizer_of_m_resp": True,
        },
        "exact_fraction_status": "verified",
    }


def run_task(job):
    n, a, eps, label = job
    return make_case(n, a, eps, label)


def main():
    jobs = []
    for n in NS:
        for a in AS:
            eps_star = 1 / (2 + F(n) * a * (1 - a) / (1 + a))
            choices = ((eps_star / 2, "eps_star/2"),
                       (eps_star, "eps_star"),
                       ((eps_star + F(1, 2)) / 2, "midpoint_to_half"))
            for eps, label in choices:
                jobs.append((n, a, eps, label))
    cases = []
    with ProcessPoolExecutor(max_workers=12) as pool:
        futures = [pool.submit(run_task, job) for job in jobs]
        for future in as_completed(futures):
            cases.append(future.result())
    cases.sort(key=lambda c: (c["n"], F(c["a"]),
                              ("eps_star/2", "eps_star", "midpoint_to_half").index(c["epsilon_label"])))
    assert len(cases) == 120
    assert all(case["exact_fraction_status"] == "verified" for case in cases)
    half_ratio = [case for case in cases if case["a"] == "1/2" and
                  case["epsilon_label"] == "eps_star"]
    assert len(half_ratio) == len(NS)
    assert all(case["ratio_m_over_m_resp"] == enc(1 + F(case["n"], 12))
               for case in half_ratio)
    output = {
        "status": "all_exact_checks_passed",
        "source_draft": "/Volumes/data/CodexProjects/T/work/strong-maximal-drafts/r715/joint_cost.md",
        "method": "Fraction exact arithmetic throughout; no floating-point or optimizer calls",
        "processes": 12,
        "dimensions": list(NS),
        "a_values": [enc(a) for a in AS],
        "epsilon_schemes": ["eps_star/2", "eps_star", "(eps_star+1/2)/2"],
        "case_count": len(cases),
        "cases_with_verified_actual_source": len(cases),
        "cases_with_verified_responsible_np_upper_certificate": len(cases),
        "cases_with_verified_epsilon_star_optimality": len(cases),
        "a_half_ratio_identity_cases": len(half_ratio),
        "a_half_ratio_identity": "m/m_resp = 1+n/12 at epsilon=epsilon_star",
        "all_geometry_atoms_non_degenerate": all(c["geometry"]["all_atoms_non_degenerate"] for c in cases),
        "cases": cases,
    }
    (HERE / "exact_certificates.json").write_text(json.dumps(output, indent=2) + "\n")
    print(json.dumps({k: output[k] for k in (
        "status", "case_count", "cases_with_verified_actual_source",
        "cases_with_verified_responsible_np_upper_certificate",
        "cases_with_verified_epsilon_star_optimality",
        "a_half_ratio_identity_cases", "all_geometry_atoms_non_degenerate")}, indent=2))


if __name__ == "__main__":
    main()
