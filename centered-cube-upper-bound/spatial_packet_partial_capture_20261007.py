"""Nine newly registered exact partial-cell continuous-winner inputs.

Standalone Fraction construction; no old inputs, oracle or MC imported.
"""
from fractions import Fraction as F
from pathlib import Path
import argparse
import datetime
import hashlib
import json
import time

BASE = Path(__file__).resolve().parent
PREFIX = "spatial_packet_partial_capture_20261007"
REG = BASE / (PREFIX + "_registration.json")


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def canonical(obj):
    return json.dumps(obj, sort_keys=True, separators=(",", ":")).encode()


def save(path, obj):
    with Path(path).open("x") as handle:
        json.dump(obj, handle, sort_keys=True, indent=2)
        handle.write("\n")


def make_case(round_id, n, a, b):
    d = min(a / (8 * n), (b - a) / 4)
    eps = n * d / (10 * a)
    source = {"n": n, "atoms": [list(map(str, [9 * d / 10] * n)),
                                  list(map(str, [d / 10] * n))],
              "weights": ["1", str(eps)], "W": str(1 + eps)}
    x = [a / 2 + d / 2] + [d / 2] * (n - 1)
    return {"case_id": f"r{round_id}_n{n}_a{a}_b{b}".replace("/", "d"),
            "round": round_id, "n": n, "a": str(a), "b": str(b),
            "d": str(d), "epsilon": str(eps), "source": source,
            "source_sha256": hashlib.sha256(canonical(source)).hexdigest(),
            "receiver": list(map(str, x))}


def register():
    cases = []
    for round_id, n in enumerate([1, 16, 256], start=1):
        for a, b in [(F(1), F(2)), (F(1), F(1) + F(1, 16 * n)), (F(2), F(3))]:
            cases.append(make_case(round_id, n, a, b))
    save(REG, {"status": "registered_before_new_winner_checks",
               "utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
               "script_sha256": digest(__file__), "cases": cases,
               "budget": "Exactly 9 new deterministic source/receiver records, 3 per n=1,16,256; q=2,4,16. No random seed because all inputs explicit rational formulas.",
               "source": "Complete two-atom positive source in one fixed half-open C=[0,d)^n, weights 1 and epsilon=n*d/(10*a), W=1+epsilon; no normalization.",
               "geometry": "D_in=a-4d/5, D_out=a+4d/5. In n>1, free-coordinate full distance 4d/5 must not dominate D_in. R+2d must remain in window; lower corner r=a+d outside original cube.",
               "winner": "Exact continuous winner from endpoints a,b and every actual arrival in (a,b]; closed capture, simultaneous equal distances; smallest-R tie rule.",
               "checks": "Exact all arrivals/full-source responses; strict unique R=a; captured labels exactly [0], full cell labels [0,1]; Bernoulli (1+4d/(5a))^n>=1+4nd/(5a)>1+epsilon; mR>tau*R^n; mR>=(R/(R+2d))^n*(1+epsilon); expanded cube captures full C labels, corner r<=R+2d; top=false.",
               "tau": "Per-row diagnostic tau=M/q for q=2,4,16; not a common-threshold Lebesgue test.",
               "scope": "Exact local partial-cell implementation guard only, not global Lebesgue integral, near-optimal source, weak counterexample, or geom qualification.",
               "stopping": "Complete all fixed 9 records once, save all outcomes; never alter old collision/age inputs or rerun old scripts."})
    print(json.dumps({"registration": str(REG), "sha256": digest(REG)}))


def process(case):
    n, a, b, d, eps = case["n"], F(case["a"]), F(case["b"]), F(case["d"]), F(case["epsilon"])
    atoms = [[F(v) for v in y] for y in case["source"]["atoms"]]
    weights = list(map(F, case["source"]["weights"]))
    x = list(map(F, case["receiver"]))
    W = sum(weights)
    assert W == 1 + eps and eps > 0
    cell_indices = [[int(coord // d) for coord in y] for y in atoms]
    assert cell_indices == [[0] * n, [0] * n]
    assert all(0 <= coord < d for y in atoms for coord in y)
    distances_per_axis = [[2 * abs(xi - yi) for xi, yi in zip(x, y)] for y in atoms]
    D = [max(axis) for axis in distances_per_axis]
    assert D == [a - 4 * d / 5, a + 4 * d / 5]
    assert d <= a / (8 * n)
    assert a - 4 * d / 5 >= 4 * d / 5
    assert D[0] < a < D[1] < b
    nodes = sorted({a, b} | {distance for distance in D if a < distance <= b})
    table = []
    for L in nodes:
        labels = [j for j, distance in enumerate(D) if distance <= L]
        mass = sum((weights[j] for j in labels), F(0))
        table.append({"L": str(L), "captured_labels": labels, "mass": str(mass),
                      "response": str(mass / L ** n)})
    M = max(F(row["response"]) for row in table)
    winners = [F(row["L"]) for row in table if F(row["response"]) == M]
    R = min(winners)
    assert winners == [a] and R == a and M == 1 / a ** n
    bernoulli_left = (1 + 4 * d / (5 * a)) ** n
    bernoulli_linear = 1 + 4 * n * d / (5 * a)
    assert bernoulli_left >= bernoulli_linear > W
    captured = [j for j, distance in enumerate(D) if distance <= R]
    assert captured == [0]
    mR = sum(weights[j] for j in captured)
    expanded = R + 2 * d
    expanded_captured = [j for j, distance in enumerate(D) if distance <= expanded]
    assert expanded_captured == [0, 1]
    assert expanded <= b and R <= b - 2 * d
    mplus = sum(weights[j] for j in expanded_captured)
    corner_r = 2 * max(abs(xi) for xi in x)
    assert corner_r == a + d and R < corner_r <= expanded
    capresidual = M * expanded ** n - mplus
    cellresidual = mR - (R / expanded) ** n * W
    assert capresidual >= 0 and cellresidual >= 0
    for row in table:
        L = F(row["L"])
        row["winner_cap_residual"] = str(M * L ** n - F(row["mass"]))
        row["strict_response_gap"] = str(M - F(row["response"]))
        assert F(row["winner_cap_residual"]) >= 0
    second = max(F(row["response"]) for row in table if F(row["L"]) != R)
    threshold_checks = []
    for q in [2, 4, 16]:
        tau = M / q
        residual = mR - tau * R ** n
        assert residual > 0
        threshold_checks.append({"q": q, "tau": str(tau), "strict_threshold_residual": str(residual), "pass": True})
    return {"case": case, "cell_indices": cell_indices,
            "D": list(map(str, D)), "full_side_distances_by_axis": [list(map(str, z)) for z in distances_per_axis],
            "arrival_table": table, "R": str(R), "M": str(M), "W": str(W),
            "actual_winner_tie_count": len(winners), "actual_winner_tie_rule": "smallest exact R",
            "relative_second_best_gap": str((M - second) / M),
            "captured_labels": captured, "expanded_captured_labels": expanded_captured,
            "partial_cell": True, "cell_full_weight": str(W), "cell_captured_mass": str(mR),
            "cell_posterior": "1", "collision": "1", "R_expanded": str(expanded),
            "expanded_mass": str(mplus), "corner": ["0"] * n, "corner_r": str(corner_r),
            "corner_outside_original": True, "top": False, "empty": False,
            "window_margin": str(b - expanded),
            "Bernoulli_left": str(bernoulli_left), "Bernoulli_linear": str(bernoulli_linear),
            "Bernoulli_nonstrict_residual": str(bernoulli_left - bernoulli_linear),
            "Bernoulli_strict_residual": str(bernoulli_linear - W),
            "expanded_winner_cap_residual": str(capresidual),
            "cell_weight_bound_residual": str(cellresidual),
            "threshold_checks": threshold_checks,
            "checks": {"positive_full_mass": True, "same_halfopen_cell": True,
                       "distances_actual_and_closed_capture": True, "free_axis_D_in_guard": True,
                       "complete_continuous_winner": True, "strict_R_equals_a": True,
                       "Bernoulli_nonstrict": True, "Bernoulli_strict": True,
                       "partial_capture": True, "expanded_full_cell_capture": True,
                       "expanded_window_eligible": True, "corner_outside_original": True,
                       "corner_expanded_bound": True, "expanded_winner_cap": True,
                       "cell_weight_bound": True}}


def run():
    start = time.monotonic()
    reg = json.loads(REG.read_text())
    assert reg["script_sha256"] == digest(__file__)
    records = []
    for case in reg["cases"]:
        assert hashlib.sha256(canonical(case["source"])).hexdigest() == case["source_sha256"]
        records.append(process(case))
    summaries = []
    for rd in [1, 2, 3]:
        rows = [row for row in records if row["case"]["round"] == rd]
        summaries.append({"round": rd, "n": rows[0]["case"]["n"], "records": len(rows),
                          "partial_cells": len(rows), "corner_outside_original": len(rows),
                          "geometry_and_winner_check_passes": sum(len(row["checks"]) for row in rows),
                          "exact_candidate_caps": sum(len(row["arrival_table"]) for row in rows),
                          "threshold_passes": sum(len(row["threshold_checks"]) for row in rows),
                          "failures": 0, "empty": 0, "top": 0,
                          "relative_gap_range": [str(min(F(row["relative_second_best_gap"]) for row in rows)),
                                                 str(max(F(row["relative_second_best_gap"]) for row in rows))]})
    resultpath = BASE / (PREFIX + "_results.json")
    save(resultpath, {"status": "complete_registered_new_partial_inputs", "records": records,
                      "summaries": summaries, "elapsed_seconds": time.monotonic() - start,
                      "registration_sha256": digest(REG), "script_sha256": digest(__file__),
                      "evidence": "All coordinates, powers, masses, responses, geometry and threshold residues are exact Fraction; no logarithm or numerical integration used.",
                      "scope": reg["scope"], "old_oracle_calls": 0, "old_input_changes": 0})
    receipt = BASE / (PREFIX + "_receipt.json")
    save(receipt, {"status": "complete", "exit_status": 0, "main_runs": 1, "old_oracle_calls": 0,
                   "old_input_changes": 0, "records": 9, "geometry_winner_check_passes": 135,
                   "exact_candidate_caps": 27, "threshold_passes": 27,
                   "registration": {"path": str(REG), "sha256": digest(REG)},
                   "script": {"path": __file__, "sha256": digest(__file__)},
                   "results": {"path": str(resultpath), "sha256": digest(resultpath)},
                   "scope": reg["scope"]})
    print(json.dumps({"summaries": summaries, "results_sha256": digest(resultpath), "receipt_sha256": digest(receipt)}))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("mode", choices=["register", "run"])
    mode = parser.parse_args().mode
    register() if mode == "register" else run()
