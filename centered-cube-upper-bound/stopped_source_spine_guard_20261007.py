"""Finite exact source-mass tilted death process and stopping guards.

Standalone Fraction tables. No maximal/oracle/old simulation imports.
"""
from fractions import Fraction as F
from pathlib import Path
import argparse
import datetime
import gzip
import hashlib
import json
import time

BASE = Path(__file__).resolve().parent
PREFIX = "stopped_source_spine_guard_20261007"
REG = BASE / (PREFIX + "_registration.json")
PS = [F(1, 2), F(1, 4), F(3, 4)]


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def canonical(obj):
    return json.dumps(obj, sort_keys=True, separators=(",", ":")).encode()


def save(path, obj):
    with Path(path).open("x") as handle:
        json.dump(obj, handle, sort_keys=True, indent=2)
        handle.write("\n")


def labels(mask, N):
    return [j for j in range(N) if mask & (1 << j)]


def subsets(mask):
    result = []
    sub = mask
    while True:
        result.append(sub)
        if sub == 0:
            return sorted(result)
        sub = (sub - 1) & mask


def harmonic(m):
    return sum((F(1, j) for j in range(1, m + 1)), F(0))


def make_case(round_id, n, raw):
    N = len(raw)
    total = sum(raw)
    weights = [F(w, total) for w in raw]
    d = F(1, 8 * n)
    positions = [[F(3 * j)] + [F(0)] * (n - 1) for j in range(N)]
    source = {"n": n, "N": N, "raw_cell_weights": raw,
              "weights": list(map(str, weights)), "W": "1",
              "atoms": [list(map(str, y)) for y in positions]}
    return {"case_id": f"r{round_id}_N{N}_n{n}", "round": round_id,
            "n": n, "N": N, "a": "1", "b": "2", "d": str(d),
            "source": source, "source_sha256": hashlib.sha256(canonical(source)).hexdigest(),
            "p": list(map(str, PS)),
            "rules": [{"id": "count_k1", "type": "count", "k": 1},
                      {"id": "count_half", "type": "count", "k": (N + 1) // 2},
                      {"id": "singleton_or_dominance", "type": "dominance", "cutoff": "3/4"}]}


def register():
    cases = [make_case(1, 1, [1, 2, 4]),
             make_case(2, 16, [1, 2, 3, 5, 11]),
             make_case(3, 256, [1, 1, 2, 3, 5, 8, 13, 55])]
    proofpaths = [BASE / "source_cell_thinning_budget_20261007.md",
                  Path("/Users/zhengzhihao/.codex/skills/math-e04-snell-envelope-optimal-stopping/SKILL.md"),
                  Path("/Users/zhengzhihao/.codex/skills/math-j01-inhomogeneous-jump-generator/SKILL.md")]
    save(REG, {
        "status": "registered_before_finite_source_spine_tables",
        "utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "script_sha256": digest(__file__), "cases": cases,
        "read_workflow_snapshots": [{"path": str(path), "sha256": digest(path)} for path in proofpaths],
        "seeds": "No random seeds or sampling; explicit fixed masses/positions and exhaustive finite state tables.",
        "spatial_scope": "Positions j*3 along first axis, half-open d=1/(8n) cells, [1,2]; only distinct complete cells verified, no cube max/receiver or dimension-order inference.",
        "Q_rates": "For all nonempty A: wA=sum retained original weights; lambda_j=1-wj/wA, total r-1, tag posterior wj/wA. Singleton total 0, no Q empty state.",
        "finite_time": "For each nonempty initial A and p=exp(-t) in {1/2,1/4,3/4}, exhaustive all B subset A including empty; P(B)=p^|B|(1-p)^(r-|B|), density=e^t*wB/wA=wB/(p*wA). Direct hidden tag chosen wj/wA survives, all other cells survive independently with p. Check both Q constructions, total probabilities and compensated mass exactly.",
        "stopping": "count<=1; count<=ceil(N/2); count<=1 OR maxweight/wA>=3/4. All rules adapted to current retained source. Threshold equality retained. Three rules evaluated at every nonempty initial state A.",
        "H": "Q mean remaining stopping time: terminal0; H(A)=1/(r-1)+sum lambda_j/(r-1)*H(A-j).",
        "G": "P expected integral of compensated source mass, E_P integral_0^T e^t*w_A_t dt: terminal0; G(A)=wA/(r-1)+sum G(A-j)/(r-1). Check G=wA*H.",
        "F": "P compensated terminal source mass E_P[e^T*w_A_T], terminal wA; F(A)=sum F(A-j)/(r-1). Check F=wA. Not raw unscaled surviving mass expectation.",
        "holding_time_justification": "P nonterminal holding Exp(r), E exp(holding)=r/(r-1), E integral_0^holding exp(t)dt=1/(r-1), uniform deleted label; r>=2 until stop. At most N-1 jumps. No unbounded optional sampling assumption.",
        "count_formula": "H(A)=0 for |A|<=k, else harmonic(|A|-1)-harmonic(k-1); full initial A gives harmonic(N-1)-harmonic(k-1).",
        "fixed_budget": {"cases": 3, "nonempty_states": 293, "death_rate_entries": 1116,
                         "finite_time_laws": 879, "finite_time_transition_rows_including_empty": 20484,
                         "stopping_state_rule_entries": 879},
        "evidence": "All probability/rate/budget values exact Fraction; full tables saved; no logarithms or floating integration.",
        "scope": "Finite source-thinning/spine/stopping implementation guard. No Phi computation, geom/history/nearmax qualification, global receiver budget, or asymptotic dimension fit.",
    })
    print(json.dumps({"registration": str(REG), "sha256": digest(REG)}))


def state_mass(mask, weights):
    return sum((weights[j] for j in labels(mask, len(weights))), F(0))


def process(case):
    N, n = case["N"], case["n"]
    weights = list(map(F, case["source"]["weights"]))
    atoms = [[F(v) for v in y] for y in case["source"]["atoms"]]
    d = F(case["d"])
    assert all(w > 0 for w in weights) and sum(weights) == 1
    cells = [tuple(int(coord // d) for coord in y) for y in atoms]
    assert len(set(cells)) == N
    assert all(all(d * index <= coord < d * (index + 1)
                   for index, coord in zip(cell, y)) for cell, y in zip(cells, atoms))
    masks = sorted(range(1, 2 ** N), key=lambda A: (A.bit_count(), A))
    masses = {0: F(0), **{A: state_mass(A, weights) for A in masks}}
    state_rows, rates = [], {}
    rate_entries, transition_entries = 0, 0
    for A in masks:
        js = labels(A, N)
        r, wA = len(js), masses[A]
        tags = {j: weights[j] / wA for j in js}
        rates[A] = {j: 1 - tags[j] for j in js}
        assert sum(tags.values()) == 1
        assert sum(rates[A].values()) == r - 1
        assert all(0 <= rate < 1 for rate in rates[A].values())
        assert wA - sum(weights[j] for j in js) == 0
        rate_entries += r
        time_laws = []
        for p in PS:
            transitions = []
            psum, qsum, weighted_mass = F(0), F(0), F(0)
            for B in subsets(A):
                Pprob = F(1)
                for j in js:
                    Pprob *= p if B & (1 << j) else 1 - p
                density = masses[B] / (p * wA)
                Qprob = Pprob * density
                direct = F(0)
                tag_components = []
                for tag in js:
                    conditional = F(1) if B & (1 << tag) else F(0)
                    for j in js:
                        if j != tag:
                            conditional *= p if B & (1 << j) else 1 - p
                    contribution = tags[tag] * conditional
                    direct += contribution
                    tag_components.append({"tag": tag, "conditional_survival_probability": str(conditional),
                                           "tag_weighted_probability": str(contribution)})
                assert Qprob == direct
                assert B != 0 or Qprob == 0
                psum += Pprob
                qsum += Qprob
                weighted_mass += Pprob * masses[B] / p
                transitions.append({"B_mask": B, "B_labels": labels(B, N), "wB": str(masses[B]),
                                    "P_probability": str(Pprob), "mass_tilt_density": str(density),
                                    "Q_probability": str(Qprob), "direct_tag_mixture": str(direct),
                                    "tag_components": tag_components})
            assert psum == qsum == 1 and weighted_mass == wA
            transition_entries += len(transitions)
            time_laws.append({"p": str(p), "P_sum": "1", "Q_sum": "1",
                              "E_P_compensated_mass": str(weighted_mass), "transitions": transitions})
        state_rows.append({"A_mask": A, "A_labels": js, "count": r, "wA": str(wA),
                           "tag_posterior": [{"j": j, "posterior": str(tags[j])} for j in js],
                           "Q_death_rates": [{"j": j, "lambda": str(rates[A][j]), "child_mask": A ^ (1 << j)} for j in js],
                           "Q_total_rate": r - 1, "Q_constant_generator_residual": "0",
                           "P_mass_compensation_residual": "0", "finite_time_laws": time_laws})
    stopping_tables = []
    count_formula_checks = 0
    for rule in case["rules"]:
        H, G, FF = {}, {}, {}
        table = []
        for A in masks:
            js, r, wA = labels(A, N), A.bit_count(), masses[A]
            max_share = max(weights[j] for j in js) / wA
            terminal = r <= rule["k"] if rule["type"] == "count" else r <= 1 or max_share >= F(rule["cutoff"])
            contributions = []
            if terminal:
                H[A], G[A], FF[A] = F(0), F(0), wA
            else:
                assert r >= 2
                H[A], G[A], FF[A] = F(1, r - 1), wA / (r - 1), F(0)
                for j in js:
                    child = A ^ (1 << j)
                    assert child > 0
                    hc = rates[A][j] / (r - 1) * H[child]
                    gc = G[child] / (r - 1)
                    fc = FF[child] / (r - 1)
                    H[A] += hc
                    G[A] += gc
                    FF[A] += fc
                    contributions.append({"deleted_j": j, "child_mask": child,
                                          "Q_embedded_probability": str(rates[A][j] / (r - 1)),
                                          "H_contribution": str(hc), "G_contribution": str(gc),
                                          "F_contribution": str(fc)})
            assert G[A] == wA * H[A] and FF[A] == wA
            count_expected = None
            if rule["type"] == "count":
                k = rule["k"]
                count_expected = F(0) if r <= k else harmonic(r - 1) - harmonic(k - 1)
                assert H[A] == count_expected
                count_formula_checks += 1
            table.append({"A_mask": A, "A_labels": js, "count": r, "wA": str(wA),
                          "terminal": terminal, "maxweight_share": str(max_share),
                          "dominance_at_equality": rule["type"] == "dominance" and max_share == F(3, 4),
                          "H_Q_expected_time": str(H[A]),
                          "G_P_integrated_compensated_mass": str(G[A]),
                          "F_P_compensated_terminal_mass": str(FF[A]),
                          "G_minus_wA_H": "0", "F_minus_wA": "0",
                          "count_harmonic_expected": str(count_expected) if count_expected is not None else None,
                          "Q_mean_first_hold": "0" if terminal else str(F(1, r - 1)),
                          "P_expected_growth_before_jump": "1" if terminal else str(F(r, r - 1)),
                          "P_first_hold_integrated_mass": "0" if terminal else str(wA / (r - 1)),
                          "contributions": contributions})
        stopping_tables.append({"rule": rule, "states": table,
                                "full_initial_H": str(H[2 ** N - 1]),
                                "full_initial_G": str(G[2 ** N - 1]),
                                "full_initial_F": str(FF[2 ** N - 1]),
                                "terminal_states": sum(row["terminal"] for row in table),
                                "dominance_equality_states": sum(row["dominance_at_equality"] for row in table)})
    fullH = {table["rule"]["id"]: table["full_initial_H"] for table in stopping_tables}
    summary = {"case_id": case["case_id"], "round": case["round"], "N": N, "n": n,
               "nonempty_states": len(masks), "distinct_source_cells": N,
               "death_rate_entries": rate_entries, "rate_total_checks": len(masks), "tag_sum_checks": len(masks),
               "finite_time_laws": len(masks) * len(PS), "finite_time_transition_rows": transition_entries,
               "tilt_tag_identity_checks": transition_entries,
               "stopping_state_rule_entries": len(masks) * len(case["rules"]),
               "G_equals_wA_H_checks": len(masks) * len(case["rules"]),
               "F_equals_wA_checks": len(masks) * len(case["rules"]),
               "count_harmonic_checks": count_formula_checks, "full_initial_H": fullH,
               "dominance_equality_states": stopping_tables[-1]["dominance_equality_states"], "failures": 0}
    return {"case": case, "source_cell_indices": [list(cell) for cell in cells],
            "states": state_rows, "stopping_tables": stopping_tables, "summary": summary}


def run():
    start = time.monotonic()
    reg = json.loads(REG.read_text())
    assert reg["script_sha256"] == digest(__file__)
    manifests, summaries = [], []
    for case in reg["cases"]:
        assert hashlib.sha256(canonical(case["source"])).hexdigest() == case["source_sha256"]
        result = process(case)
        path = BASE / (PREFIX + "_" + case["case_id"] + ".json.gz")
        with path.open("xb") as raw:
            with gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as zipped:
                zipped.write(canonical(result))
        manifests.append({"path": str(path), "sha256": digest(path)})
        summaries.append(result["summary"])
        print(json.dumps({"done": case["case_id"], "states": result["summary"]["nonempty_states"],
                          "transitions": result["summary"]["finite_time_transition_rows"]}), flush=True)
    resultpath = BASE / (PREFIX + "_results.json")
    save(resultpath, {"status": "complete_registered_exact_source_spine_guard", "summaries": summaries,
                      "profiles": manifests, "elapsed_seconds": time.monotonic() - start,
                      "registration_sha256": digest(REG), "script_sha256": digest(__file__),
                      "scope": reg["scope"], "old_oracle_calls": 0, "Phi_evaluations": 0,
                      "evidence": "Fraction complete finite transition/first-jump tables; n only confirms disjoint fixed source cells."})
    receiptpath = BASE / (PREFIX + "_receipt.json")
    save(receiptpath, {"status": "complete", "exit_status": 0, "main_runs": 1,
                       "registration": {"path": str(REG), "sha256": digest(REG)},
                       "script": {"path": __file__, "sha256": digest(__file__)},
                       "results": {"path": str(resultpath), "sha256": digest(resultpath)},
                       "profiles": manifests, "old_oracle_calls": 0, "Phi_evaluations": 0,
                       "scope": reg["scope"]})
    print(json.dumps({"results_sha256": digest(resultpath), "receipt_sha256": digest(receiptpath)}))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("mode", choices=["register", "run"])
    mode = parser.parse_args().mode
    register() if mode == "register" else run()
