"""New deterministic full-source continuous-window posterior-age guards.

No old oracle/MC imports. Register first, then --run; all outputs exclusive-create.
"""
from __future__ import annotations
import argparse
from collections import defaultdict
from decimal import Decimal, localcontext
from fractions import Fraction as F
from itertools import product
from math import isqrt
from pathlib import Path
import datetime
import gzip
import hashlib
import json
import time

BASE = Path(__file__).resolve().parent
PREFIX = "posterior_radial_age_probe_20261007"
REG = BASE / (PREFIX + "_registration.json")
UNIT = 2 ** 128


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def canonical(value):
    return json.dumps(value, sort_keys=True, separators=(",", ":")).encode()


def save(path, obj):
    with Path(path).open("x") as handle:
        json.dump(obj, handle, indent=2, sort_keys=True)
        handle.write("\n")


def make_case(round_id, model, n, N, a, b, seed, receivers):
    import numpy as np
    rng = np.random.default_rng(seed)
    if model == "small_random":
        atoms = [[F(int(x), 4) for x in row]
                 for row in rng.integers(-8, 9, size=(N, n))]
    elif model == "correlated_ray":
        direction = [F((-1) ** i * (i % 3 + 1), 3) for i in range(n)]
        atoms = [[F(j - N // 2, 4) * z for z in direction] for j in range(N)]
    elif model == "nested_sign_layers":
        atoms = [[F(2 ** (j % 5), 16) * (-1) ** (j // 5)
                  * (1 if i < n // 2 else (-1) ** (j % 3))
                  for i in range(n)] for j in range(N)]
    elif model == "separated_clusters":
        atoms = [[F(3 * (j % 3 - 1)) + F(int(x), 16)
                  for x in rng.integers(-1, 2, size=n)] for j in range(N)]
    elif model == "finite_lattice":
        catalogue = list(product([-1, 0, 1], repeat=n))
        indices = rng.choice(len(catalogue), size=N, replace=False)
        atoms = [[F(x, 2) for x in catalogue[int(j)]] for j in indices]
    elif model == "large_block_layers":
        centers = rng.integers(-4, 5, size=N)
        signs = rng.choice([-1, 1], size=(N, n))
        atoms = [[F(int(centers[j]) * (1 if i < n // 2 else -1), 4)
                  + F((j % 4 + 1) * int(signs[j, i]), 64)
                  for i in range(n)] for j in range(N)]
    else:
        raise ValueError(model)
    atoms[-1] = atoms[0][:]
    raw_weights = [int(x) for x in rng.integers(1, 10, size=N)]
    weights = [F(x, sum(raw_weights)) for x in raw_weights]
    ar, br = F(a), F(b)
    rr = np.random.default_rng(seed + 500000)
    xs = []
    for j in range(receivers):
        if j == 0:
            x, kind = atoms[0][:], "source_center"
        elif j == 1:
            dx = [ar * F(int(v), 16) for v in rr.integers(-7, 8, size=n)]
            dx[0] = ar / 2
            x = [z + d for z, d in zip(atoms[0], dx)]
            kind = "arrival_at_a_face"
        elif j == 2:
            r = (ar + br) / 2
            dx = [r * F(int(v), 16) for v in rr.integers(-7, 8, size=n)]
            dx[-1] = -r / 2
            x = [z + d for z, d in zip(atoms[1], dx)]
            kind = "middle_radius_face"
        elif j == 3:
            x = [z + (-1) ** i * br / 2 for i, z in enumerate(atoms[-1])]
            kind = "closed_b_corner"
        elif j == 4:
            x, kind = [F(0)] * n, "origin"
        else:
            x = [max(z[i] for z in atoms) + 10 * br for i in range(n)]
            kind = "outside_all_response"
        xs.append({"x": list(map(str, x)), "kind": kind})
    source = {"n": n, "atoms": [list(map(str, y)) for y in atoms],
              "weights": list(map(str, weights))}
    return {"case_id": f"r{round_id}_{model}_n{n}_N{N}_a{ar}_b{br}".replace("/", "d"),
            "round": round_id, "model": model, "n": n, "N": N,
            "a": str(ar), "b": str(br), "seed": seed,
            "receiver_seed": seed + 500000,
            "source": source, "source_sha256": hashlib.sha256(canonical(source)).hexdigest(),
            "receivers": xs}


def register():
    plans = []
    for n, N, seed in [(1, 4, 852001), (2, 7, 852002), (3, 10, 852003)]:
        plans.append(make_case(1, "small_random", n, N, "1", "2", seed, 6))
    for model, n, N, seed in [
        ("correlated_ray", 4, 12, 8530412),
        ("nested_sign_layers", 8, 16, 8530816),
        ("separated_clusters", 8, 24, 8530824),
        ("finite_lattice", 4, 32, 8530432),
    ]:
        plans.append(make_case(2, model, n, N, "1", "2", seed, 6))
    for n, N, a, b, seed in [
        (16, 8, "1", "2", 8541608), (16, 32, "1", "2", 8541632),
        (64, 8, "1", "2", 8546408), (64, 32, "1", "2", 8546432),
        (64, 32, "1/2", "3/2", 8546432),
        (256, 8, "1", "2", 85425608), (256, 32, "1", "2", 85425632),
        (256, 32, "1", "3", 85425632),
    ]:
        plans.append(make_case(3, "large_block_layers", n, N, a, b, seed, 4))
    payload = {
        "status": "registered_before_response_oracle_execution",
        "utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "script_sha256": digest(__file__), "plans": plans,
        "fixed_budget": {"cases": 15, "receivers": 74, "coarse_counterexamples": 1},
        "source": "Complete positive labelled measure of mass 1; coincident labels retained, never componentwise maxima.",
        "winner": "Exact Fraction: endpoints a,b and all full-side l-infinity arrivals in (a,b]; cumulative closed-cube ties updated simultaneously; smallest R among exact maximizers.",
        "candidate_L": "All endpoints/arrivals, adjacent midpoints and (a+b)/2, deduplicated; zero captured mass marked undefined.",
        "primary_check": "For every positive-mass L and s=a,L plus each real arrival in [a,L], m_s <= M*s**n, exactly; equivalent posterior tail relation recorded. Between arrivals mass constant and M*s**n nondecreasing.",
        "age_diagnostics": "g and E A at Decimal precision 160, not logarithm interval certification.",
        "sigma_half": "128-bit outward rational square-root bounds of (L/max(a,D))**n and M*L**n/mL; compare upper posterior moment to lower 2*sqrt(c). Crossing intervals means unresolved, not failure of theorem.",
        "nearwinner": "Exact Fraction relative second-best gap and exact tie count, no fuzzy or near-threshold replacement of definition.",
        "save": "Full input sources, receiver coordinates, D labels, arrival cumulative masses, actual M/R and every L posterior; sufficient for later separate periodic-coordinate postprocessing.",
        "counterexample": {"n": 1, "atoms": [["0"]], "weights": ["1"], "a": "1", "b": "2", "x": ["3/4"], "coarse_scales": ["1", "2"], "omitted_arrival": "3/2", "scope": "Finite-grid M cannot cap omitted continuous arrival; not a weak-maximal counterexample."},
        "stopping": "Fixed complete inputs: finish all 15 cases and all 74 receivers; no parameter retuning, no old oracle or MC execution.",
        "scope": "Finite deterministic receiver guards only, not Lebesgue probability/integral, geom qualification, actual FIRST or general asymptotic evidence.",
    }
    save(REG, payload)
    print(json.dumps({"registration": str(REG), "sha256": digest(REG), "cases": len(plans)}))


def sqrt_interval(value):
    assert value >= 0
    z = isqrt((value.numerator * UNIT * UNIT) // value.denominator)
    lo = F(z, UNIT)
    hi = lo if lo * lo == value else F(z + 1, UNIT)
    assert lo * lo <= value <= hi * hi
    return lo, hi


def dec(value):
    return Decimal(value.numerator) / Decimal(value.denominator)


def probe_receiver(case, receiver):
    n, a, b = case["n"], F(case["a"]), F(case["b"])
    atoms = [[F(v) for v in y] for y in case["source"]["atoms"]]
    weights = list(map(F, case["source"]["weights"]))
    x = list(map(F, receiver["x"]))
    assert sum(weights) == 1 and all(w > 0 for w in weights)
    D = [2 * max(abs(xi - yi) for xi, yi in zip(x, y)) for y in atoms]
    groups = defaultdict(list)
    for j, distance in enumerate(D):
        groups[distance].append(j)
    nodes = sorted({a, b} | {s for s in D if a < s <= b})
    table = []
    captured = [j for j, distance in enumerate(D) if distance <= a]
    mass = sum((weights[j] for j in captured), F(0))
    for s in nodes:
        new = groups[s] if s > a else captured
        if s > a:
            mass += sum((weights[j] for j in new), F(0))
        table.append({"s": str(s), "mass": str(mass), "response": str(mass / s ** n),
                      "arriving_labels": list(new)})
    M = max(F(row["response"]) for row in table)
    if M == 0:
        return {"kind": receiver["kind"], "x": receiver["x"], "D": list(map(str, D)),
                "M": "0", "R": None, "empty_response": True, "arrival_table": table,
                "candidate_L": [], "primary_checks": 0}
    winners = [row for row in table if F(row["response"]) == M]
    R = F(winners[0]["s"])
    other = [F(row["response"]) for row in table if F(row["s"]) != R]
    second = max(other, default=F(0))
    gap = (M - second) / M
    for row in table:
        s = F(row["s"])
        row["global_cap_residual"] = str(M * s ** n - F(row["mass"]))
        assert F(row["global_cap_residual"]) >= 0
    Ls = sorted(set(nodes) | {(u + v) / 2 for u, v in zip(nodes, nodes[1:])} | {(a + b) / 2})
    candidates, checks = [], 0
    logcache = {}
    def logq(value):
        if value not in logcache:
            logcache[value] = dec(value).ln()
        return logcache[value]
    with localcontext() as context:
        context.prec = 160
        for L in Ls:
            labels = [j for j, distance in enumerate(D) if distance <= L]
            mL = sum((weights[j] for j in labels), F(0))
            if mL == 0:
                candidates.append({"L": str(L), "mass": "0", "undefined_posterior": True})
                continue
            posterior = [weights[j] / mL if j in labels else F(0) for j in range(len(weights))]
            assert sum(posterior) == 1
            c = M * L ** n / mL
            assert c >= 1
            checked = sorted({a, L} | {s for s in D if a <= s <= L})
            for s in checked:
                ms = sum((weights[j] for j, distance in enumerate(D) if distance <= s), F(0))
                assert ms <= M * s ** n
                assert ms / mL <= c * (s / L) ** n
                checks += 1
            g = logq(c)
            EA = sum((dec(posterior[j]) * n * (logq(L) - logq(max(a, D[j])))
                      for j in labels), Decimal(0))
            mlower, mupper = F(0), F(0)
            for j in labels:
                slo, shi = sqrt_interval((L / max(a, D[j])) ** n)
                mlower += posterior[j] * slo
                mupper += posterior[j] * shi
            clo, chi = sqrt_interval(c)
            rhslo, rhshi = 2 * clo, 2 * chi
            candidates.append({
                "L": str(L), "mass": str(mL), "captured_labels": labels,
                "posterior": list(map(str, posterior)), "c_exp_g": str(c),
                "age_reference_max_a_D": list(map(str, [max(a, d) for d in D])),
                "g_decimal160": str(g), "EA_decimal160": str(EA),
                "EA_bound_residual_decimal160": str(g + 1 - EA),
                "EA_diagnostic_pass": EA <= g + 1,
                "tail_checked_s": list(map(str, checked)),
                "tail_exact_checks": len(checked),
                "L_cap_residual": str(M * L ** n - mL),
                "sigma_half_moment_interval": [str(mlower), str(mupper)],
                "sigma_half_bound_interval": [str(rhslo), str(rhshi)],
                "sigma_half_enclosure_pass": mupper <= rhslo,
                "L_at_or_above_winner": L >= R,
                "winner_L": L == R,
            })
        logarithmic_gap = str(logq(M / second)) if second > 0 else None
    return {"kind": receiver["kind"], "x": receiver["x"], "D": list(map(str, D)),
            "M": str(M), "R": str(R), "empty_response": False,
            "exact_winner_tie_count": len(winners), "tied_R": [row["s"] for row in winners],
            "second_best_response": str(second), "relative_nearwinner_gap": str(gap),
            "log_nearwinner_gap_decimal160": logarithmic_gap,
            "arrival_table": table, "candidate_L": candidates, "primary_checks": checks}


def coarse_counterexample():
    # Exact separate registered example; no random search or parameter choice.
    D, a, b = F(3, 2), F(1), F(2)
    true_M, true_R, coarse_M = F(2, 3), D, F(1, 2)
    mass, L = F(1), b
    coarse_c = coarse_M * L / mass
    e_minus_t = D / L
    assert mass > coarse_M * D
    assert mass == true_M * D
    assert coarse_c * e_minus_t == F(3, 4)
    return {"n": 1, "a": str(a), "b": str(b), "source": {"atoms": [["0"]], "weights": ["1"]},
            "x": ["3/4"], "D": str(D), "true_R": str(true_R), "true_M": str(true_M),
            "coarse_scales": ["1", "2"], "coarse_M": str(coarse_M),
            "cap_violation_positive": str(mass - coarse_M * D),
            "L": str(L), "t": "log(4/3)", "posterior_tail": "1",
            "coarse_exp_g_minus_t": "3/4", "true_cap_residual": "0",
            "scope": "The coarse finite family is valid for its own max, but cannot replace the continuous max in an omitted-arrival cap."}


def run():
    start = time.monotonic()
    reg = json.loads(REG.read_text())
    assert reg["script_sha256"] == digest(__file__), "Script changed after registration"
    summaries, profile_files = [], []
    for case in reg["plans"]:
        assert hashlib.sha256(canonical(case["source"])).hexdigest() == case["source_sha256"]
        records = [probe_receiver(case, receiver) for receiver in case["receivers"]]
        filename = BASE / (PREFIX + "_" + case["case_id"] + ".json.gz")
        with filename.open("xb") as raw:
            with gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as zipped:
                zipped.write(canonical({"case": case, "records": records}))
        profile_files.append({"path": str(filename), "sha256": digest(filename)})
        positive = [record for record in records if not record["empty_response"]]
        Cs = [c for record in positive for c in record["candidate_L"] if c["mass"] != "0"]
        EAs = [Decimal(c["EA_decimal160"]) for c in Cs]
        winner_EAs = [Decimal(c["EA_decimal160"]) for c in Cs if c["winner_L"]]
        summary = {
            "case_id": case["case_id"], "round": case["round"], "n": case["n"], "N": case["N"],
            "a": case["a"], "b": case["b"], "model": case["model"],
            "seed": case["seed"], "source_sha256": case["source_sha256"],
            "source_weight_min": str(min(map(F, case["source"]["weights"]))),
            "source_weight_max": str(max(map(F, case["source"]["weights"]))),
            "distinct_physical_atoms": len(set(map(tuple, case["source"]["atoms"]))),
            "receiver_count": len(records), "positive_response_receivers": len(positive),
            "empty_response_receivers": len(records) - len(positive),
            "positive_L_count": len(Cs),
            "primary_exact_checks": sum(record["primary_checks"] for record in records),
            "primary_exact_failures": 0,
            "sigma_half_enclosure_passes": sum(c["sigma_half_enclosure_pass"] for c in Cs),
            "sigma_half_unresolved": sum(not c["sigma_half_enclosure_pass"] for c in Cs),
            "EA_diagnostic_passes": sum(c["EA_diagnostic_pass"] for c in Cs),
            "EA_diagnostic_failures": sum(not c["EA_diagnostic_pass"] for c in Cs),
            "EA_range_decimal160": [str(min(EAs)), str(max(EAs))] if EAs else None,
            "winner_EA_range_decimal160": [str(min(winner_EAs)), str(max(winner_EAs))] if winner_EAs else None,
            "winner_g_zero_exact_count": sum(c["c_exp_g"] == "1" for c in Cs if c["winner_L"]),
            "exact_multiscale_tie_receivers": sum(r["exact_winner_tie_count"] > 1 for r in positive),
            "relative_nearwinner_gap_range": [str(min(F(r["relative_nearwinner_gap"]) for r in positive)), str(max(F(r["relative_nearwinner_gap"]) for r in positive))] if positive else None,
        }
        summaries.append(summary)
        print(json.dumps({"done": case["case_id"], "checks": summary["primary_exact_checks"], "L": len(Cs)}), flush=True)
    result = {
        "status": "complete_fixed_registered_batch", "elapsed_seconds": time.monotonic() - start,
        "registration_sha256": digest(REG), "script_sha256": digest(__file__),
        "summaries": summaries, "profiles": profile_files,
        "coarse_counterexample": coarse_counterexample(),
        "evidence": "Exact Fraction continuous-winner/tail algebra; outward rational sigma-half enclosures for finite saved cases; EA/g logs exploratory Decimal160 only.",
        "scope": reg["scope"], "phase_postprocessing": "Not executed; complete source, posterior, true R, all L saved for separately registered follow-up.",
    }
    output = BASE / (PREFIX + "_results.json")
    save(output, result)
    receipt = BASE / (PREFIX + "_receipt.json")
    save(receipt, {"status": "complete", "exit_status": 0, "reruns": 0,
                   "registration": {"path": str(REG), "sha256": digest(REG)},
                   "script": {"path": __file__, "sha256": digest(__file__)},
                   "results": {"path": str(output), "sha256": digest(output)},
                   "profiles": profile_files, "elapsed_seconds": result["elapsed_seconds"],
                   "scope": reg["scope"]})
    print(json.dumps({"results": str(output), "results_sha256": digest(output), "receipt_sha256": digest(receipt)}))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("mode", choices=["register", "run"])
    args = parser.parse_args()
    register() if args.mode == "register" else run()
