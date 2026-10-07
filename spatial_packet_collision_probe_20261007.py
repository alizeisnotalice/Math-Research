"""Registered postprocessing of saved continuous winners; no maximization.

Fraction geometry, complete half-open cell source allocation, exact Bernoulli
subset-count DP, and outward rational logarithm bounds. Only new prefix writes.
"""
from collections import defaultdict
from fractions import Fraction as F
from math import gcd, lcm
from pathlib import Path
import argparse
import datetime
import gzip
import hashlib
import json
import time

BASE = Path(__file__).resolve().parent
PREFIX = "spatial_packet_collision_probe_20261007"
REG = BASE / (PREFIX + "_registration.json")
AGE = "posterior_radial_age_probe_20261007"
LOG_BITS = 192
LOG_TERMS = 100
UNIT = 2 ** LOG_BITS


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def canonical(obj):
    return json.dumps(obj, sort_keys=True, separators=(",", ":")).encode()


def save(path, obj):
    with Path(path).open("x") as handle:
        json.dump(obj, handle, indent=2, sort_keys=True)
        handle.write("\n")


def register():
    oldresult = BASE / (AGE + "_results.json")
    oldreg = BASE / (AGE + "_registration.json")
    payload = json.loads(oldresult.read_text())
    manifest = payload["profiles"]
    assert len(manifest) == 15
    for item in manifest:
        assert digest(item["path"]) == item["sha256"]
    save(REG, {
        "status": "registered_before_new_postprocessing", "utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "script_sha256": digest(__file__),
        "input_results": {"path": str(oldresult), "sha256": digest(oldresult)},
        "input_registration": {"path": str(oldreg), "sha256": digest(oldreg)},
        "input_profiles": manifest,
        "budget": {"cases": 15, "receivers": 74, "q_values": [2, 4, 16]},
        "seeds": "Original frozen source/receiver seeds retained; no new random samples; independent cell Bernoulli law integrated exactly by DP.",
        "d": "min(a/(8*n),(b-a)/4)",
        "cells": "Fixed half-open axis grid floor(y_i/d), lower corner d*cell_index; complete cell weight includes all unselected labels and fractional coincident allocation.",
        "geometry": "Saved actual R/M only: eligible R<=b-2*d; closed saved D determines mR,m_(R+2d). All labels in each captured cell captured at R+2d; corner r<=R+2d; mR>tau*R^n; eligible mR>= (R/(R+2d))^n*wC and m_(R+2d)<=M*(R+2d)^n.",
        "tau": "Per-receiver diagnostics tau=M/q, q=2,4,16. Not a common-threshold Lebesgue experiment.",
        "posterior": "piC=mu(C intersect Q_R)/mR, collision chi=sum piC^2; source allocation never renormalized cellwise.",
        "bernoulli": "One independent fair etaC per complete cell, nu=2*sum etaC*muC; captured cells determine V=2*sum piC*etaC-1. Exact subset integer-mass DP; unselected cells integrate out.",
        "local_contract": "E[log(q)-log_+(q*(1+V))] <= (2+4*log(q))*chi.",
        "log_interval": {"range_reduction": "x=2^k*t, 1<=t<2; log t=2*atanh((t-1)/(t+1))", "series_terms": LOG_TERMS, "tail_upper": "2*z^(2K+1)/((2K+1)*(1-z^2))", "outward_dyadic_bits": LOG_BITS, "check": "Expectation upper <= RHS lower, exact Fraction; unresolved retained if intervals overlap."},
        "stopping": "One complete postprocessing pass through all frozen cases; all empty/top/failure branches saved; no tuning and no old oracle execution.",
        "scope": "Finite local geometry/Bernoulli guards only, not numerical validation of global 6W radial integral budget or Lebesgue coverage; no asymptotic claim.",
    })
    print(json.dumps({"registration": str(REG), "sha256": digest(REG)}))


def outward(lo, hi):
    assert lo <= hi
    low = F((lo.numerator * UNIT) // lo.denominator, UNIT)
    upper_integer = -((-hi.numerator * UNIT) // hi.denominator)
    high = F(upper_integer, UNIT)
    assert low <= lo <= hi <= high
    return low, high


def atanh_log_interval(t):
    assert 1 <= t <= 2
    z = (t - 1) / (t + 1)
    power = z
    total = F(0)
    for j in range(LOG_TERMS):
        total += 2 * power / (2 * j + 1)
        power *= z * z
    remainder = 2 * power / ((2 * LOG_TERMS + 1) * (1 - z * z))
    return outward(total, total + remainder)


LOG_TWO = atanh_log_interval(F(2))
LOG_CACHE = {F(1): (F(0), F(0)), F(2): LOG_TWO}


def log_interval(value):
    assert value >= 1
    if value in LOG_CACHE:
        return LOG_CACHE[value]
    k = value.numerator.bit_length() - value.denominator.bit_length()
    k = max(k, 0)
    while value < 2 ** k:
        k -= 1
    while value >= 2 ** (k + 1):
        k += 1
    t = value / 2 ** k
    lo, hi = atanh_log_interval(t)
    answer = outward(k * LOG_TWO[0] + lo, k * LOG_TWO[1] + hi)
    LOG_CACHE[value] = answer
    return answer


def log_plus_interval(value):
    return (F(0), F(0)) if value <= 1 else log_interval(value)


def posterior_subset_dp(pi):
    assert sum(pi) == 1 and all(p > 0 for p in pi)
    denominator = lcm(*(p.denominator for p in pi))
    units = [int(p * denominator) for p in pi]
    assert sum(units) == denominator
    dp = {0: 1}
    for weight in units:
        updated = defaultdict(int)
        for mass, multiplicity in dp.items():
            updated[mass] += multiplicity
            updated[mass + weight] += multiplicity
        dp = dict(updated)
    outcomes = 2 ** len(pi)
    assert sum(dp.values()) == outcomes
    assert sum(F(score * count, denominator * outcomes) for score, count in dp.items()) == F(1, 2)
    assert all(dp.get(denominator - score) == count for score, count in dp.items())
    return denominator, units, outcomes, sorted(dp.items())


def process_row(case, row):
    n, a, b = case["n"], F(case["a"]), F(case["b"])
    d = min(a / (8 * n), (b - a) / 4)
    source = case["source"]
    atoms = [[F(v) for v in atom] for atom in source["atoms"]]
    weights = list(map(F, source["weights"]))
    x, D = list(map(F, row["x"])), list(map(F, row["D"]))
    assert len(D) == len(atoms) == len(weights) and sum(weights) == 1
    cells = defaultdict(list)
    for j, y in enumerate(atoms):
        cells[tuple(int(v // d) for v in y)].append(j)
    fullcells = []
    for index, labels in sorted(cells.items()):
        assert all(all(d * h <= coord < d * (h + 1) for h, coord in zip(index, atoms[j])) for j in labels)
        fullcells.append({"index": list(index), "labels": labels,
                          "weight": str(sum((weights[j] for j in labels), F(0))),
                          "corner": list(map(str, [d * h for h in index]))})
    basic = {"x": row["x"], "kind": row["kind"], "D": row["D"], "d": str(d),
             "M": row["M"], "R": row["R"], "full_cells": fullcells}
    if row["empty_response"]:
        return dict(basic, branch="empty", tau=None, posterior=None, local_checks=[])
    M, R = F(row["M"]), F(row["R"])
    mR = sum((weights[j] for j, distance in enumerate(D) if distance <= R), F(0))
    assert mR > 0 and M * R ** n == mR
    expanded = R + 2 * d
    mplus = sum((weights[j] for j, distance in enumerate(D) if distance <= expanded), F(0))
    eligible = expanded <= b
    branch = "eligible" if eligible else "top"
    pi, capturedcells = [], []
    for cell in fullcells:
        labels = cell["labels"]
        captured = [j for j in labels if D[j] <= R]
        if not captured:
            continue
        wC = F(cell["weight"])
        captured_mass = sum((weights[j] for j in captured), F(0))
        p = captured_mass / mR
        pi.append(p)
        corner = list(map(F, cell["corner"]))
        corner_r = 2 * max(abs(xi - yi) for xi, yi in zip(x, corner))
        geometry = {
            "index": cell["index"], "labels": labels, "captured_labels": captured,
            "full_weight": str(wC), "captured_mass": str(captured_mass), "pi": str(p),
            "corner": cell["corner"], "corner_r": str(corner_r),
            "all_cell_labels_in_expanded_cube": all(D[j] <= expanded for j in labels),
            "corner_bound": corner_r <= expanded,
            "cell_weight_bound_residual": str(mR - (R / expanded) ** n * wC),
            "cell_weight_bound_applicable": eligible,
        }
        assert geometry["all_cell_labels_in_expanded_cube"] and geometry["corner_bound"]
        if eligible:
            assert F(geometry["cell_weight_bound_residual"]) >= 0
        capturedcells.append(geometry)
    assert sum(pi) == 1
    collision = sum((p * p for p in pi), F(0))
    assert 0 < collision <= 1
    capresidual = M * expanded ** n - mplus
    if eligible:
        assert capresidual >= 0
    denominator, units, outcomes, dp = posterior_subset_dp(pi)
    localchecks = []
    for q in [2, 4, 16]:
        zlo, zhi = log_interval(F(q))
        expected_log_lo = expected_log_hi = F(0)
        for score, count in dp:
            value = F(2 * q * score, denominator)
            lo, hi = log_plus_interval(value)
            probability = F(count, outcomes)
            expected_log_lo += probability * lo
            expected_log_hi += probability * hi
        loss_lo, loss_hi = zlo - expected_log_hi, zhi - expected_log_lo
        rhs_lo, rhs_hi = (2 + 4 * zlo) * collision, (2 + 4 * zhi) * collision
        tau = M / q
        threshold_residual = mR - tau * R ** n
        assert threshold_residual > 0
        localchecks.append({"q": q, "tau": str(tau), "strict_threshold_residual": str(threshold_residual),
                            "log_q_interval": [str(zlo), str(zhi)],
                            "expected_log_plus_interval": [str(expected_log_lo), str(expected_log_hi)],
                            "expected_loss_interval": [str(loss_lo), str(loss_hi)],
                            "bound_interval": [str(rhs_lo), str(rhs_hi)],
                            "gap_interval": [str(rhs_lo - loss_hi), str(rhs_hi - loss_lo)],
                            "local_inequality_interval_pass": loss_hi <= rhs_lo})
    return dict(basic, branch=branch, mR=str(mR), expanded_R=str(expanded), expanded_mass=str(mplus),
                expanded_cap_residual=str(capresidual), expanded_cap_applicable=eligible,
                captured_cells=capturedcells, posterior_cell=list(map(str, pi)), collision=str(collision),
                subset_dp={"denominator": denominator, "integer_cell_mass": units,
                           "outcomes": outcomes, "states": [[score, count] for score, count in dp],
                           "mean_subset_mass": "1/2", "V_mean": "0"}, local_checks=localchecks)


def run():
    start = time.monotonic()
    reg = json.loads(REG.read_text())
    assert reg["script_sha256"] == digest(__file__)
    summaries, manifest = [], []
    for item in reg["input_profiles"]:
        assert digest(item["path"]) == item["sha256"]
        old = json.loads(gzip.decompress(Path(item["path"]).read_bytes()))
        case = old["case"]
        assert hashlib.sha256(canonical(case["source"])).hexdigest() == case["source_sha256"]
        records = [process_row(case, row) for row in old["records"]]
        output = BASE / (PREFIX + "_" + case["case_id"] + ".json.gz")
        with output.open("xb") as raw:
            with gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as zipped:
                zipped.write(canonical({"case": case, "input_profile": item, "records": records}))
        manifest.append({"path": str(output), "sha256": digest(output)})
        positive = [row for row in records if row["branch"] != "empty"]
        cells = [cell for row in positive for cell in row["captured_cells"]]
        eligiblecells = [cell for cell in cells if cell["cell_weight_bound_applicable"]]
        local = [check for row in positive for check in row["local_checks"]]
        summaries.append({"case_id": case["case_id"], "round": case["round"], "n": case["n"], "N": case["N"],
                          "a": case["a"], "b": case["b"], "seed": case["seed"],
                          "source_sha256": case["source_sha256"], "receiver_count": len(records),
                          "empty": sum(row["branch"] == "empty" for row in records),
                          "top": sum(row["branch"] == "top" for row in records),
                          "eligible": sum(row["branch"] == "eligible" for row in records),
                          "captured_cell_incidents": len(cells), "eligible_cell_weight_checks": len(eligiblecells),
                          "all_cell_capture_checks": len(cells), "corner_checks": len(cells),
                          "local_q_checks": len(local), "local_interval_passes": sum(check["local_inequality_interval_pass"] for check in local),
                          "local_unresolved": sum(not check["local_inequality_interval_pass"] for check in local),
                          "collision_range": [str(min(F(row["collision"]) for row in positive)), str(max(F(row["collision"]) for row in positive))] if positive else None,
                          "max_DP_denominator": max((row["subset_dp"]["denominator"] for row in positive), default=0),
                          "max_DP_states": max((len(row["subset_dp"]["states"]) for row in positive), default=0),
                          "max_captured_cells": max((len(row["captured_cells"]) for row in positive), default=0)})
        print(json.dumps({"done": case["case_id"], "local_checks": len(local)}), flush=True)
    result = {"status": "complete_registered_postprocessing", "elapsed_seconds": time.monotonic() - start,
              "registration_sha256": digest(REG), "script_sha256": digest(__file__),
              "summaries": summaries, "profiles": manifest, "log_two_interval": list(map(str, LOG_TWO)),
              "scope": reg["scope"], "source_receiver_threshold": reg["tau"],
              "no_rem maximization": True, "old_oracle_calls": 0}
    resultpath = BASE / (PREFIX + "_results.json")
    save(resultpath, result)
    receipt = BASE / (PREFIX + "_receipt.json")
    save(receipt, {"status": "complete", "exit_status": 0, "main_runs": 1,
                   "registration": {"path": str(REG), "sha256": digest(REG)},
                   "script": {"path": __file__, "sha256": digest(__file__)},
                   "results": {"path": str(resultpath), "sha256": digest(resultpath)},
                   "input_profiles": reg["input_profiles"], "profiles": manifest,
                   "old_oracle_calls": 0, "new_maximizations": 0,
                   "elapsed_seconds": result["elapsed_seconds"], "scope": reg["scope"]})
    print(json.dumps({"results_sha256": digest(resultpath), "receipt_sha256": digest(receipt)}))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("mode", choices=["register", "run"])
    mode = parser.parse_args().mode
    register() if mode == "register" else run()
