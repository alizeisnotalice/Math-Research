"""Fixed physical-cell heterogeneous retention patterns, exact scalar DP.

AST-only reuse of audited logarithm helpers; no old oracle/max execution.
"""
from fractions import Fraction as F
from pathlib import Path
import argparse
import ast
import datetime
import gzip
import hashlib
import json
import time

BASE = Path(__file__).resolve().parent
PREFIX = "source_cell_heterogeneous_clock_guard_20261007"
REG = BASE / (PREFIX + "_registration.json")
OLD = "spatial_packet_collision_probe_20261007"
HELPER = BASE / (OLD + ".py")
RETENTIONS = [F(1), F(1, 2), F(1, 4)]
QS = [2, 4, 16]


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def canonical(obj):
    return json.dumps(obj, sort_keys=True, separators=(",", ":")).encode()


def save(path, obj):
    with Path(path).open("x") as handle:
        json.dump(obj, handle, indent=2, sort_keys=True)
        handle.write("\n")


names = ["outward", "atanh_log_interval", "log_interval", "log_plus_interval"]
selected = [node for node in ast.parse(HELPER.read_text()).body
            if isinstance(node, ast.FunctionDef) and node.name in names]
assert [node.name for node in selected] == names
helper_ast = ast.Module(body=selected, type_ignores=[])
HELPER_AST_SHA = hashlib.sha256(ast.dump(helper_ast, include_attributes=False).encode()).hexdigest()
namespace = {"F": F, "LOG_BITS": 192, "LOG_TERMS": 100, "UNIT": 2 ** 192}
exec(compile(helper_ast, str(HELPER), "exec"), namespace)
namespace["LOG_TWO"] = namespace["atanh_log_interval"](F(2))
namespace["LOG_CACHE"] = {F(1): (F(0), F(0)), F(2): namespace["LOG_TWO"]}
log_interval = namespace["log_interval"]
log_plus_interval = namespace["log_plus_interval"]


def register():
    oldresult = BASE / (OLD + "_results.json")
    result = json.loads(oldresult.read_text())
    assert digest(HELPER) == result["script_sha256"]
    plans = []
    for item in result["profiles"]:
        assert digest(item["path"]) == item["sha256"]
        data = json.loads(gzip.decompress(Path(item["path"]).read_bytes()))
        fullcells = data["records"][0]["full_cells"]
        assert all(row["full_cells"] == fullcells for row in data["records"])
        assignments = []
        for shift in range(3):
            mapping = [{"index": cell["index"], "labels": cell["labels"],
                        "full_weight": cell["weight"], "residue": cell["index"][0] % 3,
                        "p": str(RETENTIONS[(cell["index"][0] + shift) % 3])}
                       for cell in fullcells]
            assert sum(F(cell["full_weight"]) for cell in mapping) == 1
            assignments.append({"pattern": shift, "assignment": mapping})
        plans.append({"case": data["case"], "input_profile": item,
                      "original_full_cells": fullcells, "patterns": assignments})
    save(REG, {"status": "registered_before_heterogeneous_clock_DP",
               "utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
               "script_sha256": digest(__file__), "plans": plans,
               "input_results": {"path": str(oldresult), "sha256": digest(oldresult)},
               "helper": {"path": str(HELPER), "sha256": digest(HELPER), "selected_AST_sha256": HELPER_AST_SHA,
                          "functions": names, "atanh_terms": 100, "outward_bits": 192},
               "pattern_rule": "For shift=0,1,2, pC=[1,1/2,1/4][(first full cell index+shift) mod3]; Python floor/mod convention handles negative cells. Complete source-cell assignment is frozen once per case and identical for all receivers.",
               "q": QS, "rounds": "Original r1 n1/2/3, r2 n4/8, r3 n16/64/256; all 15 cases/74 receivers retained.",
               "coin_law": "Independent complete-cell etaC~Bernoulli(pC); nu=sum etaC*muC/pC; original source and original winner untouched. Given fixed pi, V=sum piC*etaC/pC-1.",
               "exact_DP": "Original shares uC/T, each accepted response shift uC/pC is integer; probability update accept numerator/reject denominator-minus-numerator, with common product denominator. Marginalize only unselected coins, retain full assignment.",
               "checks": "EV=0; EV2=sum piC^2*(1/pC-1); Eell<= (2+4logq)*EV2, ell=logq-log_+(q*(1+V)).",
               "variance_zero": "All captured pC=1 implies V=ell=0 exactly; certify equality algebraically, not by subtraction of unrelated log intervals. For each outcome with q*(1+V)=q, ell=0 exact also.",
               "budget": {"cases": 15, "receivers": 74, "nonempty_expected": 64, "DP_patterns": 192, "q_checks": 576},
               "seeds": "No new random draws; source seeds inherited; exact coin law integrated.",
               "scope": "Local fixed-pattern scalar gate only; not adaptive policy simulation, not 22-mass-work global integration, original geom qualification, or closed upper bound.",
               "stopping": "One complete fixed pass, no retuning or old oracle; empty/top/zero variance/negative signed loss all retained."})
    print(json.dumps({"registration_sha256": digest(REG)}))


def pattern_row(row, mapping, shift):
    pi = list(map(F, row["posterior_cell"]))
    cells = row["captured_cells"]
    assert [F(cell["pi"]) for cell in cells] == pi
    ps = [mapping[tuple(cell["index"])] for cell in cells]
    T, units = row["subset_dp"]["denominator"], row["subset_dp"]["integer_cell_mass"]
    assert [F(u, T) for u in units] == pi
    dp, probden = {0: 1}, 1
    shifts = []
    for unit, p in zip(units, ps):
        shift_units = F(unit) / p
        assert shift_units.denominator == 1
        offset = int(shift_units)
        shifts.append(offset)
        A, B = p.numerator, p.denominator
        updated = {}
        for score, count in dp.items():
            if B > A:
                updated[score] = updated.get(score, 0) + (B - A) * count
            updated[score + offset] = updated.get(score + offset, 0) + A * count
        dp, probden = updated, probden * B
    assert sum(dp.values()) == probden
    mean = sum((F(count, probden) * (F(score, T) - 1) for score, count in dp.items()), F(0))
    variance = sum((F(count, probden) * (F(score, T) - 1) ** 2 for score, count in dp.items()), F(0))
    expected_variance = sum((pc * pc * (1 / p - 1) for pc, p in zip(pi, ps)), F(0))
    assert mean == 0 and variance == expected_variance
    if variance == 0:
        assert all(p == 1 for p in ps) and dp == {T: probden}
    checks = []
    for q in QS:
        zlo, zhi = log_interval(F(q))
        elllo, ellhi = F(0), F(0)
        exact_identity_outcome_probability = F(0)
        for score, count in dp.items():
            arg = F(q * score, T)
            probability = F(count, probden)
            if arg == q:
                lo, hi = F(0), F(0)
                exact_identity_outcome_probability += probability
            else:
                lowlog, highlog = log_plus_interval(arg)
                lo, hi = zlo - highlog, zhi - lowlog
            elllo += probability * lo
            ellhi += probability * hi
        rhslo, rhshi = (2 + 4 * zlo) * variance, (2 + 4 * zhi) * variance
        checks.append({"q": q, "expected_signed_ell_interval": [str(elllo), str(ellhi)],
                       "RHS_interval": [str(rhslo), str(rhshi)],
                       "gap_interval": [str(rhslo - ellhi), str(rhshi - elllo)],
                       "exact_identity_outcome_probability": str(exact_identity_outcome_probability),
                       "interval_pass": ellhi <= rhslo, "algebraic_zero_variance": variance == 0})
    return {"pattern": shift, "captured_cell_indices": [cell["index"] for cell in cells],
            "captured_pC": list(map(str, ps)), "posterior_pi": row["posterior_cell"],
            "EV": "0", "EV_squared": str(variance), "variance_formula": str(expected_variance),
            "zero_variance": variance == 0, "share_denominator": T,
            "accepted_response_integer_shifts": shifts,
            "DP_probability_denominator": probden, "DP_states": sorted([score, count] for score, count in dp.items()),
            "q_checks": checks}


def run():
    start = time.monotonic()
    reg = json.loads(REG.read_text())
    assert reg["script_sha256"] == digest(__file__)
    assert reg["helper"]["sha256"] == digest(HELPER)
    assert reg["helper"]["selected_AST_sha256"] == HELPER_AST_SHA
    assert digest(reg["input_results"]["path"]) == reg["input_results"]["sha256"]
    manifest, summaries = [], []
    for plan in reg["plans"]:
        item, case = plan["input_profile"], plan["case"]
        assert digest(item["path"]) == item["sha256"]
        old = json.loads(gzip.decompress(Path(item["path"]).read_bytes()))
        assert old["case"] == case
        maps = [{tuple(cell["index"]): F(cell["p"]) for cell in pattern["assignment"]} for pattern in plan["patterns"]]
        records = []
        for j, row in enumerate(old["records"]):
            assert row["full_cells"] == plan["original_full_cells"]
            identity = {"row": j, "source_row_sha256": hashlib.sha256(canonical(row)).hexdigest(),
                        "branch": row["branch"], "x": row["x"], "D": row["D"], "M": row["M"], "R": row["R"]}
            if row["branch"] == "empty":
                records.append(dict(identity, patterns=[]))
            else:
                records.append(dict(identity, patterns=[pattern_row(row, mapping, shift) for shift, mapping in enumerate(maps)]))
        path = BASE / (PREFIX + "_" + case["case_id"] + ".json.gz")
        with path.open("xb") as raw:
            with gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as zipped:
                zipped.write(canonical({"case": case, "input_profile": item, "complete_cell_patterns": plan["patterns"], "records": records}))
        manifest.append({"path": str(path), "sha256": digest(path)})
        patterns = [p for row in records for p in row["patterns"]]
        checks = [q for p in patterns for q in p["q_checks"]]
        nondeg = [q for p in patterns if not p["zero_variance"] for q in p["q_checks"]]
        summaries.append({"case_id": case["case_id"], "round": case["round"], "n": case["n"], "N": case["N"],
                          "receivers": len(records), "empty": sum(row["branch"] == "empty" for row in records),
                          "top": sum(row["branch"] == "top" for row in records),
                          "pattern_laws": len(patterns), "EV_and_variance_identity_passes": len(patterns),
                          "zero_variance_patterns": sum(p["zero_variance"] for p in patterns),
                          "q_checks": len(checks), "q_interval_passes": sum(q["interval_pass"] for q in checks),
                          "q_unresolved": sum(not q["interval_pass"] for q in checks),
                          "minimum_nondegenerate_gap_lower": str(min(F(q["gap_interval"][0]) for q in nondeg)) if nondeg else None,
                          "negative_signed_expectations": sum(F(q["expected_signed_ell_interval"][1]) < 0 for q in checks),
                          "max_DP_states": max((len(p["DP_states"]) for p in patterns), default=0)})
        print(json.dumps({"done": case["case_id"], "checks": len(checks)}), flush=True)
    resultpath = BASE / (PREFIX + "_results.json")
    save(resultpath, {"status": "complete_registered_fixed_heterogeneous_patterns", "summaries": summaries, "profiles": manifest,
                      "elapsed_seconds": time.monotonic() - start, "registration_sha256": digest(REG), "script_sha256": digest(__file__),
                      "scope": reg["scope"], "old_oracle_calls": 0, "new_maximizations": 0})
    receipt = BASE / (PREFIX + "_receipt.json")
    save(receipt, {"status": "complete", "exit_status": 0, "main_runs": 1,
                   "registration": {"path": str(REG), "sha256": digest(REG)}, "script": {"path": __file__, "sha256": digest(__file__)},
                   "results": {"path": str(resultpath), "sha256": digest(resultpath)}, "input_profiles": [p["input_profile"] for p in reg["plans"]],
                   "profiles": manifest, "old_oracle_calls": 0, "new_maximizations": 0, "scope": reg["scope"]})
    print(json.dumps({"results_sha256": digest(resultpath), "receipt_sha256": digest(receipt)}))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("mode", choices=["register", "run"])
    mode = parser.parse_args().mode
    register() if mode == "register" else run()
