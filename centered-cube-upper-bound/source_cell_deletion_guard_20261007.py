"""New general-p coin and deletion scalar guards on frozen cell posteriors.

Only four audited log helper definitions are AST-loaded from the previous file;
its module/main, maximization, and geometry routines are never executed.
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
PREFIX = "source_cell_deletion_guard_20261007"
REG = BASE / (PREFIX + "_registration.json")
OLD = "spatial_packet_collision_probe_20261007"
HELPER_SOURCE = BASE / (OLD + ".py")
PS = [F(1, 4), F(3, 4), F(7, 8)]
QS = [2, 4, 16]
HELPER_NAMES = ["outward", "atanh_log_interval", "log_interval", "log_plus_interval"]


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def canonical(obj):
    return json.dumps(obj, sort_keys=True, separators=(",", ":")).encode()


def save(path, obj):
    with Path(path).open("x") as handle:
        json.dump(obj, handle, sort_keys=True, indent=2)
        handle.write("\n")


def load_only_log_helpers():
    tree = ast.parse(HELPER_SOURCE.read_text())
    definitions = [node for node in tree.body
                   if isinstance(node, ast.FunctionDef) and node.name in HELPER_NAMES]
    assert [node.name for node in definitions] == HELPER_NAMES
    helper_ast = ast.Module(body=definitions, type_ignores=[])
    helper_ast_hash = hashlib.sha256(ast.dump(helper_ast, include_attributes=False).encode()).hexdigest()
    namespace = {"F": F, "LOG_BITS": 192, "LOG_TERMS": 100, "UNIT": 2 ** 192}
    exec(compile(helper_ast, str(HELPER_SOURCE), "exec"), namespace)
    namespace["LOG_TWO"] = namespace["atanh_log_interval"](F(2))
    namespace["LOG_CACHE"] = {F(1): (F(0), F(0)), F(2): namespace["LOG_TWO"]}
    return namespace, helper_ast_hash


HELPERS, HELPER_AST_SHA = load_only_log_helpers()
log_interval = HELPERS["log_interval"]
log_plus_interval = HELPERS["log_plus_interval"]


def register():
    oldresults = BASE / (OLD + "_results.json")
    oldreg = BASE / (OLD + "_registration.json")
    result = json.loads(oldresults.read_text())
    assert digest(HELPER_SOURCE) == result["script_sha256"]
    for item in result["profiles"]:
        assert digest(item["path"]) == item["sha256"]
    save(REG, {
        "status": "registered_before_general_p_and_deletion_postprocessing",
        "utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "script_sha256": digest(__file__),
        "helper": {"path": str(HELPER_SOURCE), "sha256": digest(HELPER_SOURCE),
                   "selected_functions": HELPER_NAMES, "selected_AST_sha256": HELPER_AST_SHA,
                   "reuse_scope": "Only selected log definitions AST-executed in a new Fraction namespace; old main/module/response/geometry not executed.",
                   "log_parameters": {"atanh_terms": 100, "outward_dyadic_bits": 192,
                                      "strict_tail": "2*z^(201)/(201*(1-z^2))"}},
        "inputs": [{"path": str(oldresults), "sha256": digest(oldresults)},
                   {"path": str(oldreg), "sha256": digest(oldreg)}],
        "input_profiles": result["profiles"],
        "rounds": [{"round": 1, "n": [1, 2, 3], "receiver_count": 18},
                   {"round": 2, "n": [4, 8], "receiver_count": 24},
                   {"round": 3, "n": [16, 64, 256], "receiver_count": 32}],
        "p": list(map(str, PS)), "q": QS,
        "source": "Original full complete cell measures muC, all source labels retained; use frozen piC=muC(Q_R)/mR and chi=sum piC^2, never remaximize.",
        "coin_law": "Independent etaC~Bernoulli(p) per complete fixed cell, nu=(1/p)*sum etaC*muC. U=sum piC etaC, V=U/p-1, ell=log(q)-log_+(q*(1+V)), log_+(0)=0.",
        "weighted_DP": "For p=A/B and saved pi integer shares uC/T: reject multiplier B-A, accept A; probability numerator per score over B^K. Integrate uncaptured cells out, complete source still retained. Check E U=p, E V=0, E V^2=((1-p)/p)*chi exactly.",
        "coin_contract": "E ell <= (2+4*log(q))*((1-p)/p)*chi; compare expectation upper with RHS lower using audited rational log enclosures.",
        "deletion_contract": "sum_C min(log(q),-log(1-piC)) <= 1+(1+4*log(q))*chi. Handle pi=1 explicitly as log(q); otherwise use exact branch criterion pi<=1-1/q, no ambiguous floating log branch.",
        "budget": {"cases": 15, "receiver_count": 74, "nonempty_expected": 64,
                   "weighted_DP_laws": 192, "coin_q_checks": 576, "deletion_q_checks": 192},
        "seeds": "Frozen original construction seeds retained; no new random sampling, DP integrates full prescribed coin law exactly.",
        "stopping": "One fixed complete pass through all 15 frozen profiles; retain empty/top/unresolved/failures, no input tuning or old oracle execution.",
        "scope": "Finite posterior scalar general-p/deletion guards, not numerical proof of general lemma, 16W full-space integral, geom/near-optimal qualification, or dimension order.",
    })
    print(json.dumps({"registration": str(REG), "sha256": digest(REG)}))


def weighted_dp(units, denominator, p, collision):
    assert sum(units) == denominator
    A, B = p.numerator, p.denominator
    dp = {0: 1}
    for unit in units:
        updated = {}
        for score, weight in dp.items():
            updated[score] = updated.get(score, 0) + (B - A) * weight
            updated[score + unit] = updated.get(score + unit, 0) + A * weight
        dp = updated
    probability_denominator = B ** len(units)
    assert sum(dp.values()) == probability_denominator
    expected_U = sum((F(score * count, denominator * probability_denominator)
                      for score, count in dp.items()), F(0))
    expected_V = sum((F(count, probability_denominator) * (F(score, denominator) / p - 1)
                      for score, count in dp.items()), F(0))
    expected_V2 = sum((F(count, probability_denominator) * (F(score, denominator) / p - 1) ** 2
                       for score, count in dp.items()), F(0))
    assert expected_U == p and expected_V == 0
    assert expected_V2 == (1 - p) / p * collision
    return sorted(dp.items()), probability_denominator, expected_V2


def process_row(row, row_number):
    identity = {"source_row_number": row_number,
                "source_row_sha256": hashlib.sha256(canonical(row)).hexdigest(),
                "x": row["x"], "D": row["D"], "kind": row["kind"],
                "M": row["M"], "R": row["R"], "d": row["d"], "branch": row["branch"],
                "full_cells": row["full_cells"]}
    if row["branch"] == "empty":
        return dict(identity, posterior=None, p_checks=[], deletion_checks=[])
    pi = list(map(F, row["posterior_cell"]))
    assert sum(pi) == 1 and all(0 < value <= 1 for value in pi)
    collision = F(row["collision"])
    assert collision == sum(value * value for value in pi)
    T, units = row["subset_dp"]["denominator"], row["subset_dp"]["integer_cell_mass"]
    assert [F(unit, T) for unit in units] == pi
    pchecks = []
    for p in PS:
        dp, probability_denominator, variance = weighted_dp(units, T, p, collision)
        qchecks = []
        for q in QS:
            zlo, zhi = log_interval(F(q))
            loglo, loghi = F(0), F(0)
            for score, numerator in dp:
                low, high = log_plus_interval(F(q * score, T) / p)
                probability = F(numerator, probability_denominator)
                loglo += probability * low
                loghi += probability * high
            losslo, losshi = zlo - loghi, zhi - loglo
            rhslo, rhshi = (2 + 4 * zlo) * variance, (2 + 4 * zhi) * variance
            qchecks.append({"q": q, "original_row_tau": row["local_checks"][QS.index(q)]["tau"],
                            "log_q_interval": [str(zlo), str(zhi)],
                            "expected_log_plus_interval": [str(loglo), str(loghi)],
                            "expected_ell_interval": [str(losslo), str(losshi)],
                            "RHS_interval": [str(rhslo), str(rhshi)],
                            "strict_gap_interval": [str(rhslo - losshi), str(rhshi - losslo)],
                            "interval_pass": losshi <= rhslo})
        pchecks.append({"p": str(p), "variance_factor": str((1 - p) / p),
                        "expected_U": str(p), "expected_V": "0", "expected_V_squared": str(variance),
                        "DP_probability_denominator": probability_denominator,
                        "DP_states": [[score, numerator] for score, numerator in dp],
                        "q_checks": qchecks})
    deletion = []
    for q in QS:
        zlo, zhi = log_interval(F(q))
        sumlo, sumhi = F(0), F(0)
        terms = []
        for value in pi:
            if value == 1:
                branch, argument, low, high = "pi_equals_one", None, zlo, zhi
            elif value <= 1 - F(1, q):
                argument = 1 / (1 - value)
                low, high = log_interval(argument)
                branch = "untruncated_deletion_log"
            else:
                branch, argument, low, high = "truncated_log_q", None, zlo, zhi
            sumlo += low
            sumhi += high
            terms.append({"pi": str(value), "branch": branch,
                          "log_argument": str(argument) if argument is not None else None,
                          "term_interval": [str(low), str(high)]})
        rhslo, rhshi = 1 + (1 + 4 * zlo) * collision, 1 + (1 + 4 * zhi) * collision
        deletion.append({"q": q, "terms": terms, "deletion_sum_interval": [str(sumlo), str(sumhi)],
                         "RHS_interval": [str(rhslo), str(rhshi)],
                         "strict_gap_interval": [str(rhslo - sumhi), str(rhshi - sumlo)],
                         "interval_pass": sumhi <= rhslo})
    return dict(identity, captured_cells=row["captured_cells"], posterior_cell=row["posterior_cell"],
                collision=row["collision"], DP_share_denominator=T, integer_cell_shares=units,
                p_checks=pchecks, deletion_checks=deletion)


def run():
    start = time.monotonic()
    reg = json.loads(REG.read_text())
    assert reg["script_sha256"] == digest(__file__)
    assert reg["helper"]["sha256"] == digest(HELPER_SOURCE)
    assert reg["helper"]["selected_AST_sha256"] == HELPER_AST_SHA
    for item in reg["inputs"]:
        assert digest(item["path"]) == item["sha256"]
    summaries, manifests = [], []
    for item in reg["input_profiles"]:
        assert digest(item["path"]) == item["sha256"]
        old = json.loads(gzip.decompress(Path(item["path"]).read_bytes()))
        case = old["case"]
        records = [process_row(row, j) for j, row in enumerate(old["records"])]
        path = BASE / (PREFIX + "_" + case["case_id"] + ".json.gz")
        with path.open("xb") as raw:
            with gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as zipped:
                zipped.write(canonical({"case": case, "input_profile": item, "records": records}))
        manifests.append({"path": str(path), "sha256": digest(path)})
        nonempty = [row for row in records if row["branch"] != "empty"]
        coins = [c for row in nonempty for p in row["p_checks"] for c in p["q_checks"]]
        deletions = [c for row in nonempty for c in row["deletion_checks"]]
        summaries.append({"case_id": case["case_id"], "round": case["round"], "n": case["n"], "N": case["N"],
                          "seed": case["seed"], "source_sha256": case["source_sha256"],
                          "receiver_count": len(records), "nonempty": len(nonempty),
                          "empty": len(records) - len(nonempty), "top": sum(row["branch"] == "top" for row in records),
                          "coin_checks": len(coins), "coin_interval_passes": sum(c["interval_pass"] for c in coins),
                          "coin_unresolved": sum(not c["interval_pass"] for c in coins),
                          "deletion_checks": len(deletions), "deletion_interval_passes": sum(c["interval_pass"] for c in deletions),
                          "deletion_unresolved": sum(not c["interval_pass"] for c in deletions),
                          "pi_equals_one_receivers": sum(any(F(v) == 1 for v in row["posterior_cell"]) for row in nonempty),
                          "minimum_coin_gap_lower": str(min(F(c["strict_gap_interval"][0]) for c in coins)) if coins else None,
                          "minimum_deletion_gap_lower": str(min(F(c["strict_gap_interval"][0]) for c in deletions)) if deletions else None})
        print(json.dumps({"done": case["case_id"], "coin_checks": len(coins), "deletion_checks": len(deletions)}), flush=True)
    resultpath = BASE / (PREFIX + "_results.json")
    save(resultpath, {"status": "complete_registered_general_p_and_deletion_guard", "summaries": summaries,
                      "profiles": manifests, "elapsed_seconds": time.monotonic() - start,
                      "registration_sha256": digest(REG), "script_sha256": digest(__file__),
                      "helper": reg["helper"], "scope": reg["scope"],
                      "old_maximizations": 0, "old_oracle_calls": 0,
                      "evidence": "Exact weighted integer Bernoulli DP and 192bit outward rational log bounds; finite frozen posteriors only."})
    receipt = BASE / (PREFIX + "_receipt.json")
    save(receipt, {"status": "complete", "exit_status": 0, "main_runs": 1,
                   "registration": {"path": str(REG), "sha256": digest(REG)},
                   "script": {"path": __file__, "sha256": digest(__file__)},
                   "results": {"path": str(resultpath), "sha256": digest(resultpath)},
                   "input_profiles": reg["input_profiles"], "profiles": manifests,
                   "old_oracle_calls": 0, "new_maximizations": 0, "scope": reg["scope"]})
    print(json.dumps({"results_sha256": digest(resultpath), "receipt_sha256": digest(receipt)}))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("mode", choices=["register", "run"])
    mode = parser.parse_args().mode
    register() if mode == "register" else run()
