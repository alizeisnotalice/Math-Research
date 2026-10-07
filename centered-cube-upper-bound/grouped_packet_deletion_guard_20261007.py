"""Actual complete packet deletion: exact full-space finite-J cell integrals.

Frozen arrangement/captures only; new counterfactual J winners, no old oracle.
"""
from collections import defaultdict
from fractions import Fraction as F
from pathlib import Path
from math import prod
import argparse
import ast
import datetime
import gzip
import hashlib
import json
import time

BASE = Path(__file__).resolve().parent
PREFIX = "grouped_packet_deletion_guard_20261007"
REG = BASE / (PREFIX + "_registration.json")
OLD = "log_max_variation_probe_20261007"
HELPER = BASE / "spatial_packet_collision_probe_20261007.py"


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


def add(coefs, arg, coef):
    if arg > 1 and coef != 0:
        coefs[arg] = coefs.get(arg, F(0)) + coef
        if coefs[arg] == 0:
            del coefs[arg]


def combine(target, source, factor=F(1)):
    for arg, coef in source.items():
        add(target, arg, factor * coef)


def expression(coefs, constant=F(0)):
    lo = hi = constant
    for arg, coef in sorted(coefs.items()):
        low, high = log_interval(arg)
        lo += coef * (low if coef > 0 else high)
        hi += coef * (high if coef > 0 else low)
    assert lo <= hi
    return {"constant": str(constant),
            "signed_log_coefficients": [{"argument": str(a), "coefficient": str(c)} for a, c in sorted(coefs.items())],
            "interval": [str(lo), str(hi)],
            "sign": "exact_zero" if lo == hi == 0 else "strict_negative" if hi < 0 else "strict_positive" if lo > 0 else "unresolved"}


def groupings(source):
    atoms = [tuple(map(F, y)) for y in source["atoms"]]
    weights = list(map(F, source["weights"]))
    locations = defaultdict(list)
    for j, y in enumerate(atoms):
        locations[y].append(j)
    packets, packet_ids = [], {}
    def identity(labels):
        key = tuple(sorted(labels))
        if key not in packet_ids:
            keyid = f"B{len(packets)}"
            packet_ids[key] = keyid
            packets.append({"id": keyid, "labels": list(key),
                            "full_packet_mass": str(sum((weights[j] for j in key), F(0)))})
        return packet_ids[key]
    families = [
        {"id": "physical_atoms", "packet_ids": [identity(js) for point, js in sorted(locations.items())]},
        {"id": "first_coordinate_left_right", "cut": "2", "packet_ids": [identity([j for j, y in enumerate(atoms) if y[0] < 2]), identity([j for j, y in enumerate(atoms) if y[0] >= 2])]},
        {"id": "all_source", "packet_ids": [identity(range(len(atoms)))]},
    ]
    byid = {packet["id"]: packet for packet in packets}
    for family in families:
        labels = [j for keyid in family["packet_ids"] for j in byid[keyid]["labels"]]
        assert sorted(labels) == list(range(len(atoms)))
        for js in locations.values():
            assert sum(set(js) <= set(byid[keyid]["labels"]) for keyid in family["packet_ids"]) == 1
    return packets, families


def cuts_for(source):
    atoms = [[F(v) for v in y] for y in source["atoms"]]
    scales = list(map(F, source["J"]))
    return [sorted({y[i] + sign * R / 2 for y in atoms for R in scales for sign in [-1, 1]}) for i in range(source["n"])]


def register():
    receiptpath = BASE / (OLD + "_receipt.json")
    receipt = json.loads(receiptpath.read_text())
    plans = []
    for n in [1, 2, 3]:
        path = BASE / f"{OLD}_n{n}_cell_profiles.json.gz"
        assert digest(path) == receipt["file_sha256"][path.name]
        data = json.loads(gzip.decompress(path.read_bytes()))
        source = {key: data["input"][key] for key in ["n", "atoms", "weights", "J", "tau"]}
        packets, families = groupings(source)
        cuts = cuts_for(source)
        plans.append({"n": n, "input_profile": {"path": str(path), "sha256": digest(path)},
                      "source": source, "source_sha256": hashlib.sha256(canonical(source)).hexdigest(),
                      "packets": packets, "groupings": families, "cell_count": len(data["cells"]),
                      "coordinate_cuts": [list(map(str, axis)) for axis in cuts],
                      "bounding_rectangle": [[str(axis[0]), str(axis[-1])] for axis in cuts]})
    notes = list(BASE.glob("grouped_source_clock*.md"))
    save(REG, {"status": "registered_before_new_counterfactual_packet_winners",
               "utc": datetime.datetime.now(datetime.timezone.utc).isoformat(), "script_sha256": digest(__file__),
               "plans": plans, "old_receipt": {"path": str(receiptpath), "sha256": digest(receiptpath)},
               "helper": {"path": str(HELPER), "sha256": digest(HELPER), "selected_AST_sha256": HELPER_AST_SHA,
                          "functions": names, "atanh_terms": 100, "outward_bits": 192},
               "grouped_clock_notes_present": [{"path": str(path), "sha256": digest(path)} for path in notes],
               "group_rule": "Three fixed complete partitions: unique physical atom (coincident labels together), first coordinate <2 vs >=2, and all source. Uncaptured/external source mass retained except the chosen complete deletion. No renormalization after deletion.",
               "spatial_coverage": "Verify saved bounds form the complete Cartesian product of every atom+-R/2 cut, exact cell volumes, whole support rectangle. Deleted source subset adds no cut; both kernels zero exterior. Coordinate faces Lebesgue-null.",
               "winner": "Read frozen original M/R and full original captured labels per finite J. For each new complete packet deletion, subtract packet masses from all J queries, choose actual deleted-source finite-J max with smallest-R exact tie. No old oracle/main/MC calls.",
               "compression": "Share exact outcomes between identical full-capture-by-J patterns; save every original cell volume/bounds and pattern reference, and all packet counterfactual arrays in each pattern.",
               "tau": "Each old case tau frozen, common for original and every grouping/deletion; full original denominator and posterior always used.",
               "primary": "Phi=tau integral log_+(M/tau), I=tau|M>tau|. E(group)=sum complete deletion Phi losses-I. Compare with C=tau integral_E sum_B(1+4z)pi_B^2, z=log(M/tau). Preserve signed E and packet exact-zero losses.",
               "interval": "Compile signed rational log coefficient maps, merge identical arguments before audited outward log integration; negative coefficients reverse endpoints. Certify C-E via its own combined map. Preserve unresolved intervals and strict M=tau plateaus.",
               "budget": {"rounds": 3, "n": [1, 2, 3], "saved_cells": 10914, "grouping_integrals": 9,
                          "unique_packet_integrals": 23},
               "seeds": "No new samples or seeds; complete saved Lebesgue arrangement and exact deterministic counterfactuals.",
               "scope": "True full-space finite common-J deletion integrals, not continuous-window proof, free receiver gate, nearmax source, grouped-clock global budget, or sqrt order evidence.",
               "stopping": "One complete fixed pass through three old profiles; no parameter/source/group tuning, no extra old oracle rerun."})
    print(json.dumps({"registration_sha256": digest(REG)}))


def process(plan):
    data = json.loads(gzip.decompress(Path(plan["input_profile"]["path"]).read_bytes()))
    source = plan["source"]
    n, scales, tau = source["n"], list(map(F, source["J"])), F(source["tau"])
    weights = list(map(F, source["weights"]))
    assert all(w > 0 for w in weights) and sum(weights) == 1
    assert sorted(scales) == scales and scales[0] == 1
    cuts = [list(map(F, axis)) for axis in plan["coordinate_cuts"]]
    interval_maps = [{(axis[i], axis[i + 1]): i for i in range(len(axis) - 1)} for axis in cuts]
    expected_cells = prod(len(axis) - 1 for axis in cuts)
    assert expected_cells == len(data["cells"]) == plan["cell_count"]
    signatures, patterns, cellrows, seen = {}, [], [], set()
    areas, occurrences = defaultdict(F), defaultdict(int)
    for cellindex, cell in enumerate(data["cells"]):
        bounds = [tuple(map(F, pair)) for pair in cell["cell_bounds"]]
        key = tuple(interval_maps[i][pair] for i, pair in enumerate(bounds))
        assert key not in seen
        seen.add(key)
        volume = F(cell["volume"])
        assert volume == prod((right - left for left, right in bounds), start=F(1))
        assert all(left < F(mid) < right for (left, right), mid in zip(bounds, cell["midpoint"]))
        signature = tuple(tuple(js) for js in cell["full_capture_labels_by_scale"])
        if signature not in signatures:
            pid = len(patterns)
            signatures[signature] = pid
            masses = [sum((weights[j] for j in js), F(0)) for js in signature]
            responses = [m / R ** n for m, R in zip(masses, scales)]
            assert list(map(str, masses)) == cell["full_mass_by_scale"]
            assert list(map(str, responses)) == cell["full_response_by_scale"]
            M, R = F(cell["M"]), F(cell["winner_R"])
            win = scales.index(R)
            assert responses[win] == M and all(value <= M for value in responses)
            assert win == next(j for j, value in enumerate(responses) if value == M)
            assert cell["E"] == (M > tau) and cell["threshold_plateau"] == (M == tau)
            counter = {}
            for packet in plan["packets"]:
                removed = set(packet["labels"])
                remaining = [[j for j in js if j not in removed] for js in signature]
                newmasses = [sum((weights[j] for j in js), F(0)) for js in remaining]
                newresponse = [m / r ** n for m, r in zip(newmasses, scales)]
                Md = max(newresponse)
                newwin = next(j for j, value in enumerate(newresponse) if value == Md)
                assert Md <= M
                pi = sum((weights[j] for j in signature[win] if j in removed), F(0)) / masses[win] if M > 0 else None
                counter[packet["id"]] = {
                    "remaining_capture_labels_by_J": remaining, "mass_by_J": list(map(str, newmasses)),
                    "response_by_J": list(map(str, newresponse)), "M_deleted": str(Md),
                    "winner_R_deleted": str(scales[newwin]), "winner_tie_count": sum(v == Md for v in newresponse),
                    "zero_response": Md == 0, "E_deleted": Md > tau, "threshold_plateau_deleted": Md == tau,
                    "winner_switch": scales[newwin] != R,
                    "positive_response_winner_switch": Md > 0 and scales[newwin] != R,
                    "original_winner_pi_B": str(pi) if pi is not None else None,
                }
            families = []
            for family in plan["groupings"]:
                pis = [F(counter[keyid]["original_winner_pi_B"]) for keyid in family["packet_ids"]] if M > 0 else []
                if M > 0:
                    assert sum(pis) == 1 and all(0 <= p <= 1 for p in pis)
                families.append({"id": family["id"], "pi": list(map(str, pis)),
                                 "collision": str(sum((p * p for p in pis), F(0))) if pis else None})
            patterns.append({"id": pid, "capture_labels_by_J": [list(js) for js in signature],
                             "original_mass_by_J": list(map(str, masses)), "original_response_by_J": list(map(str, responses)),
                             "M": str(M), "R": str(R), "E": M > tau, "threshold_plateau": M == tau,
                             "original_winner_tie_count": cell["maximum_tie_count"],
                             "packet_counterfactuals": counter, "group_posteriors": families})
        pid = signatures[signature]
        pattern = patterns[pid]
        assert pattern["M"] == cell["M"] and pattern["R"] == cell["winner_R"]
        assert pattern["E"] == cell["E"]
        areas[pid] += volume
        occurrences[pid] += 1
        cellrows.append({"original_cell_index": cellindex, "bounds": cell["cell_bounds"],
                         "volume": str(volume), "midpoint": cell["midpoint"], "pattern_id": pid})
    assert len(seen) == expected_cells
    volume_sum = sum(areas.values(), F(0))
    support_rect_volume = prod((axis[-1] - axis[0] for axis in cuts), start=F(1))
    assert volume_sum == support_rect_volume
    phi, packet_phi = {}, {packet["id"]: {} for packet in plan["packets"]}
    I = F(0)
    for pattern in patterns:
        area, M = areas[pattern["id"]], F(pattern["M"])
        pattern["Lebesgue_volume"] = str(area)
        pattern["cell_count"] = occurrences[pattern["id"]]
        if M > tau:
            I += tau * area
            add(phi, M / tau, tau * area)
        for keyid, outcome in pattern["packet_counterfactuals"].items():
            Md = F(outcome["M_deleted"])
            if Md > tau:
                add(packet_phi[keyid], Md / tau, tau * area)
    packet_results, loss_maps = [], {}
    for packet in plan["packets"]:
        keyid = packet["id"]
        loss = dict(phi)
        combine(loss, packet_phi[keyid], -F(1))
        loss_maps[keyid] = loss
        expr = expression(loss)
        packet_results.append({"packet": packet, "remaining_total_mass": str(1 - F(packet["full_packet_mass"])),
                               "Phi_deleted": expression(packet_phi[keyid]), "actual_deletion_loss": expr,
                               "nonnegative_loss_certified": F(expr["interval"][0]) >= 0,
                               "original_E_winner_switch_cell_count": sum(occurrences[p["id"]] for p in patterns if p["E"] and p["packet_counterfactuals"][keyid]["winner_switch"]),
                               "positive_deleted_response_E_winner_switch_cell_count": sum(occurrences[p["id"]] for p in patterns if p["E"] and p["packet_counterfactuals"][keyid]["positive_response_winner_switch"]),
                               "original_E_zero_response_cell_count": sum(occurrences[p["id"]] for p in patterns if p["E"] and p["packet_counterfactuals"][keyid]["zero_response"]),
                               "original_E_pi_zero_cells": sum(occurrences[p["id"]] for p in patterns if p["E"] and p["packet_counterfactuals"][keyid]["original_winner_pi_B"] == "0"),
                               "original_E_pi_one_cells": sum(occurrences[p["id"]] for p in patterns if p["E"] and p["packet_counterfactuals"][keyid]["original_winner_pi_B"] == "1")})
    family_results = []
    for family in plan["groupings"]:
        total_loss = {}
        for keyid in family["packet_ids"]:
            combine(total_loss, loss_maps[keyid])
        scalar, scalar_const = {}, F(0)
        for pattern in patterns:
            if not pattern["E"]:
                continue
            collision = F(next(f for f in pattern["group_posteriors"] if f["id"] == family["id"])["collision"])
            coefficient = tau * areas[pattern["id"]] * collision
            scalar_const += coefficient
            add(scalar, F(pattern["M"]) / tau, 4 * coefficient)
        excess = expression(total_loss, -I)
        scalar_expr = expression(scalar, scalar_const)
        gapmap = dict(scalar)
        combine(gapmap, total_loss, -F(1))
        gap = expression(gapmap, scalar_const + I)
        family_results.append({"grouping": family, "sum_actual_deletion_loss": expression(total_loss),
                               "I": str(I), "actual_excess_sum_D_minus_I": excess,
                               "scalar_tau_integral_sum_1plus4z_pi_squared": scalar_expr,
                               "scalar_minus_actual_excess": gap,
                               "upper_comparison_interval_pass": F(gap["interval"][0]) >= 0})
    summary = {"n": n, "a": "1", "J": source["J"], "tau": str(tau), "W": "1",
               "cells": len(cellrows), "patterns": len(patterns), "support_rectangle_volume": str(volume_sum),
               "original_E_cells": sum(p["cell_count"] for p in patterns if p["E"]),
               "threshold_plateau_cells": sum(p["cell_count"] for p in patterns if p["threshold_plateau"]),
               "original_positive_response_tie_cells": sum(p["cell_count"] for p in patterns if F(p["M"]) > 0 and p["original_winner_tie_count"] > 1),
               "unique_packets": len(plan["packets"]), "new_counterfactual_patterns": len(patterns) * len(plan["packets"]),
               "counterfactual_scale_responses": len(patterns) * len(plan["packets"]) * len(scales),
               "I": str(I), "Phi_original": expression(phi),
               "packet_nonnegative_loss_passes": sum(p["nonnegative_loss_certified"] for p in packet_results),
               "grouping_upper_comparison_passes": sum(f["upper_comparison_interval_pass"] for f in family_results),
               "family_results": family_results}
    return {"plan": plan, "coverage": {"whole_cartesian_support_certified": True, "exact_volume_sum": str(volume_sum),
                                       "exterior_original_and_deleted_responses": "0", "faces": "Lebesgue-null coordinate hyperplanes"},
            "cells": cellrows, "patterns": patterns, "packet_results": packet_results, "summary": summary}


def run():
    start = time.monotonic()
    reg = json.loads(REG.read_text())
    assert reg["script_sha256"] == digest(__file__)
    assert reg["helper"]["sha256"] == digest(HELPER) and reg["helper"]["selected_AST_sha256"] == HELPER_AST_SHA
    manifests, summaries = [], []
    for plan in reg["plans"]:
        assert digest(plan["input_profile"]["path"]) == plan["input_profile"]["sha256"]
        assert hashlib.sha256(canonical(plan["source"])).hexdigest() == plan["source_sha256"]
        result = process(plan)
        path = BASE / f"{PREFIX}_n{plan['n']}_cells.json.gz"
        with path.open("xb") as raw:
            with gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as zipped:
                zipped.write(canonical(result))
        manifests.append({"path": str(path), "sha256": digest(path)})
        summaries.append(result["summary"])
        print(json.dumps({"done_n": plan["n"], "patterns": result["summary"]["patterns"],
                          "grouping_passes": result["summary"]["grouping_upper_comparison_passes"]}), flush=True)
    output = BASE / (PREFIX + "_results.json")
    save(output, {"status": "complete_registered_full_space_finite_J_packet_deletion", "summaries": summaries,
                  "profiles": manifests, "elapsed_seconds": time.monotonic() - start,
                  "registration_sha256": digest(REG), "script_sha256": digest(__file__), "scope": reg["scope"],
                  "old_oracle_calls": 0, "new_packet_counterfactual_J_winners": True})
    receipt = BASE / (PREFIX + "_receipt.json")
    save(receipt, {"status": "complete", "exit_status": 0, "main_runs": 1,
                   "script": {"path": __file__, "sha256": digest(__file__)},
                   "registration": {"path": str(REG), "sha256": digest(REG)},
                   "results": {"path": str(output), "sha256": digest(output)},
                   "input_profiles": [p["input_profile"] for p in reg["plans"]], "profiles": manifests,
                   "old_oracle_calls": 0, "scope": reg["scope"]})
    print(json.dumps({"results_sha256": digest(output), "receipt_sha256": digest(receipt)}))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("mode", choices=["register", "run"])
    mode = parser.parse_args().mode
    register() if mode == "register" else run()
