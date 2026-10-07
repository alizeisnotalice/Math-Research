"""Read-only reconstruction of the new saved DP/geometry records, no oracle."""
from fractions import Fraction as F
from pathlib import Path
import gzip
import hashlib
import json

BASE = Path(__file__).resolve().parent
PREFIX = "spatial_packet_collision_probe_20261007"
resultpath = BASE / (PREFIX + "_results.json")
result = json.loads(resultpath.read_text())
counts = {
    "profile_hashes": 0, "complete_cell_weight_sums": 0,
    "posterior_cell_sums": 0, "exact_DP_V_mean_zero": 0,
    "exact_DP_variance_equals_collision": 0, "local_interval_checks": 0,
    "source_labels_on_grid_face": 0, "negative_cell_coordinate_labels": 0,
    "corner_outside_original_cube": 0, "partial_capture_cell": 0,
    "uncaptured_labels_in_captured_cells": 0, "max_DP_outcomes": 0,
}
min_gap = None
for item in result["profiles"]:
    path = Path(item["path"])
    assert hashlib.sha256(path.read_bytes()).hexdigest() == item["sha256"]
    counts["profile_hashes"] += 1
    data = json.loads(gzip.decompress(path.read_bytes()))
    case = data["case"]
    d = F(data["records"][0]["d"])
    for y in case["source"]["atoms"]:
        yy = list(map(F, y))
        counts["source_labels_on_grid_face"] += any((v / d).denominator == 1 for v in yy)
        counts["negative_cell_coordinate_labels"] += any(v // d < 0 for v in yy)
    for row in data["records"]:
        assert sum(F(cell["weight"]) for cell in row["full_cells"]) == 1
        labels = [j for cell in row["full_cells"] for j in cell["labels"]]
        assert sorted(labels) == list(range(case["N"]))
        counts["complete_cell_weight_sums"] += 1
        if row["branch"] == "empty":
            continue
        pi = list(map(F, row["posterior_cell"]))
        assert sum(pi) == 1
        counts["posterior_cell_sums"] += 1
        R = F(row["R"])
        for cell in row["captured_cells"]:
            counts["corner_outside_original_cube"] += F(cell["corner_r"]) > R
            counts["partial_capture_cell"] += len(cell["labels"]) != len(cell["captured_labels"])
            counts["uncaptured_labels_in_captured_cells"] += len(cell["labels"]) - len(cell["captured_labels"])
        dp = row["subset_dp"]
        denominator, outcomes = dp["denominator"], dp["outcomes"]
        assert sum(count for score, count in dp["states"]) == outcomes
        EV = sum(F(count, outcomes) * (F(2 * score, denominator) - 1) for score, count in dp["states"])
        EV2 = sum(F(count, outcomes) * (F(2 * score, denominator) - 1) ** 2 for score, count in dp["states"])
        assert EV == 0
        assert EV2 == F(row["collision"]) == sum(p * p for p in pi)
        counts["exact_DP_V_mean_zero"] += 1
        counts["exact_DP_variance_equals_collision"] += 1
        counts["max_DP_outcomes"] = max(counts["max_DP_outcomes"], outcomes)
        for check in row["local_checks"]:
            losslo, losshi = map(F, check["expected_loss_interval"])
            rhslo, rhshi = map(F, check["bound_interval"])
            assert losslo <= losshi <= rhslo <= rhshi
            gap = rhslo - losshi
            assert gap > 0
            min_gap = gap if min_gap is None else min(gap, min_gap)
            counts["local_interval_checks"] += 1
output = {
    "status": "saved_DP_reconstruction_passed", "counts": counts,
    "minimum_strict_local_gap_lower": str(min_gap),
    "input_results_sha256": hashlib.sha256(resultpath.read_bytes()).hexdigest(),
    "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    "scope": "Only saved Fraction/hash/DP/cell arithmetic, no old or new maximization/oracle.",
}
path = BASE / (PREFIX + "_saved_dp_audit.json")
with path.open("x") as handle:
    json.dump(output, handle, indent=2)
    handle.write("\n")
print(json.dumps({"counts": counts, "audit_sha256": hashlib.sha256(path.read_bytes()).hexdigest(), "min_gap_approx_diagnostic": float(min_gap)}))
