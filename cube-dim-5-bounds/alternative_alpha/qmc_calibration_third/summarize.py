#!/usr/bin/env python3
"""Read-only combination of six completed third-calibration files."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np


Q = 10.0 / 7.0
SEEDS = (2026102801, 2026102802)
SOURCE_DIR = Path(__file__).parents[1] / "qmc_face_i_second"


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    labels = ("old_sample_min", "old_sample_max", "old_sample_median")
    calibration = {}
    source_rows = []
    for label in labels:
        rows = [json.loads(Path(__file__).with_name(f"{label}_conditional_seed_{seed}.json").read_text())
                for seed in SEEDS]
        if not all(row["status"] == "complete_noncertified" for row in rows):
            raise RuntimeError(f"incomplete calibration for {label}")
        source = rows[0]["source_selection"]
        if any(row["source_selection"] != source for row in rows[1:]):
            raise RuntimeError(f"source disagreement for {label}")
        g48 = np.asarray([row["green_48_proxy"] for row in rows])
        g24 = np.asarray([row["green_24_subsample_proxy"] for row in rows])
        calibration[label] = {
            "source_selection": source,
            "old_12_proxy": source["truncated_green_proxy"],
            "green_48_blocks": [float(x) for x in g48], "green_24_blocks": [float(x) for x in g24],
            "green_48_mean": float(g48.mean()), "green_48_block_difference": float(abs(g48[0] - g48[1])),
            "green_24_mean": float(g24.mean()),
            "48_minus_24_mean": float(g48.mean() - g24.mean()),
            "48_minus_old12_mean": float(g48.mean() - source["truncated_green_proxy"]),
            "48_over_old12": float(g48.mean() / source["truncated_green_proxy"]),
            "tail_density_evaluations_by_block": [row["tail_density_evaluations"] for row in rows],
            "elapsed_seconds_by_block": [row["elapsed_seconds"] for row in rows],
        }
    for number in (1, 2):
        batch_path = SOURCE_DIR / f"face_batch_{number}_results.json"
        batch = json.loads(batch_path.read_text())
        source_rows.extend(batch["points"])
    all_old = np.asarray([row["truncated_green_proxy"] for row in source_rows])
    old_i = float(np.mean(all_old ** Q))
    old_max = calibration["old_sample_max"]["old_12_proxy"]
    new_max = calibration["old_sample_max"]["green_48_mean"]
    replacement_i = float(old_i + (new_max ** Q - old_max ** Q) / len(all_old))
    result = {
        "status": "read_only_noncertified_third_calibration_summary",
        "q": "10/7", "point_calibrations": calibration,
        "max_point_I_diagnostic": {
            "old_128_point_I_proxy": old_i,
            "old_max_Gq_share_of_sum": float(old_max ** Q / np.sum(all_old ** Q)),
            "replace_only_max_with_48_mean_I_proxy": replacement_i,
            "I_change": replacement_i - old_i,
            "relative_I_change": replacement_i / old_i - 1.0,
        },
        "min_stability_statement": "The old sampled-min point changed under time/sample calibration. This says nothing rigorous about g: all values are truncated/interpolated/QMC proxies with no certified ordering relative to essential infimum.",
        "source_face_batch_hashes": {f"face_batch_{number}": sha(SOURCE_DIR / f"face_batch_{number}_results.json")
                                     for number in (1, 2)},
        "limitations": ["No tail in t has been bounded.", "No stable-density interpolation or QMC error is bounded.",
                        "24-node result is a correlated subsample of 48 nodes, not an independent quadrature.",
                        "No global g, strict I, or C5 conclusion follows."],
        "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    }
    raw = json.dumps(result, indent=2, sort_keys=True).encode() + b"\n"
    Path(__file__).with_name("calibration_summary.json").write_bytes(raw)
    print(raw.decode(), end="")
    print("sha256", hashlib.sha256(raw).hexdigest())


if __name__ == "__main__":
    main()
