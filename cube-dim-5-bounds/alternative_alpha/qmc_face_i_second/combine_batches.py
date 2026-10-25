#!/usr/bin/env python3
"""Combine the two already-recorded non-certified face-QMC batches; no rerun."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np


Q = 10.0 / 7.0


def main():
    files = [Path(__file__).with_name(f"face_batch_{i}_results.json") for i in (1, 2)]
    batches = [json.loads(path.read_text()) for path in files]
    if not all(item["status"] == "complete_noncertified" and len(item["points"]) == 64 for item in batches):
        raise RuntimeError("both fixed batches must be complete before combining")
    green = np.asarray([point["truncated_green_proxy"] for batch in batches for point in batch["points"]])
    per_batch_i = np.asarray([batch["aggregate"]["I_proxy"] for batch in batches])
    combined_i = float(np.mean(green ** Q))
    sampled_g = float(np.min(green))
    cost = float((10.0 / 3.0) ** Q * combined_i / sampled_g ** Q)
    result = {
        "status": "complete_noncertified_combination_of_two_fixed_face_batches",
        "batch_files": [path.name for path in files],
        "point_count": int(len(green)),
        "q": "10/7",
        "I_proxy": combined_i,
        "sampled_minimum_G_proxy": sampled_g,
        "conditional_cost_proxy": cost,
        "independent_batch_I_proxies": [float(x) for x in per_batch_i],
        "two_batch_I_sample_std": float(per_batch_i.std(ddof=1)),
        "two_batch_I_standard_error_diagnostic": float(per_batch_i.std(ddof=1) / np.sqrt(2)),
        "independent_batch_cost_proxies": [batch["aggregate"]["conditional_cost_proxy"] for batch in batches],
        "limitations": [
            "Both I and g are finite-QMC, finite-time, interpolated/tailed numerical proxies.",
            "The truncated/interpolated/QMC sampled minimum has no certified ordering relative to the unknown essential infimum g.",
            "The two-batch variation is not a rigorous confidence interval.",
            "This is neither a global g certificate nor a strict I or C5 estimate.",
        ],
        "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    }
    raw = json.dumps(result, indent=2, sort_keys=True).encode() + b"\n"
    Path(__file__).with_name("combined_results.json").write_bytes(raw)
    print(raw.decode(), end="")
    print("sha256", hashlib.sha256(raw).hexdigest())


if __name__ == "__main__":
    main()
