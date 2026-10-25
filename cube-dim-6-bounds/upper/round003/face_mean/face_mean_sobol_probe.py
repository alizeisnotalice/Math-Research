#!/usr/bin/env python3
"""Small deterministic feasibility screen for E[1/||X_1||_infty].

X_1 is the unit-time body32 Cauchy sum.  This is not an interval proof and
is not used as an inequality certificate.
"""
from __future__ import annotations

from hashlib import sha256
from itertools import product
from pathlib import Path
from time import perf_counter
import json
import numpy as np
from scipy.stats import qmc


def main() -> None:
    started = perf_counter()
    level = 15
    seeds = [20261025, 20261026, 20261027, 20261028]
    V = np.array([(1,) + s for s in product((-1, 1), repeat=5)], dtype=float)
    results = []
    for seed in seeds:
        u = qmc.Sobol(d=32, scramble=True, seed=seed).random_base2(level)
        z = np.tan(np.pi * (u - 0.5))
        x = z @ V
        radius = np.max(np.abs(x), axis=1)
        values = 1.0 / radius
        results.append({"seed": seed, "points": len(values), "mean": float(values.mean()), "standard_error_iid_display": float(values.std(ddof=1)/np.sqrt(len(values)))})
    means = np.array([r["mean"] for r in results])
    payload = {
        "status": "numerical_feasibility_probe_only",
        "model": "body32, X_1=sum_v v C_v",
        "batch": {"sobol_dimension": 32, "points_per_scramble": 2 ** level, "scrambles": seeds},
        "replicates": results,
        "replicate_mean": float(means.mean()),
        "replicate_standard_error": float(means.std(ddof=1)/np.sqrt(len(means))),
        "face_mean_estimate": float(means.mean()/384.0),
        "warning": "Scrambled-Sobol variation is diagnostic only, not a rigorous lower or upper bound.",
        "elapsed_seconds": perf_counter()-started,
        "source_sha256": sha256(Path(__file__).read_bytes()).hexdigest(),
    }
    Path(__file__).with_name("face_mean_sobol_probe.json").write_text(json.dumps(payload, indent=2, sort_keys=True)+"\n")
    print(json.dumps(payload, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
