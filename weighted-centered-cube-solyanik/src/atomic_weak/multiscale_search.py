#!/usr/bin/env python3
"""Exploratory search beyond the archived <=8-atom bounded-family searches.

The candidates are symmetric geometric chains
    sum_{k=0}^{K-1} m_k (delta_{-exp(-kd)} + delta_{exp(-kd)}).
The outer scale is fixed at one by dilation invariance.  Thus K=1 is the
old two-atom family, while K>=6 introduces a genuine multi-scale structure
that the archived search did not cover.  The score is the exact atomic halo
formula evaluated in double precision; all results here are candidates only.

Each job writes a JSON checkpoint atomically.  Signed caches are separate from
the pre-signature historical records, and a source/dependency/parameter change
quarantines instead of silently reusing an old checkpoint.
"""
import argparse
import json
import math
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

for key in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
            "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[key] = "1"
sys.dont_write_bytecode = True

import numpy as np
from scipy.optimize import differential_evolution
from run_provenance import load_matching_or_quarantine, signature, write_record

ROOT = Path(__file__).resolve().parent
CACHE_ROOT = ROOT / "signed_cache" / "multiscale"


def merge(intervals):
    result = []
    for left, right in sorted(intervals):
        if left >= right:
            continue
        if result and left <= result[-1][1]:
            result[-1] = (result[-1][0], max(result[-1][1], right))
        else:
            result.append((left, right))
    return result


def block_halo(atoms):
    """Strict {M_c sigma>1} modulo endpoints, from all contiguous blocks."""
    atoms = sorted(atoms)
    x = [point for point, mass in atoms]
    intervals = []
    for left in range(len(atoms)):
        for right in range(left, len(atoms)):
            # Do not subtract global prefix sums: a micro-scale block can be
            # lost after being added to an unrelated macroscopic mass.
            mass = math.fsum(mass for _, mass in atoms[left:right + 1])
            a = x[right] - mass / 2
            b = x[left] + mass / 2
            if a < b:
                intervals.append((a, b))
    return merge(intervals)


def weighted_mass(intervals, tau):
    total = 0.0
    for a, b in intervals:
        if b <= 0:
            total += ((-a) ** tau - (-b) ** tau) / tau
        elif a >= 0:
            total += (b ** tau - a ** tau) / tau
        else:
            total += ((-a) ** tau + b ** tau) / tau
    return total


def score(atoms, tau):
    denominator = sum(mass * abs(point) ** (tau - 1) for point, mass in atoms)
    return weighted_mass(block_halo(atoms), tau) / denominator


def chain(vector, k):
    """Parameters: d in [0.04,1], followed by K log2 masses."""
    d = vector[0]
    masses = [2.0 ** exponent for exponent in vector[1:]]
    assert len(masses) == k
    positive = [(math.exp(-d * j), masses[j]) for j in range(k)]
    return sorted([(-x, mass) for x, mass in positive] + positive)


def job(item):
    tau, k, seed, maxiter, popsize = item
    name = f"tau{int(100*tau):02d}_K{k:02d}_seed{seed}"
    output = CACHE_ROOT / f"{name}.json"
    expected = signature(Path(__file__),
                         {"tau": tau, "K": k, "seed": seed,
                          "maxiter": maxiter, "popsize": popsize},
                         {"family": "symmetric_geometric_chain",
                          "d_bounds": [0.04, 1.0],
                          "log2_mass_bounds": [-22.0, 5.0],
                          "outer_scale": 1.0},
                         [ROOT / "run_provenance.py"])
    cached = load_matching_or_quarantine(output, expected)
    if cached is not None:
        return cached
    started = time.time()
    rng = np.random.default_rng(seed)
    bounds = [(0.04, 1.0)] + [(-22.0, 5.0)] * k
    dimension = len(bounds)
    population = np.array([[rng.uniform(a, b) for a, b in bounds]
                           for _ in range(popsize * dimension)])
    # The old symmetric pair is present as a seed; remaining scale masses are
    # tiny, so a genuine improvement cannot be due to omitted initialization.
    pair_seed = np.array([0.4] + [math.log2(4.0 / 3.0)] + [-20.0] * (k - 1))
    geometric_seed = np.array([0.25] + [math.log2(4.0 / 3.0) - 0.55*j
                                         for j in range(k)])
    population[0] = pair_seed
    population[1] = geometric_seed
    calls = 0

    def objective(vector):
        nonlocal calls
        calls += 1
        value = score(chain(vector, k), tau)
        return -value if math.isfinite(value) else 1e100

    initial = max(-objective(v) for v in population)
    result = differential_evolution(objective, bounds, init=population,
                                    seed=seed, maxiter=maxiter, popsize=popsize,
                                    tol=0.0, atol=0.0, polish=False,
                                    updating="immediate")
    atoms = chain(result.x, k)
    halo = block_halo(atoms)
    final = score(atoms, tau)
    record = {
        "status": "completed", "tau": tau, "K": k, "seed": seed,
        "maxiter": maxiter, "popsize": popsize, "evaluations": calls,
        "seconds": time.time() - started, "initial_score": initial,
        "final_score": final, "d": float(result.x[0]),
        "log2_masses": [float(y) for y in result.x[1:]],
        "atoms": [[float(x), float(m)] for x, m in atoms],
        "halo": [[float(a), float(b)] for a, b in halo],
        "message": str(result.message),
    }
    record = write_record(output, record, expected)
    return record


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--workers", type=int, default=15)
    parser.add_argument("--maxiter", type=int, default=900)
    parser.add_argument("--popsize", type=int, default=12)
    parser.add_argument("--seeds", type=int, default=2)
    parser.add_argument("--ks", type=int, nargs="+", default=[6, 10, 14])
    args = parser.parse_args()
    CACHE_ROOT.mkdir(parents=True, exist_ok=True)
    tasks = [(tau, k, 250000 + int(100*tau)*1000 + 31*k + rep,
              args.maxiter, args.popsize)
             for tau in (0.25, 0.50) for k in args.ks for rep in range(args.seeds)]
    print(json.dumps({"event": "start", "tasks": len(tasks), "workers": args.workers,
                      "maxiter": args.maxiter, "popsize": args.popsize,
                      "ks": args.ks}, sort_keys=True), flush=True)
    rows = []
    with ProcessPoolExecutor(max_workers=args.workers) as executor:
        futures = {executor.submit(job, task): task for task in tasks}
        for future in as_completed(futures):
            row = future.result()
            rows.append(row)
            print(json.dumps({"event": "done", "tau": row["tau"], "K": row["K"],
                              "seed": row["seed"], "initial": row["initial_score"],
                              "final": row["final_score"], "evals": row["evaluations"],
                              "seconds": row["seconds"]}, sort_keys=True), flush=True)
    rows.sort(key=lambda row: (row["tau"], row["K"], row["seed"]))
    (CACHE_ROOT / "summary.json").write_text(json.dumps(rows, indent=2))
    print(json.dumps({"event": "complete", "rows": len(rows)}, sort_keys=True))


if __name__ == "__main__":
    main()
