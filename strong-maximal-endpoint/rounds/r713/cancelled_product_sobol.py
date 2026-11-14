#!/usr/bin/env python3
"""Scrambled-Sobol numerical probe for the signed product cancellation model.

This is numerical exploration only (no confidence or rigorous error claim).
Each job draws n independent uniforms V_i on [-r,1] and an independent
Exp(1) variable Z, evaluating log|P| with stable log-domain formulas.
"""
from __future__ import annotations

import json
import math
import os
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

for _name in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS",
              "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[_name] = "1"

import numpy as np
from scipy.stats import qmc

HERE = Path(__file__).resolve().parent
NS = (4, 8, 16, 32, 64, 128, 256, 512)
BETAS = (-4, -2, 0, 2, 4)
SIGNS = (-1, 1)
LEVELS = ((15, 71315001), (17, 71317003))
CHUNK = 4096


def solve_mu():
    # Monotone fixed-point iteration for the unique root mu=1+exp(-mu).
    mu = 1.3
    for _ in range(1000):
        nxt = 1.0 + math.exp(-mu)
        if abs(nxt - mu) <= 1e-15:
            return nxt
        mu = nxt
    raise RuntimeError("mu iteration did not converge")


MU = solve_mu()
R = MU - 1.0


def one_job(job):
    n, beta, theta_sign, power, seed = job
    count = 1 << power
    logtheta = -(MU + math.log(MU) - 1.0) * n + beta * math.sqrt(n)
    engine = qmc.Sobol(d=n + 1, scramble=True, seed=seed)
    y_chunks, yz_chunks = [], []
    tail_counts = np.zeros(3, dtype=np.int64)
    badfinite = 0
    completed = 0
    while completed < count:
        take = min(CHUNK, count - completed)
        u = engine.random(take)
        v = -R + MU * u[:, :n]
        abs_v = np.abs(v)
        s = -np.log(abs_v).sum(axis=1)
        b = np.log1p(-v).sum(axis=1)
        t = logtheta + b + s

        # sign((-1)^(n-1) theta prod V_i)
        negative_parity = (np.count_nonzero(v < 0.0, axis=1) & 1)
        signprod = np.where(negative_parity == 1, -1.0, 1.0)
        sgn = ((-1.0) ** (n - 1)) * theta_sign * signprod

        logfactor = np.empty_like(t)
        plus = sgn > 0
        logfactor[plus] = np.logaddexp(0.0, t[plus])
        minus = ~plus
        if np.any(minus):
            tm = t[minus]
            # log|1-exp(T)|, stable both sides of T=0.
            with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
                logfactor[minus] = np.maximum(0.0, tm) + np.log(-np.expm1(-np.abs(tm)))

        # Y=-log|P|=S-log|1+sgn exp(T)|. Exact cancellation is a nonfinite
        # sample and is counted explicitly; it is not silently discarded.
        y = s - logfactor
        z = -np.log1p(-u[:, n])
        yz = y + z
        correction = y - s
        badfinite += int(np.count_nonzero(~(np.isfinite(y) & np.isfinite(yz) &
                                            np.isfinite(s) & np.isfinite(t))))
        for j, threshold in enumerate((1.0, 2.0, 4.0)):
            tail_counts[j] += np.count_nonzero(correction > threshold)
        y_chunks.append(y)
        yz_chunks.append(yz)
        completed += take

    ys = np.concatenate(y_chunks)
    yzs = np.concatenate(yz_chunks)
    tail = tail_counts.astype(np.float64) / count
    return {
        "n": n, "beta": beta, "theta_sign": theta_sign,
        "log_abs_theta": logtheta, "sobol_power": power, "samples": count,
        "seed": seed, "chunk_size": CHUNK,
        "median_Y_plus_Z_minus_mu_n": float(np.median(yzs) - MU * n),
        "mean_Y": float(np.mean(ys)),
        "mean_Y_minus_mu_n": float(np.mean(ys) - MU * n),
        "positive_correction_tail_times_sqrt_n": {
            str(th): float(prob * math.sqrt(n))
            for th, prob in zip((1, 2, 4), tail)
        },
        "positive_correction_tail_probability": {
            str(th): float(prob) for th, prob in zip((1, 2, 4), tail)
        },
        "badfinite_count": badfinite,
        "badfinite_rate": badfinite / count,
        "numerical_status": "ok" if badfinite == 0 else "badfinite_samples_observed",
        "interpretation": "numerical estimate only; scrambled Sobol replicates are not confidence or certification bounds",
    }


def main():
    jobs = []
    for n in NS:
        for beta in BETAS:
            for theta_sign in SIGNS:
                for power, seed_base in LEVELS:
                    seed = seed_base + n * 1000 + (beta + 4) * 31 + (theta_sign + 1) * 7
                    jobs.append((n, beta, theta_sign, power, seed))
    start = time.time()
    raw = []
    with ProcessPoolExecutor(max_workers=12) as pool:
        futures = [pool.submit(one_job, job) for job in jobs]
        for future in as_completed(futures):
            raw.append(future.result())
    raw.sort(key=lambda x: (x["n"], x["beta"], x["theta_sign"], x["sobol_power"]))

    by_key = {}
    for result in raw:
        key = (result["n"], result["beta"], result["theta_sign"])
        by_key.setdefault(key, {})[result["sobol_power"]] = result
    comparisons = []
    for key in sorted(by_key):
        low = by_key[key][15]
        high = by_key[key][17]
        comparisons.append({
            "n": key[0], "beta": key[1], "theta_sign": key[2],
            "seed_2pow15": low["seed"], "seed_2pow17": high["seed"],
            "high_minus_low_median_Y_plus_Z_minus_mu_n":
                high["median_Y_plus_Z_minus_mu_n"] - low["median_Y_plus_Z_minus_mu_n"],
            "high_minus_low_mean_Y": high["mean_Y"] - low["mean_Y"],
            "high_minus_low_tail_times_sqrt_n": {
                th: high["positive_correction_tail_times_sqrt_n"][th]
                     - low["positive_correction_tail_times_sqrt_n"][th]
                for th in ("1", "2", "4")
            },
            "badfinite_2pow15": low["badfinite_count"],
            "badfinite_2pow17": high["badfinite_count"],
            "difference_note": "two independently scrambled deterministic-seed estimates; difference is diagnostic, not an error bar",
        })

    output = {
        "status": "complete" if len(raw) == len(jobs) else "incomplete",
        "model": "P=prod(V_i)+(-1)^(n-1)*theta*prod(1-V_i), V_i iid Unif[-r,1], r=mu-1",
        "mu_equation": "mu=1+exp(-mu)",
        "mu_numerical": MU,
        "r_numerical": R,
        "Y_definition": "Y=-log(abs(P))",
        "stable_log_formula": "S=-sum(log|V_i|), B=sum(log(1-V_i)), T=log|theta|+B+S; log|P|=-S+log|1+sgn*exp(T)|",
        "stable_branches": {
            "plus": "logaddexp(0,T)",
            "minus": "max(0,T)+log(-expm1(-abs(T)))",
        },
        "theta_log_scale": "-(mu+log(mu)-1)*n + beta*sqrt(n)",
        "n_values": list(NS), "beta_values": list(BETAS), "theta_signs": list(SIGNS),
        "sobol_levels": [{"power": power, "samples": 1 << power, "seed_base": seed}
                         for power, seed in LEVELS],
        "scramble": "SciPy Sobol scramble=True; one deterministic seed per case and precision level",
        "replicate_limit": "Two independent scrambles are compared for numerical stability only; no confidence interval or rigorous error certification is claimed.",
        "chunk_size": CHUNK,
        "workers": 12,
        "blas_threads_per_worker": 1,
        "case_count": len(by_key),
        "sampling_run_count": len(raw),
        "badfinite_total": sum(row["badfinite_count"] for row in raw),
        "wall_seconds": time.time() - start,
        "precision_comparisons": comparisons,
        "runs": raw,
    }
    (HERE / "cancelled_product_sobol.json").write_text(json.dumps(output, indent=2) + "\n")
    print(json.dumps({key: output[key] for key in (
        "status", "case_count", "sampling_run_count", "badfinite_total",
        "wall_seconds", "workers", "blas_threads_per_worker")}, indent=2))


if __name__ == "__main__":
    main()
