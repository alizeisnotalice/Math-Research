#!/usr/bin/env python3
"""Independent high-precision computation of the micro-energy slope."""
from __future__ import annotations

import importlib.util
import json
import os
import time
import traceback
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

for _name in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS",
              "VECLIB_MAXIMUM_THREADS", "NUMEXPR_MAXIMUM_THREADS"):
    os.environ[_name] = "1"

import mpmath as mp

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent
R722_FILE = ROOT / "rounds/r722/curvature.py"
_spec = importlib.util.spec_from_file_location("r722_curvature_definitions", R722_FILE)
_r722 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_r722)

NS = (4, 8, 16, 32, 64, 128, 256, 512, 1024)
PRECISIONS = (50, 80)
WORKERS = 12
MAXTERMS = 20000


def fmt(x, dps):
    return mp.nstr(x, n=dps)


def conditional_second_moment(N, x, c, p):
    """Conditional E[Z_1^2 | sum_{i<=N} Z_i=x] by the exact finite mixture."""
    if x <= 0:
        raise ValueError("conditional total must be positive")
    density = _r722.mixture_pdf(N, x, c, p)
    if density <= 0:
        raise ArithmeticError("nonpositive conditional mixture density")
    terms = []
    for j in range(N + 1):
        t = x - c * j
        if t <= 0:
            continue
        joint = _r722.binom_weight(N, j, p) * _r722.gamma_pdf(N, t)
        second_j = (2 * t**2 / (N * (N + 1)) +
                    2 * c * t * j / N**2 + c**2 * j / N)
        terms.append(joint * second_j)
    second = mp.fsum(terms) / density
    mean = x / N
    variance_scaled = mp.mpf(N) / (N - 1) * (second - mean**2)
    return second, variance_scaled, density


def solve_baseline(task):
    n, dps = task
    start = time.monotonic()
    try:
        b = _r722.baseline(n, dps)
        rec = {k: (fmt(v, dps) if isinstance(v, mp.mpf) else v)
               for k, v in b.items()}
        cache = {key: rec[key] for key in
                 ("r", "c", "p", "cdf_n_residual", "cdf_nminus1_shifted_residual")}
        return {"n": n, "dps": dps, "status": "ok", "baseline": cache,
                "record": {"status": "ok", **rec},
                "elapsed_seconds": time.monotonic() - start}
    except Exception as exc:
        return {"n": n, "dps": dps, "status": "error",
                "elapsed_seconds": time.monotonic() - start,
                "error": f"{type(exc).__name__}: {exc}",
                "traceback": traceback.format_exc(limit=3)}


def compute(task):
    n, dps, b = task
    start = time.monotonic()
    mp.mp.dps = dps
    mp.mp.maxterms = MAXTERMS
    try:
        r = mp.mpf(b["r"])
        c = mp.mpf(b["c"])
        p = mp.mpf(b["p"])
        N = n - 1
        fn = _r722.mixture_pdf(n, r, c, p)
        fN_r = _r722.mixture_pdf(N, r, c, p)
        fN_shift = _r722.mixture_pdf(N, r - c, c, p)
        cdfNm2 = _r722.mixture_cdf(n - 2, r - 2 * c, c, p)
        q = p
        h = q * fN_shift / fn
        offb = p**2 * (cdfNm2 - mp.mpf("0.5")) / fn
        second_r, v_r, cond_density_r = conditional_second_moment(N, r, c, p)
        second_shift, v_shift, cond_density_shift = conditional_second_moment(N, r - c, c, p)
        a_n = ((1 - p) * fN_r / fn * v_r + h * v_shift)
        eta = (n * c - r) / (n - 1)
        denominator = h - offb
        kappa = h**2 * eta**2 / denominator
        slope = (kappa - a_n) / 4
        residuals = {
            "CDF_n_minus_half": b["cdf_n_residual"],
            "CDF_nminus1_shifted_minus_half": b["cdf_nminus1_shifted_residual"],
            "recomputed_CDF_n_minus_half": fmt(_r722.mixture_cdf(n, r, c, p) - mp.mpf("0.5"), dps),
            "recomputed_CDF_nminus1_shifted_minus_half": fmt(_r722.mixture_cdf(N, r - c, c, p) - mp.mpf("0.5"), dps),
        }
        signs = {"fn": fn > 0, "fN_r": fN_r > 0, "fN_shift": fN_shift > 0,
                 "h": h > 0, "h_minus_offb": denominator > 0,
                 "vtilde_r": v_r > 0, "vtilde_shift": v_shift > 0,
                 "a_n": a_n > 0, "kappa": kappa >= 0}
        return {"n": n, "precision_dps": dps, "status": "ok",
                "elapsed_seconds": time.monotonic() - start,
                "baseline": {"r": fmt(r, dps), "c": fmt(c, dps), "p": fmt(p, dps),
                             "fn": fmt(fn, dps), "fN_r": fmt(fN_r, dps),
                             "fN_r_minus_c": fmt(fN_shift, dps),
                             "CDF_nminus2_rminus2c": fmt(cdfNm2, dps)},
                "quantities": {"q_equals_p": fmt(q, dps), "h": fmt(h, dps),
                               "offb": fmt(offb, dps), "h_minus_offb": fmt(denominator, dps),
                               "eta": fmt(eta, dps),
                               "conditional_second_moment_N_at_r": fmt(second_r, dps),
                               "conditional_second_moment_N_at_r_minus_c": fmt(second_shift, dps),
                               "vtilde_N_at_r": fmt(v_r, dps),
                               "vtilde_N_at_r_minus_c": fmt(v_shift, dps),
                               "a_n": fmt(a_n, dps), "kappa_n": fmt(kappa, dps),
                               "micro_energy_slope": fmt(slope, dps)},
                "conditional_density_checks": {"density_N_r": fmt(cond_density_r, dps),
                                               "density_N_r_minus_c": fmt(cond_density_shift, dps)},
                "positivity_checks": signs, "all_required_positive": all(signs.values()),
                "micro_energy_slope_sign": ("positive" if slope > 0 else
                                             "negative" if slope < 0 else "zero"),
                "baseline_residuals": residuals}
    except Exception as exc:
        return {"n": n, "precision_dps": dps, "status": "error",
                "elapsed_seconds": time.monotonic() - start,
                "error": f"{type(exc).__name__}: {exc}",
                "traceback": traceback.format_exc(limit=3)}


def main():
    start_wall = time.time()
    start_mono = time.monotonic()
    baselines = {}
    baseline_records = []
    results = []
    timings = {}
    with ProcessPoolExecutor(max_workers=WORKERS) as pool:
        # Every baseline (n,dps) is its own worker task; this stage was the
        # dominant cost in the initial serial-parent implementation.
        stage = time.monotonic()
        bfutures = {pool.submit(solve_baseline, (n, dps)): (n, dps)
                    for n in NS for dps in PRECISIONS}
        for fut in as_completed(bfutures):
            n, dps = bfutures[fut]
            item = {"n": n, "dps": dps, "status": "worker_error"}
            try:
                item = fut.result()
                if item["status"] == "ok":
                    baselines[(n, dps)] = item["baseline"]
                    baseline_records.append(item["record"] | {
                        "worker_elapsed_seconds": item["elapsed_seconds"]})
                else:
                    baseline_records.append({k: v for k, v in item.items()
                                             if k not in ("baseline", "record")})
            except Exception as exc:
                baseline_records.append({"n": n, "dps": dps, "status": "error",
                                         "error": f"{type(exc).__name__}: {exc}",
                                         "traceback": traceback.format_exc(limit=3)})
            print(json.dumps({"stage": "baseline", "completed": len(baseline_records),
                              "total": len(NS) * len(PRECISIONS), "n": n,
                              "dps": dps, "status": item.get("status", "worker_error")}), flush=True)
        timings["parallel_baseline_wall_seconds"] = time.monotonic() - stage
        timings["parallel_baseline_max_worker_seconds"] = max(
            (x.get("worker_elapsed_seconds", x.get("elapsed_seconds", 0))
             for x in baseline_records), default=0)
        stage = time.monotonic()
        futures = {}
        for n in NS:
            for dps in PRECISIONS:
                b = baselines.get((n, dps))
                if b is None:
                    results.append({"n": n, "precision_dps": dps, "status": "baseline_error"})
                else:
                    futures[pool.submit(compute, (n, dps, b))] = (n, dps)
        for fut in as_completed(futures):
            try:
                results.append(fut.result())
            except Exception as exc:
                n, dps = futures[fut]
                results.append({"n": n, "precision_dps": dps, "status": "worker_error",
                                "error": f"{type(exc).__name__}: {exc}"})
            write_json(baseline_records, results, start_wall, partial=True,
                       stage_timings=timings)
        timings["parallel_moment_wall_seconds"] = time.monotonic() - stage
    results.sort(key=lambda x: (x.get("n", 0), x.get("precision_dps", 0)))
    bykey = {(x.get("n"), x.get("precision_dps")): x for x in results}
    diffs = []
    for n in NS:
        lo, hi = bykey.get((n, 50)), bykey.get((n, 80))
        if lo and hi and lo.get("status") == hi.get("status") == "ok":
            for field in ("h", "offb", "eta", "vtilde_N_at_r", "vtilde_N_at_r_minus_c",
                          "a_n", "kappa_n", "micro_energy_slope"):
                x = mp.mpf(lo["quantities"][field]); y = mp.mpf(hi["quantities"][field])
                diffs.append({"n": n, "field": field,
                              "abs_difference_50_vs_80": mp.nstr(abs(x - y), 20)})
    timings["total_wall_seconds"] = time.monotonic() - start_mono
    write_json(baseline_records, results, start_wall, partial=False, differences=diffs,
               stage_timings=timings)


def write_json(baselines, results, start, partial, differences=None, stage_timings=None):
    output = {"status": "partial" if partial else
              ("complete" if len(results) == len(NS) * len(PRECISIONS) else "incomplete"),
              "purpose": "Independent finite-mixture calculation of the micro-energy slope; numerical verification only.",
              "dimensions": list(NS), "precisions_dps": list(PRECISIONS),
              "processes": WORKERS, "blas_threads_per_worker": 1,
              "cdf_pdf_implementation": "Imported definitions only from rounds/r722/curvature.py; no r722 main execution.",
              "baseline_cache": "Each (n,dps) baseline was solved in its own worker task in the same 12-process pool, then cached and passed to its moment task.",
              "conditional_moment_formula": "For t=x-cj, E[Z1^2|j,total=x]=2t^2/[N(N+1)]+2ctj/N^2+c^2 j/N; posterior j weights proportional to Bin(N,j)p^j(1-p)^(N-j) GammaPDF_N(t).",
              "slope_formula": "(kappa_n-a_n)/4 with a_n, eta, kappa_n as provided in the task.",
              "numerical_scope": "mpmath finite sums and root solves; 50/80 dps differences are diagnostics, not interval certification.",
              "baseline_solutions": baselines, "results": results,
              "stage_timings": stage_timings or {},
              "precision_differences": differences or [],
              "status_counts": {s: sum(1 for x in results if x.get("status") == s)
                                for s in ("ok", "error", "worker_error", "baseline_error")},
              "wall_seconds": time.time() - start}
    (HERE / "micro.json").write_text(json.dumps(output, indent=2) + "\n")
    print(json.dumps({"status": output["status"], "completed": len(results),
                      "status_counts": output["status_counts"],
                      "wall_seconds": output["wall_seconds"]}, indent=2), flush=True)


if __name__ == "__main__":
    main()
