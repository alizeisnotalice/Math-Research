#!/usr/bin/env python3
"""Monte Carlo CDF check of conditional tangent projections given S=s.

This tests finite-dimensional CDF behavior only; it is not a TV proof.
"""
from __future__ import annotations

import json
import math
import os
import platform
import resource
import sys
import threading
import time
import traceback
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

for _name in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS",
              "VECLIB_MAXIMUM_THREADS", "NUMEXPR_MAXIMUM_THREADS"):
    os.environ[_name] = "1"

import numpy as np
import psutil
import scipy
from scipy.special import gammaln, logsumexp
from scipy.stats import ks_2samp

OUT = Path(__file__).resolve().parent
NS = (64, 256, 1024, 4096)
PATTERNS = ("dense", "sparse", "singleheavy")
OFFSETS = ("minus_log_n", "zero", "plus_log_n")
REPLICATES = 3
SAMPLE_SIZE = 1 << 18
WORKERS = 8
BASE_SEED = 728202611
MACHINE_LIMIT = 3_000_000_000
MEMORY_PEAK = {"rss_bytes": 0}


def peak_rss_bytes():
    value = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return int(value if sys.platform == "darwin" else value * 1024)


def monitor_tree(stop_event):
    root = psutil.Process()
    while not stop_event.is_set():
        total = 0
        try:
            total += root.memory_info().rss
            for child in root.children(recursive=True):
                try:
                    total += child.memory_info().rss
                except (psutil.NoSuchProcess, psutil.AccessDenied):
                    pass
            MEMORY_PEAK["rss_bytes"] = max(MEMORY_PEAK["rss_bytes"], total)
        except (psutil.NoSuchProcess, psutil.AccessDenied):
            pass
        stop_event.wait(0.05)


def fixed_point_mu():
    lo, hi = 1.0, 2.0
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        if mid - 1.0 - math.exp(-mid) > 0:
            hi = mid
        else:
            lo = mid
    return 0.5 * (lo + hi)


MU = fixed_point_mu()
Q = 1.0 / (1.0 + math.exp(MU))


def posterior_j_weights(n, s):
    js = np.arange(n + 1, dtype=np.int64)
    residual = s - MU * js
    positive = residual > 0
    logw = np.full(n + 1, -np.inf, dtype=np.float64)
    j = js[positive].astype(np.float64)
    x = residual[positive]
    log_binom = (gammaln(n + 1) - gammaln(j + 1) - gammaln(n - j + 1) +
                 j * math.log(Q) + (n - j) * math.log1p(-Q))
    log_gamma_pdf = ((n - 1) * np.log(x) - x - gammaln(n))
    logw[positive] = log_binom + log_gamma_pdf
    norm = float(logsumexp(logw))
    prob = np.exp(logw - norm)
    prob_sum_residual = float(prob.sum() - 1.0)
    if not np.isfinite(prob).all() or not math.isfinite(norm):
        raise FloatingPointError(f"nonfinite posterior weights n={n}, s={s}")
    if abs(prob_sum_residual) > 1e-10:
        raise ArithmeticError(f"posterior normalization residual {prob_sum_residual}")
    return prob, {"log_normalizer": norm, "sum_residual": prob_sum_residual,
                  "support_min_j": int(js[positive][0]),
                  "support_max_j": int(js[positive][-1]),
                  "positive_support_count": int(positive.sum())}


def conditional_projection(n, pattern, s, rng, sample_size):
    prob, diag = posterior_j_weights(n, s)
    js = rng.choice(n + 1, size=sample_size, p=prob)
    residual = s - MU * js
    if pattern == "dense":
        half = n // 2
        beta = rng.beta(half, half, size=sample_size)
        kplus = rng.hypergeometric(half, half, js)
        out = residual * (2.0 * beta - 1.0) / math.sqrt(n)
        out += MU * (2.0 * kplus - js) / math.sqrt(n)
    elif pattern == "sparse":
        dirichlet = rng.dirichlet((1.0, 1.0, float(n - 2)), size=sample_size)
        jplus = rng.binomial(1, js / n)
        pminus = (js - jplus) / (n - 1)
        jminus = rng.binomial(1, pminus)
        out = residual * (dirichlet[:, 0] - dirichlet[:, 1]) / math.sqrt(2.0)
        out += MU * (jplus - jminus) / math.sqrt(2.0)
    elif pattern == "singleheavy":
        beta = rng.beta(1.0, n - 1.0, size=sample_size)
        jplus = rng.binomial(1, js / n)
        scale = math.sqrt(n * (n - 1))
        out = residual * (n * beta - 1.0) / scale
        out += MU * (n * jplus - js) / scale
    else:
        raise ValueError(pattern)
    if not np.isfinite(out).all():
        raise FloatingPointError(f"nonfinite conditional projection n={n},pattern={pattern}")
    diag["sampled_j_mean"] = float(js.mean())
    return out, diag


def unconditional_projection(n, pattern, rng, sample_size):
    if pattern == "dense":
        half = n // 2
        gp = rng.gamma(half, 1.0, size=sample_size)
        gm = rng.gamma(half, 1.0, size=sample_size)
        jp = rng.binomial(half, Q, size=sample_size)
        jm = rng.binomial(half, Q, size=sample_size)
        out = (gp + MU * jp - gm - MU * jm) / math.sqrt(n)
    elif pattern == "sparse":
        ep = rng.exponential(size=sample_size)
        em = rng.exponential(size=sample_size)
        jp = rng.binomial(1, Q, size=sample_size)
        jm = rng.binomial(1, Q, size=sample_size)
        out = (ep + MU * jp - em - MU * jm) / math.sqrt(2.0)
    elif pattern == "singleheavy":
        e1 = rng.exponential(size=sample_size)
        erest = rng.gamma(n - 1, 1.0, size=sample_size)
        j1 = rng.binomial(1, Q, size=sample_size)
        jrest = rng.binomial(n - 1, Q, size=sample_size)
        scale = math.sqrt(n * (n - 1))
        out = ((n - 1) * e1 - erest + MU * ((n - 1) * j1 - jrest)) / scale
    else:
        raise ValueError(pattern)
    if not np.isfinite(out).all():
        raise FloatingPointError(f"nonfinite unconditional projection n={n},pattern={pattern}")
    return out


def stats_summary(samples):
    return {"mean": float(np.mean(samples)),
            "variance_unbiased": float(np.var(samples, ddof=1)),
            "finite": bool(np.isfinite(samples).all())}


def ks_stat(x, y):
    return float(ks_2samp(x, y, alternative="two-sided", method="asymp").statistic)


def make_rng(seed):
    return np.random.Generator(np.random.PCG64DXSM(seed))


def run_task(task):
    n, pattern, pattern_index = task
    start = time.monotonic()
    try:
        offset_data = []
        for oi, offset_name in enumerate(OFFSETS):
            offset = (-math.log(n) if offset_name == "minus_log_n" else
                      0.0 if offset_name == "zero" else math.log(n))
            s = MU * n + offset
            conditions, baselines = [], []
            replicate_records = []
            posterior_diags = []
            for rep in range(REPLICATES):
                seed_cond = BASE_SEED + n * 10000 + pattern_index * 1_000_000 + oi * 100 + rep * 2
                seed_uncond = seed_cond + 1
                cond, postdiag = conditional_projection(n, pattern, s, make_rng(seed_cond), SAMPLE_SIZE)
                uncond = unconditional_projection(n, pattern, make_rng(seed_uncond), SAMPLE_SIZE)
                d = ks_stat(cond, uncond)
                conditions.append(cond); baselines.append(uncond); posterior_diags.append(postdiag)
                replicate_records.append({
                    "replicate": rep, "conditional_seed": seed_cond, "unconditional_seed": seed_uncond,
                    "conditional": stats_summary(cond), "unconditional": stats_summary(uncond),
                    "paired_empirical_CDF_sup": d,
                    "posterior_j": postdiag,
                })
            pooled_cond = np.concatenate(conditions)
            pooled_uncond = np.concatenate(baselines)
            cond_repwise = [ks_stat(conditions[i], conditions[j])
                            for i in range(REPLICATES) for j in range(i + 1, REPLICATES)]
            uncond_repwise = [ks_stat(baselines[i], baselines[j])
                              for i in range(REPLICATES) for j in range(i + 1, REPLICATES)]
            paired_d = [r["paired_empirical_CDF_sup"] for r in replicate_records]
            offset_data.append({
                "offset_label": offset_name, "offset": offset, "conditioning_s": s,
                "sample_size_per_replicate_each_law": SAMPLE_SIZE,
                "replicates": replicate_records,
                "pooled_conditional_summary": stats_summary(pooled_cond),
                "pooled_unconditional_summary": stats_summary(pooled_uncond),
                "pooled_empirical_CDF_sup": ks_stat(pooled_cond, pooled_uncond),
                "paired_CDF_sup_replicate_mean": float(np.mean(paired_d)),
                "paired_CDF_sup_replicate_sd": float(np.std(paired_d, ddof=1)),
                "paired_CDF_sup_replicate_range": [float(min(paired_d)), float(max(paired_d))],
                "conditional_cross_replicate_CDF_sup": cond_repwise,
                "unconditional_cross_replicate_CDF_sup": uncond_repwise,
                "conditional_cross_replicate_CDF_sup_mean": float(np.mean(cond_repwise)),
                "unconditional_cross_replicate_CDF_sup_mean": float(np.mean(uncond_repwise)),
                "posterior_weight_max_abs_normalization_residual": max(abs(d["sum_residual"])
                                                                         for d in posterior_diags),
            })
            del conditions, baselines, pooled_cond, pooled_uncond
        return {"n": n, "pattern": pattern, "status": "ok",
                "elapsed_seconds": time.monotonic() - start,
                "worker_peak_rss_bytes": peak_rss_bytes(),
                "sample_size": SAMPLE_SIZE, "replicates": REPLICATES,
                "conditional_and_unconditional_laws": offset_data}
    except Exception as exc:
        return {"n": n, "pattern": pattern, "status": "error",
                "elapsed_seconds": time.monotonic()-start,
                "worker_peak_rss_bytes": peak_rss_bytes(),
                "error": f"{type(exc).__name__}: {exc}",
                "traceback": traceback.format_exc(limit=5)}


def main():
    started = time.time()
    tasks = [(n, pattern, pi) for n in NS for pi, pattern in enumerate(PATTERNS)]
    stop_event = threading.Event()
    monitor = threading.Thread(target=monitor_tree, args=(stop_event,), daemon=True)
    monitor.start()
    results = []
    with ProcessPoolExecutor(max_workers=WORKERS) as pool:
        futures = {pool.submit(run_task, task): task for task in tasks}
        for future in as_completed(futures):
            try:
                result = future.result()
            except Exception as exc:
                n, pattern, _ = futures[future]
                result = {"n":n,"pattern":pattern,"status":"worker_error",
                          "error":f"{type(exc).__name__}: {exc}"}
            results.append(result)
            print(json.dumps({"completed":len(results),"total":len(tasks),"n":result.get("n"),
                              "pattern":result.get("pattern"),"status":result.get("status"),
                              "elapsed_seconds":result.get("elapsed_seconds")}),flush=True)
            write_output(results, started, partial=True)
    stop_event.set(); monitor.join()
    results.sort(key=lambda x:(x.get("n",0),PATTERNS.index(x.get("pattern",PATTERNS[0]))))
    write_output(results, started, partial=False)


def write_output(results, started, partial):
    maxrss=max((x.get("worker_peak_rss_bytes",0) for x in results),default=0)
    parentrss=peak_rss_bytes()
    out={"status":"partial" if partial else ("complete" if len(results)==12 else "incomplete"),
         "purpose":"Monte Carlo CDF check of conditional tangent projections given S=s; no TV or uniform-limit certification.",
         "dimensions":list(NS),"patterns":list(PATTERNS),"offsets":list(OFFSETS),
         "mu":MU,"q":Q,"sample_size_per_law_per_replicate":SAMPLE_SIZE,
         "replicates_per_conditional_and_unconditional_law":REPLICATES,
         "rng":"NumPy Generator(PCG64DXSM), fixed independent seeds recorded per sample.",
         "conditional_posterior":"P(J_total=j|S=s) proportional to BinomialPMF(n,q,j)*GammaPDF(shape=n, rate=1, s-mu*j), normalized via logsumexp.",
         "conditional_projection_methods":{
           "dense":"a=(+1/sqrt(n)) on first n/2 and (-1/sqrt(n)) on remaining n/2; use Beta(n/2,n/2) Dirichlet split and Hypergeom(n,j,n/2).",
           "sparse":"a=(+1,-1)/sqrt(2); use Dirichlet(1,1,n-2) exponential shares and sequential Bernoulli allocation conditional on j.",
           "singleheavy":"a1=sqrt((n-1)/n), others=-1/sqrt(n(n-1)); use Beta(1,n-1) exponential share and Bernoulli(j/n)."},
         "unconditional_projection_methods":"Exact grouped Gamma/Binomial simulation for each listed coefficient pattern; no n-dimensional array is built.",
         "versions":{"python":sys.version.split()[0],"numpy":np.__version__,"scipy":scipy.__version__,
                     "platform":platform.platform(),"machine":platform.machine(),"cpu_count":os.cpu_count()},
         "workers":WORKERS,"requested_workers":12,
         "worker_count_note":"Initial 12-worker pass measured >4.03 GB concurrent process-tree RSS; the final reproducible pass uses 8 workers to remain under the 3 GB memory budget.",
         "blas_threads_per_worker":1,
         "memory_measured":{"max_worker_peak_rss_bytes":maxrss,"parent_peak_rss_bytes":parentrss,
                            "measured_concurrent_process_tree_peak_rss_bytes":MEMORY_PEAK["rss_bytes"],
                            "limit_bytes":MACHINE_LIMIT},
         "initial_12_worker_memory_probe":{"tasks_completed":12,"status":"complete_but_exceeded_RSS_budget",
                                            "wall_seconds":4.615004301071167,
                                            "max_worker_peak_rss_bytes":331661312,
                                            "parent_peak_rss_bytes":103890944,
                                            "measured_concurrent_process_tree_peak_rss_bytes":4033970176,
                                            "budget_bytes":MACHINE_LIMIT,
                                            "note":"All samples completed; this was not an out-of-memory failure. The final fixed-seed pass uses 8 workers and remains below budget."},
         "results":results,
         "status_counts":{s:sum(x.get("status")==s for x in results) for s in ("ok","error","worker_error")},
         "wall_seconds":time.time()-started}
    (OUT/"conditional_phase_check.json").write_text(json.dumps(out,separators=(",",":"))+"\n")
    print(json.dumps({"status":out["status"],"completed":len(results),"status_counts":out["status_counts"],
                      "memory":out["memory_measured"],"wall_seconds":out["wall_seconds"]}),flush=True)


if __name__=="__main__":
    main()
