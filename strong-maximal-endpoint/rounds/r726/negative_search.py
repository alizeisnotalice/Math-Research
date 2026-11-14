#!/usr/bin/env python3
"""Negative-control extension of the r726 finite-angle proxy search."""
from __future__ import annotations

import importlib.util
import json
import math
import os
import resource
import sys
import time
import traceback
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

for _name in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS",
              "VECLIB_MAXIMUM_THREADS", "NUMEXPR_MAXIMUM_THREADS"):
    os.environ[_name] = "1"

import numpy as np
from scipy import special
from scipy.optimize import minimize

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent
MAIN_PATH = HERE / "nonlinear_search.py"
spec = importlib.util.spec_from_file_location("r726_nonlinear_definitions", MAIN_PATH)
mainmod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mainmod)

WORKERS = 12
BOUND = (-8.0, 8.0)
BOUND_TOL = 1e-7


def peak_rss_bytes():
    value = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return int(value if sys.platform == "darwin" else value * 1024)


def optimize_extended(n, eps, a, b, z, frequency, theta, oldcp, oldcm):
    mu = mainmod.MU
    starts = ((oldcp, oldcm), (-2.0, mu), (mu, -2.0), (-2.0, -2.0),
              (0.0, 0.0), (mu, mu), (4.0, 4.0))
    candidates = []
    status_records = []
    for start in starts:
        result = minimize(lambda x: -mainmod.proxy_value(x[0], x[1], theta, eps,
                                                         a, b, z, frequency),
                          x0=np.array(start, dtype=float), method="L-BFGS-B",
                          bounds=(BOUND, BOUND),
                          options={"ftol": 1e-13, "gtol": 1e-9, "maxiter": 400,
                                   "maxfun": 2000})
        cp, cm = map(float, result.x)
        mask = ((1 if cp <= -8 + BOUND_TOL else 0) |
                (2 if cp >= 8 - BOUND_TOL else 0) |
                (4 if cm <= -8 + BOUND_TOL else 0) |
                (8 if cm >= 8 - BOUND_TOL else 0))
        zero = ((1 if abs(cp) <= BOUND_TOL else 0) |
                (2 if abs(cm) <= BOUND_TOL else 0))
        score = -float(result.fun)
        candidates.append((score, cp, cm, bool(result.success), int(result.status),
                           int(result.nfev), mask, zero))
        status_records.append([int(result.success), int(result.status), int(result.nfev), mask, zero])
    finite = [i for i, x in enumerate(candidates) if math.isfinite(x[0])]
    if not finite:
        raise ArithmeticError(f"all extended proxy starts nonfinite at n={n},theta={theta}")
    winner = max(finite, key=lambda i: candidates[i][0])
    best = candidates[winner]
    return best, winner, status_records


def run_parameter(task):
    n, pattern, eps, original, r0 = task
    start = time.monotonic()
    try:
        a, b, z = original["counts_plus_minus_zero"]
        frequency = float(original["frequency_t"])
        angles = [float(x) for x in original["angle_nodes"]]
        oldcp = [float(x) for x in original["best_proxy_cplus"]]
        oldcm = [float(x) for x in original["best_proxy_cminus"]]
        oldproxy = [float(x) for x in original["best_proxy_values"]]
        newcp, newcm, newproxy, old_feasible_proxy, winners, statuses = [], [], [], [], [], []
        for theta, cp, cm in zip(angles, oldcp, oldcm):
            best, winner, status_records = optimize_extended(
                n, eps, a, b, z, frequency, theta, cp, cm)
            old_feasible = mainmod.proxy_value(cp, cm, theta, eps, a, b, z, frequency)
            if best[0] < old_feasible:
                # Keep the actual original feasible control if all local
                # extended-domain runs descend from their starts.
                cp_mask = ((1 if cp <= -8 + BOUND_TOL else 0) |
                           (2 if cp >= 8 - BOUND_TOL else 0))
                cm_mask = ((4 if cm <= -8 + BOUND_TOL else 0) |
                           (8 if cm >= 8 - BOUND_TOL else 0))
                zero_mask = ((1 if abs(cp) <= BOUND_TOL else 0) |
                             (2 if abs(cm) <= BOUND_TOL else 0))
                best = (old_feasible, cp, cm, True, -1, 0, cp_mask | cm_mask, zero_mask)
                winner = -1
            old_feasible_proxy.append(old_feasible)
            newproxy.append(best[0]); newcp.append(best[1]); newcm.append(best[2])
            winners.append(winner); statuses.append(status_records)

        proxy_deltas = [new - old for new, old in zip(newproxy, old_feasible_proxy)]
        # Reuse the exact per-parameter/per-replicate Sobol seeds. Recompute
        # actual medians because the proxy control field may have changed.
        reps = []
        for old_rep in original["replicates"]:
            power, seed = old_rep["sample_power"], old_rep["seed"]
            medians, pmf_resid, finite = mainmod.sample_medians(
                n, pattern, eps, angles, newcp, newcm, power, seed)
            logi0 = math.log(float(special.i0e(eps))) + abs(eps)
            newscore = mainmod.logmeanexp(np.asarray(medians) - mainmod.MU * n) - logi0
            even_score = mainmod.logmeanexp(np.asarray(medians)[::2] - mainmod.MU * n) - logi0
            oldscore = float(old_rep["logscore"])
            reps.append({
                "sample_power": power, "sample_size": 1 << power, "seed": seed,
                "new_actual_medians": [mainmod.compact_float(x) for x in medians],
                "old_logscore": old_rep["logscore"],
                "new_logscore": mainmod.compact_float(newscore),
                "new_logscore_even32": mainmod.compact_float(even_score),
                "new_even32_minus_full64": mainmod.compact_float(even_score-newscore),
                "new_minus_old_logscore": mainmod.compact_float(newscore-oldscore),
                "new_logscore_minus_mu_n": mainmod.compact_float(newscore),
                "new_normalized_logscore_minus_mu_n_over_logn": mainmod.compact_float(newscore/math.log(n)),
                "new_relative_logscore_minus_baseline_r0": mainmod.compact_float(newscore+mainmod.MU*n-r0),
                "new_relative_to_r0_positive": bool(newscore+mainmod.MU*n > r0),
                "max_binomial_pmf_normalization_residual": mainmod.compact_float(pmf_resid),
                "all_samples_finite": finite,
            })
        boundary_counts = {
            "best_cplus_minus8": sum(x <= -8 + BOUND_TOL for x in newcp),
            "best_cplus_zero": sum(abs(x) <= BOUND_TOL for x in newcp),
            "best_cplus_plus8": sum(x >= 8 - BOUND_TOL for x in newcp),
            "best_cminus_minus8": sum(x <= -8 + BOUND_TOL for x in newcm),
            "best_cminus_zero": sum(abs(x) <= BOUND_TOL for x in newcm),
            "best_cminus_plus8": sum(x >= 8 - BOUND_TOL for x in newcm),
        }
        winner_failures = sum(winners[i] >= 0 and not bool(statuses[i][winners[i]][0])
                              for i in range(len(winners)))
        non_success_starts = sum(not bool(s[0]) for row in statuses for s in row)
        positive = [x["sample_power"] for x in reps if x["new_relative_to_r0_positive"]]
        return {
            "n": n, "pattern": pattern, "epsilon": eps,
            "status": "ok", "elapsed_seconds": time.monotonic()-start,
            "worker_peak_rss_bytes": peak_rss_bytes(),
            "search_domain": [-8, 8], "counts_plus_minus_zero": [a,b,z],
            "frequency_t": mainmod.compact_float(frequency),
            "additional_starts": [["old_best", "old_best"], ["-2", "mu"], ["mu", "-2"],
                                  ["-2", "-2"], ["0", "0"], ["mu", "mu"], ["4", "4"]],
            "optimizer_method": "SciPy L-BFGS-B, multistart local search; not a global optimization certificate.",
            "winner_start_indices_zero_based": winners,
            "per_start_status_success_code_nfev_boundmask_zeromask": statuses,
            "winner_non_success_angle_count": winner_failures,
            "retained_original_feasible_control_count": sum(w == -1 for w in winners),
            "non_success_start_count": non_success_starts,
            "original_proxy_values": [mainmod.compact_float(x) for x in oldproxy],
            "original_proxy_recomputed_from_saved_controls": [mainmod.compact_float(x) for x in old_feasible_proxy],
            "extended_proxy_values": [mainmod.compact_float(x) for x in newproxy],
            "proxy_new_minus_original": [mainmod.compact_float(x) for x in proxy_deltas],
            "max_proxy_improvement": mainmod.compact_float(max(proxy_deltas, default=0.0)),
            "max_proxy_worsening": mainmod.compact_float(max(0.0, max((-d for d in proxy_deltas), default=0.0))),
            "extended_best_cplus": [mainmod.compact_float(x) for x in newcp],
            "extended_best_cminus": [mainmod.compact_float(x) for x in newcm],
            "boundary_hit_counts": boundary_counts,
            "replicates": reps,
            "positive_candidate_replicates": positive,
        }
    except Exception as exc:
        return {"n": n, "pattern": pattern, "epsilon": eps, "status": "error",
                "elapsed_seconds": time.monotonic()-start,
                "worker_peak_rss_bytes": peak_rss_bytes(),
                "error": f"{type(exc).__name__}: {exc}",
                "traceback": traceback.format_exc(limit=4)}


def main():
    started = time.time()
    source = json.loads((HERE / "nonlinear_search.json").read_text())
    base = {x["n"]: x for x in source["baseline_medians"]}
    tasks = []
    for original in source["results"]:
        n = original["n"]
        tasks.append((n, original["pattern"], original["epsilon"], original,
                      float(base[n]["r0"])))
    results = []
    with ProcessPoolExecutor(max_workers=WORKERS) as pool:
        futures = {pool.submit(run_parameter, task): task[:3] for task in tasks}
        for future in as_completed(futures):
            try:
                result = future.result()
            except Exception as exc:
                n, pattern, eps = futures[future]
                result = {"n":n,"pattern":pattern,"epsilon":eps,"status":"worker_error",
                          "error":f"{type(exc).__name__}: {exc}"}
            results.append(result)
            print(json.dumps({"completed":len(results),"total":len(tasks),"n":result.get("n"),
                              "pattern":result.get("pattern"),"epsilon":result.get("epsilon"),
                              "status":result.get("status"),"elapsed_seconds":result.get("elapsed_seconds")}),
                  flush=True)
            write_output(source, results, started, partial=True)
    results.sort(key=lambda x:(x.get("n",0),mainmod.PATTERNS.index(x.get("pattern",mainmod.PATTERNS[0])),
                               mainmod.EPSILONS.index(x.get("epsilon",mainmod.EPSILONS[0]))))
    positives=[{"n":x["n"],"pattern":x["pattern"],"epsilon":x["epsilon"],
                "positive_sample_powers":x.get("positive_candidate_replicates",[])}
               for x in results if x.get("positive_candidate_replicates")]
    write_output(source, results, started, partial=False, positives=positives)


def write_output(source, results, started, partial, positives=None):
    max_worker=max((x.get("worker_peak_rss_bytes",0) for x in results),default=0)
    parent=peak_rss_bytes()
    output={"status":"partial" if partial else
            ("complete" if len(results)==source["expected_parameter_count"] else "incomplete"),
            "purpose":"Negative-control extension of r726 proxy search; finite multistart and randomized-Sobol diagnostics only.",
            "source_main_batch":"rounds/r726/nonlinear_search.json is read-only and remains unchanged.",
            "search_domain":[-8,8],"theta_nodes":64,"sample_powers":[15,17],
            "optimizer_note":"Added old [0,8] winner plus (-2,mu),(mu,-2),(-2,-2),(0,0),(mu,mu),(4,4). L-BFGS-B multistart is not a global optimization certificate.",
            "sampling_note":"Reuses original parameter/replicate Sobol seeds, gamma design, and common uniforms; recomputes actual medians for extended-control fields.",
            "negative_side_results":results,"positive_candidate_summaries":positives or [],
            "status_counts":{s:sum(x.get("status")==s for x in results) for s in ("ok","error","worker_error")},
            "memory_diagnostics":{"max_measured_single_worker_peak_rss_bytes":max_worker,
                                  "parent_peak_rss_bytes":parent,
                                  "conservative_12_worker_plus_parent_bytes":12*max_worker+parent,
                                  "limit_bytes":3000000000},
            "wall_seconds":time.time()-started}
    (HERE/"negative_search.json").write_text(json.dumps(output,separators=(",",":"))+"\n")
    print(json.dumps({"status":output["status"],"completed":len(results),
                      "status_counts":output["status_counts"],
                      "positive_candidate_count":len(output["positive_candidate_summaries"]),
                      "memory":output["memory_diagnostics"],"wall_seconds":output["wall_seconds"]}),
          flush=True)


if __name__ == "__main__":
    main()
