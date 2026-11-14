#!/usr/bin/env python3
"""High-dimensional independent-scramble stability check for r726 witnesses."""
from __future__ import annotations

import importlib.util
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
from scipy import special

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent
MAIN_FILE = HERE / "nonlinear_search.py"
spec = importlib.util.spec_from_file_location("r726_sampling_definitions", MAIN_FILE)
mainmod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mainmod)

NS = (1024, 8192)
POWERS = (19, 20)
REPLICATES = 3
BASE_SEED = 726202900
EPSILON = 1.0
PATTERN = "dense_energy"
WORKERS_IF_SMALL = 12
WORKERS_IF_LARGE = 8
RSS_THRESHOLD_BYTES = 220_000_000
MEMORY_LIMIT_BYTES = 3_000_000_000
MEMORY_OBS = {"process_tree_peak_rss_bytes": 0}


def peak_rss_bytes():
    value = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return int(value if sys.platform == "darwin" else value * 1024)


def monitor_process_tree(stop_event):
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
            MEMORY_OBS["process_tree_peak_rss_bytes"] = max(
                MEMORY_OBS["process_tree_peak_rss_bytes"], total)
        except (psutil.NoSuchProcess, psutil.AccessDenied):
            pass
        stop_event.wait(0.05)


def run_sample(task):
    n, power, replicate, seed, cps, cms, r0 = task
    start = time.monotonic()
    try:
        angles = np.arange(64, dtype=float) * (2 * math.pi / 64)
        medians, pmf_resid, finite = mainmod.sample_medians(
            n, PATTERN, EPSILON, angles, cps, cms, power, seed)
        mu = mainmod.MU
        logi0 = math.log(float(special.i0e(EPSILON))) + EPSILON
        score = mainmod.logmeanexp(np.asarray(medians) - mu * n) - logi0
        score32 = mainmod.logmeanexp(np.asarray(medians)[::2] - mu * n) - logi0
        return {"n": n, "pattern": PATTERN, "epsilon": EPSILON,
                "sample_power": power, "sample_size": 1 << power,
                "replicate": replicate, "seed": seed,
                "status": "ok", "elapsed_seconds": time.monotonic()-start,
                "worker_peak_rss_bytes": peak_rss_bytes(),
                "actual_medians_64": [mainmod.compact_float(x) for x in medians],
                "logscore_minus_mu_n": mainmod.compact_float(score),
                "normalized_logscore_minus_mu_n_over_logn": mainmod.compact_float(score / math.log(n)),
                "relative_logscore_minus_baseline_r0": mainmod.compact_float(score + mu*n - r0),
                "relative_to_r0_positive": bool(score + mu*n > r0),
                "even32_logscore_minus_mu_n": mainmod.compact_float(score32),
                "even32_minus_full64": mainmod.compact_float(score32-score),
                "max_binomial_pmf_normalization_residual": mainmod.compact_float(pmf_resid),
                "all_samples_finite": bool(finite)}
    except Exception as exc:
        return {"n": n, "pattern": PATTERN, "epsilon": EPSILON,
                "sample_power": power, "sample_size": 1 << power,
                "replicate": replicate, "seed": seed,
                "status": "error", "elapsed_seconds": time.monotonic()-start,
                "worker_peak_rss_bytes": peak_rss_bytes(),
                "error": f"{type(exc).__name__}: {exc}",
                "traceback": traceback.format_exc(limit=4)}


def main():
    start = time.time()
    neg = json.loads((HERE / "negative_search.json").read_text())
    parent = json.loads((HERE / "nonlinear_search.json").read_text())
    negative_by_key = {(x["n"], x["pattern"], x["epsilon"]): x
                       for x in neg["negative_side_results"]}
    baseline_by_n = {x["n"]: x for x in parent["baseline_medians"]}
    controls = {}
    tasks = []
    for n in NS:
        record = negative_by_key[(n, PATTERN, EPSILON)]
        controls[n] = {"cplus": record["extended_best_cplus"],
                       "cminus": record["extended_best_cminus"]}
        r0 = float(baseline_by_n[n]["r0"])
        for power in POWERS:
            for rep in range(REPLICATES):
                seed = BASE_SEED + n * 100 + power * 10 + rep
                tasks.append((n, power, rep, seed,
                              [float(x) for x in controls[n]["cplus"]],
                              [float(x) for x in controls[n]["cminus"]], r0))

    # A transparent array-size estimate guides pool width before launching the
    # full run. At 2^20, Sobol's 6-column matrix alone is 48 MiB; gamma,
    # binomial and median work arrays add roughly another 60-100 MiB. Combined
    # with the measured ~139 MiB worker RSS at 2^17 from the preceding batch,
    # we conservatively expect >220 MB per worker and select eight workers.
    estimated_worker_bytes = 240_000_000
    workers = WORKERS_IF_LARGE if estimated_worker_bytes > RSS_THRESHOLD_BYTES else WORKERS_IF_SMALL
    results = []
    stop_event = threading.Event()
    monitor = threading.Thread(target=monitor_process_tree, args=(stop_event,), daemon=True)
    monitor.start()
    with ProcessPoolExecutor(max_workers=workers) as pool:
        futures = {pool.submit(run_sample, task): task[:4] for task in tasks}
        for future in as_completed(futures):
            try:
                result = future.result()
            except Exception as exc:
                n, power, rep, seed = futures[future]
                result = {"n": n, "sample_power": power, "replicate": rep,
                          "seed": seed, "status": "worker_error",
                          "error": f"{type(exc).__name__}: {exc}"}
            results.append(result)
            print(json.dumps({"completed":len(results),"total":len(tasks),
                              "n":result.get("n"),"power":result.get("sample_power"),
                              "replicate":result.get("replicate"),"status":result.get("status"),
                              "elapsed_seconds":result.get("elapsed_seconds")}), flush=True)
            write_json(results, controls, workers, estimated_worker_bytes, start, partial=True)
    stop_event.set()
    monitor.join()
    results.sort(key=lambda x:(x.get("n",0),x.get("sample_power",0),x.get("replicate",0)))
    write_json(results, controls, workers, estimated_worker_bytes, start, partial=False)


def write_json(results, controls, workers, estimated_worker_bytes, start, partial):
    status_counts={s:sum(x.get("status")==s for x in results) for s in ("ok","error","worker_error")}
    maxrss=max((x.get("worker_peak_rss_bytes",0) for x in results),default=0)
    parentrss=peak_rss_bytes()
    aggregates=[]
    for n in NS:
        for power in POWERS:
            group=[x for x in results if x.get("n")==n and x.get("sample_power")==power and x.get("status")=="ok"]
            if group:
                vals=[float(x["relative_logscore_minus_baseline_r0"]) for x in group]
                aggregates.append({"n":n,"sample_power":power,"sample_size":1<<power,
                                   "replicate_count":len(group),
                                   "relative_logscore_seed_values":vals,
                                   "min":min(vals),"max":max(vals),"range":max(vals)-min(vals),
                                   "mean":sum(vals)/len(vals)})
    power_differences=[]
    for n in NS:
        small=[x for x in results if x.get("n")==n and x.get("sample_power")==19 and x.get("status")=="ok"]
        large=[x for x in results if x.get("n")==n and x.get("sample_power")==20 and x.get("status")=="ok"]
        if small and large:
            sv=[float(x["relative_logscore_minus_baseline_r0"]) for x in small]
            lv=[float(x["relative_logscore_minus_baseline_r0"]) for x in large]
            power_differences.append({"n":n,"mean_2pow20_minus_mean_2pow19":sum(lv)/len(lv)-sum(sv)/len(sv),
                                      "2pow19_seed_range":[min(sv),max(sv)],
                                      "2pow20_seed_range":[min(lv),max(lv)],
                                      "ranges_overlap":not (max(sv)<min(lv) or max(lv)<min(sv))})
    out={"status":"partial" if partial else
         ("complete" if len(results)==12 else "incomplete"),
         "purpose":"Independent scrambled-Sobol stability check using the negative-search dense-energy epsilon=1 control field; numerical only.",
         "source_controls":"rounds/r726/negative_search.json, n=1024/8192 dense_energy epsilon=1; main and negative-search files are read-only.",
         "dimensions":list(NS),"pattern":PATTERN,"epsilon":EPSILON,
         "frequencies_t":{str(n):mainmod.compact_float(1/math.sqrt(n)) for n in NS},
         "sample_powers":list(POWERS),"replicates_per_power":REPLICATES,
         "sample_sizes_per_parameter":[1<<p for p in POWERS],
         "base_seed":BASE_SEED,
         "seed_formula":"base+n*100+power*10+replicate_index (0,1,2); new independent scramble for every task.",
         "sobol":"6-dimensional scrambled Sobol. Gamma inverse and binomial inverse transforms follow the grouped model; raw sample arrays are not saved.",
         "theta_nodes":64,"even32_quadrature":"Uses even-indexed nodes from the same 64 median values; difference to full 64 is recorded.",
         "versions":{"python":sys.version.split()[0],"numpy":np.__version__,"scipy":scipy.__version__,
                     "platform":platform.platform(),"machine":platform.machine(),"cpu_count":os.cpu_count()},
         "pool_selection":{"workers":workers,"requested_workers_if_under_threshold":WORKERS_IF_SMALL,
                           "workers_if_estimated_single_worker_rss_exceeds_threshold":WORKERS_IF_LARGE,
                           "estimated_single_worker_peak_rss_bytes":estimated_worker_bytes,
                           "switch_threshold_bytes":RSS_THRESHOLD_BYTES,
                           "memory_limit_bytes":MEMORY_LIMIT_BYTES,
                           "estimate_basis":"6x2^20 Sobol float64 matrix is 48 MiB, plus Gamma/binomial/median vectors and observed 2^17 worker RSS; conservative estimate exceeds 220 MB, so use 8 workers."},
         "memory_measured":{"max_single_worker_peak_rss_bytes":maxrss,"parent_peak_rss_bytes":parentrss,
                            "conservative_workers_times_max_plus_parent_bytes":workers*maxrss+parentrss,
                            "measured_concurrent_process_tree_peak_rss_bytes":MEMORY_OBS["process_tree_peak_rss_bytes"]},
         "baseline_r0_by_n":{str(n):mainmod.compact_float(float(next(x["r0"] for x in json.loads((HERE/"nonlinear_search.json").read_text())["baseline_medians"] if x["n"]==n))) for n in NS},
         "controls_by_n":controls,
         "results":results,"seed_range_aggregates":aggregates,
         "sample_power_level_differences":power_differences,
         "status_counts":status_counts,"wall_seconds":time.time()-start}
    (HERE/"stability.json").write_text(json.dumps(out,separators=(",",":"))+"\n")
    print(json.dumps({"status":out["status"],"results":len(results),"status_counts":status_counts,
                      "workers":workers,"memory":out["memory_measured"],"wall_seconds":out["wall_seconds"]}),
          flush=True)


if __name__=="__main__":
    main()
