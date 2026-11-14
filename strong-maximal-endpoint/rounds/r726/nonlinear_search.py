#!/usr/bin/env python3
"""Finite-angle proxy optimization and randomized-Sobol median witnesses."""
from __future__ import annotations

import json
import math
import os
import resource
import sys
import time
import traceback
import cmath
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

for _name in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS",
              "VECLIB_MAXIMUM_THREADS", "NUMEXPR_MAXIMUM_THREADS"):
    os.environ[_name] = "1"

import numpy as np
from scipy import optimize, special, stats
from scipy.optimize import minimize
from scipy.stats import qmc

HERE = Path(__file__).resolve().parent
NS = (16, 128, 1024, 8192)
EPSILONS = (1.0, 4.0, 16.0)
PATTERNS = ("sparse", "dense_energy")
SAMPLE_POWERS = (15, 17)
WORKERS = 12
BASE_SEED = 726202611
ANGLE_COUNT = 64
TWO_PI = 2 * math.pi
STARTS = ((None, None), (0.0, None), (None, 0.0), (4.0, 4.0))
BOUND_TOL = 1e-7


def fixed_mu():
    return optimize.brentq(lambda x: x - 1.0 - math.exp(-x), 1.0, 2.0,
                           xtol=1e-14, rtol=1e-14)


MU = fixed_mu()
P_MU = 1.0 / (1.0 + math.exp(MU))


def cdf_gamma_binomial(n, s, c, p, weights=None):
    if s <= 0:
        return 0.0
    js = np.arange(n + 1, dtype=float)
    if weights is None:
        weights = stats.binom.pmf(np.arange(n + 1), n, p)
    x = np.maximum(s - c * js, 0.0)
    return float(np.dot(weights, special.gammainc(n, x)))


def baseline_median(n):
    weights = stats.binom.pmf(np.arange(n + 1), n, P_MU)
    weight_resid = float(weights.sum() - 1.0)
    weights /= weights.sum()
    mean = n * (1.0 + MU * P_MU)
    variance = n * (1.0 + MU**2 * P_MU * (1 - P_MU))
    lo = 0.0
    hi = mean + 12.0 * math.sqrt(variance)
    root = optimize.brentq(lambda s: cdf_gamma_binomial(n, s, MU, P_MU, weights) - 0.5,
                           lo, hi, xtol=1e-12, rtol=1e-14)
    residual = cdf_gamma_binomial(n, root, MU, P_MU, weights) - 0.5
    return {"n": n, "mu": MU, "p_mu": P_MU, "r0": root,
            "cdf_residual": residual, "binomial_weight_sum_residual": weight_resid,
            "bracket": [lo, hi]}


def one_site_mean(c):
    return 1.0 + c / (1.0 + math.exp(c))


def phi(c, frequency):
    p = 1.0 / (1.0 + math.exp(c))
    return (1.0 + cmath.exp(complex(-c, frequency * c))) / (
        (1.0 + math.exp(-c)) * (1.0 - 1j * frequency))


def proxy_value(cp, cm, theta, eps, a, b, z, frequency):
    value = (a * one_site_mean(cp) + b * one_site_mean(cm) + z * MU -
             eps * np.real(np.exp(1j * theta) *
                           phi(cp, frequency)**a * phi(cm, -frequency)**b))
    return float(value)


def optimize_angles(n, pattern, eps, a, b, z, frequency):
    angles = np.arange(ANGLE_COUNT, dtype=float) * TWO_PI / ANGLE_COUNT
    plus, minus, values, winner, all_starts, baseline_margins = [], [], [], [], [], []
    for theta in angles:
        start_points = ((MU, MU), (0.0, MU), (MU, 0.0), (4.0, 4.0))
        results = []
        for start in start_points:
            # The target is the largest geometric proxy. SciPy minimizes, so
            # optimize its negative and store the corresponding proxy value.
            result = minimize(lambda x: -proxy_value(x[0], x[1], theta, eps,
                                                     a, b, z, frequency),
                              x0=np.array(start), method="L-BFGS-B",
                              bounds=((0.0, 8.0), (0.0, 8.0)),
                              options={"ftol": 1e-13, "gtol": 1e-9, "maxiter": 300,
                                       "maxfun": 1500})
            x0, x1 = map(float, result.x)
            bitmask = ((1 if x0 <= BOUND_TOL else 0) |
                       (2 if x0 >= 8.0 - BOUND_TOL else 0) |
                       (4 if x1 <= BOUND_TOL else 0) |
                       (8 if x1 >= 8.0 - BOUND_TOL else 0))
            results.append((-float(result.fun), x0, x1, bool(result.success),
                            int(result.status), int(result.nfev), bitmask))
        finite = [i for i, res in enumerate(results) if math.isfinite(res[0])]
        if not finite:
            raise ArithmeticError(f"all proxy optimizations nonfinite at theta={theta}")
        win = max(finite, key=lambda i: results[i][0])
        best = results[win]
        baseline_proxy = proxy_value(MU, MU, theta, eps, a, b, z, frequency)
        if best[0] < baseline_proxy - 1e-8:
            raise AssertionError(
                f"optimized proxy below feasible baseline at n={n}, eps={eps}, theta={theta}: "
                f"{best[0]} < {baseline_proxy}")
        baseline_margins.append(best[0] - baseline_proxy)
        plus.append(best[1]); minus.append(best[2]); values.append(best[0]); winner.append(win)
        # [success, status, nfev, boundary bitmask] for each of four starts.
        all_starts.append([[int(r[3]), r[4], r[5], r[6]] for r in results])
    return angles, plus, minus, values, winner, all_starts, baseline_margins


def gamma_inverse(shape, uniforms):
    if shape <= 0:
        return np.zeros_like(uniforms)
    clipped = np.clip(uniforms, np.finfo(float).tiny, 1.0 - np.finfo(float).eps)
    return special.gammaincinv(shape, clipped)


def binomial_quantile_from_uniforms(n, p, uniforms):
    pmf = stats.binom.pmf(np.arange(n + 1), n, p)
    total = float(pmf.sum())
    if not math.isfinite(total) or total <= 0 or not np.isfinite(pmf).all():
        raise FloatingPointError(f"invalid binomial pmf n={n},p={p},sum={total}")
    cdf = np.cumsum(pmf / total)
    cdf[-1] = 1.0
    values = np.searchsorted(cdf, uniforms, side="left")
    return values.astype(np.int32), abs(total - 1.0)


def sobol_components(n, pattern, power, seed):
    size = 1 << power
    sampler = qmc.Sobol(d=6, scramble=True, seed=seed)
    u = sampler.random_base2(power)
    if pattern == "sparse":
        a = b = 1
        z = n - 2
        frequency = 1.0
    else:
        a = b = n // 2
        z = 0
        frequency = 1.0 / math.sqrt(n)
    ga = gamma_inverse(a, u[:, 0])
    ja_u = u[:, 1]
    gb = gamma_inverse(b, u[:, 2])
    jb_u = u[:, 3]
    gz = gamma_inverse(z, u[:, 4]) if z else np.zeros(size)
    if not (np.isfinite(ga).all() and np.isfinite(gb).all() and np.isfinite(gz).all()):
        raise FloatingPointError(f"nonfinite Gamma inverse n={n}, pattern={pattern}, seed={seed}")
    jz, residual_z = binomial_quantile_from_uniforms(z, P_MU, u[:, 5]) if z else (np.zeros(size, dtype=np.int32), 0.0)
    return {"size": size, "ga": ga, "ja_u": ja_u, "gb": gb, "jb_u": jb_u,
            "gz": gz, "jz": jz, "a": a, "b": b, "z": z,
            "frequency": frequency, "pmf_norm_resid_z": residual_z}


def sample_medians(n, pattern, eps, angles, cps, cms, power, seed):
    comp = sobol_components(n, pattern, power, seed)
    medians = []
    pmf_resid_max = comp["pmf_norm_resid_z"]
    all_finite = True
    for theta, cp, cm in zip(angles, cps, cms):
        pa = 1.0 / (1.0 + math.exp(cp))
        pb = 1.0 / (1.0 + math.exp(cm))
        ja, ra = binomial_quantile_from_uniforms(comp["a"], pa, comp["ja_u"])
        jb, rb = binomial_quantile_from_uniforms(comp["b"], pb, comp["jb_u"])
        pmf_resid_max = max(pmf_resid_max, ra, rb)
        s = (comp["ga"] + cp * ja + comp["gb"] + cm * jb +
             comp["gz"] + MU * comp["jz"])
        w = comp["frequency"] * (comp["ga"] + cp * ja - comp["gb"] - cm * jb)
        sample = s - eps * np.cos(theta + w)
        finite = np.isfinite(sample).all()
        all_finite = all_finite and bool(finite)
        if not finite:
            raise FloatingPointError(f"nonfinite witness sample n={n},pattern={pattern},theta={theta}")
        # np.median uses the average of the two central order statistics for
        # these even sample sizes.
        medians.append(float(np.median(sample)))
    del comp
    return medians, pmf_resid_max, all_finite


def logmeanexp(values):
    values = np.asarray(values, dtype=float)
    shift = float(np.max(values))
    return shift + math.log(float(np.exp(values - shift).mean()))


def compact_float(x):
    # Sobol medians are only numerical diagnostics; 10 significant digits
    # retain ample precision while keeping the complete 24x64 result compact.
    return format(float(x), ".10g")


def peak_rss_bytes():
    value = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return int(value if sys.platform == "darwin" else value * 1024)


def run_parameter(task):
    n, pattern, eps, r0, param_index = task
    start = time.monotonic()
    try:
        if pattern == "sparse":
            a, b, z, frequency = 1, 1, n - 2, 1.0
        else:
            a, b, z, frequency = n // 2, n // 2, 0, 1.0 / math.sqrt(n)
        angles, cps, cms, proxy, winners, starts, baseline_margins = optimize_angles(
            n, pattern, eps, a, b, z, frequency)
        replicate_data = []
        for repidx, power in enumerate(SAMPLE_POWERS):
            seed = BASE_SEED + n * 10000 + (0 if pattern == "sparse" else 4000) + int(eps) * 100 + repidx
            medians, pmf_resid, finite = sample_medians(
                n, pattern, eps, angles, cps, cms, power, seed)
            logi0 = math.log(float(special.i0e(eps))) + abs(eps)
            score = logmeanexp(np.asarray(medians) - MU * n) - logi0
            even_score = logmeanexp(np.asarray(medians)[::2] - MU * n) - logi0
            replicate_data.append({
                "sample_power": power, "sample_size": 1 << power, "seed": seed,
                "actual_medians": [compact_float(v) for v in medians],
                "logscore": compact_float(score), "logscore_even32": compact_float(even_score),
                "even32_minus_full64": compact_float(even_score - score),
                "logscore_minus_mu_n": compact_float(score),
                "normalized_logscore_minus_mu_n_over_logn": compact_float(score / math.log(n)),
                "relative_logscore_minus_baseline_r0": compact_float(score + MU * n - r0),
                "relative_to_r0_positive": bool(score + MU * n > r0),
                "max_binomial_pmf_normalization_residual": compact_float(pmf_resid),
                "all_samples_finite": finite,
            })
        winner_start_status_counts = {
            "all_4_starts_success_at_angle_count": sum(1 for row in starts if all(x[0] for x in row)),
            "any_start_non_success_at_angle_count": sum(1 for row in starts if any(not x[0] for x in row)),
            "any_start_boundary_hit_at_angle_count": sum(1 for row in starts if any(x[3] for x in row)),
            "best_winner_boundary_hit_angle_count": sum(1 for i, row in enumerate(starts)
                                                       if row[winners[i]][3]),
            "best_winner_non_success_angle_count": sum(1 for i, row in enumerate(starts)
                                                        if not row[winners[i]][0]),
        }
        lower_cp = [i for i, v in enumerate(cps) if v <= BOUND_TOL]
        lower_cm = [i for i, v in enumerate(cms) if v <= BOUND_TOL]
        positive = [r for r in replicate_data if r["relative_to_r0_positive"]]
        return {"n": n, "pattern": pattern, "epsilon": eps,
                "counts_plus_minus_zero": [a, b, z], "frequency_t": compact_float(frequency),
                "status": "ok", "elapsed_seconds": time.monotonic() - start,
                "worker_peak_rss_bytes": peak_rss_bytes(),
                "optimizer": {"method": "SciPy L-BFGS-B", "bounds": [[0, 8], [0, 8]],
                              "starts": [[compact_float(MU), compact_float(MU)], ["0", compact_float(MU)],
                                         [compact_float(MU), "0"], ["4", "4"]],
                              "success_code": 0, "boundary_tolerance": BOUND_TOL,
                              "objective_direction": "maximize supplied proxy via minimize(-proxy)",
                              "baseline_feasible_check": {"assertion": "best proxy >= value at (mu,mu)-1e-8 for every angle",
                                                           "passed": True,
                                                           "minimum_best_minus_baseline_proxy": compact_float(min(baseline_margins))},
                              "winner_start_indices_zero_based": winners,
                              "per_start_status_success_boundarymask": [[[s[0], s[1], s[3]] for s in row]
                                                                          for row in starts],
                              "summary": winner_start_status_counts},
                "angle_nodes": [compact_float(v) for v in angles],
                "best_proxy_cplus": [compact_float(v) for v in cps],
                "best_proxy_cminus": [compact_float(v) for v in cms],
                "c_zero_boundary_diagnostics": {
                    "cplus_zero_hit_angle_indices": lower_cp,
                    "cminus_zero_hit_angle_indices": lower_cm,
                    "negative_side_unsearched": bool(lower_cp or lower_cm),
                },
                "best_proxy_values": [compact_float(v) for v in proxy],
                "replicates": replicate_data,
                "positive_candidate_replicates": [x["sample_power"] for x in positive]}
    except Exception as exc:
        return {"n": n, "pattern": pattern, "epsilon": eps,
                "status": "error", "elapsed_seconds": time.monotonic() - start,
                "worker_peak_rss_bytes": peak_rss_bytes(),
                "error": f"{type(exc).__name__}: {exc}",
                "traceback": traceback.format_exc(limit=4)}


def main():
    start = time.time()
    baselines = {n: baseline_median(n) for n in NS}
    params = [(n, pattern, eps, baselines[n]["r0"], i)
              for i, (n, pattern, eps) in enumerate(
                  ( (n, pattern, eps) for n in NS for pattern in PATTERNS for eps in EPSILONS))]
    results = []
    with ProcessPoolExecutor(max_workers=WORKERS) as pool:
        futures = {pool.submit(run_parameter, task): task[:3] for task in params}
        for fut in as_completed(futures):
            try:
                result = fut.result()
            except Exception as exc:
                n, pattern, eps = futures[fut]
                result = {"n": n, "pattern": pattern, "epsilon": eps,
                          "status": "worker_error", "error": f"{type(exc).__name__}: {exc}"}
            results.append(result)
            print(json.dumps({"completed": len(results), "total": len(params),
                              "n": result.get("n"), "pattern": result.get("pattern"),
                              "epsilon": result.get("epsilon"),
                              "status": result.get("status"),
                              "elapsed_seconds": result.get("elapsed_seconds")}), flush=True)
            write_output(baselines, results, start, partial=True)
    results.sort(key=lambda x: (x.get("n", 0), PATTERNS.index(x.get("pattern", PATTERNS[0])),
                                EPSILONS.index(x.get("epsilon", EPSILONS[0]))))
    positives = [{"n": x["n"], "pattern": x["pattern"], "epsilon": x["epsilon"],
                  "positive_sample_powers": x.get("positive_candidate_replicates", [])}
                 for x in results if x.get("positive_candidate_replicates")]
    write_output(baselines, results, start, partial=False, positives=positives)


def write_output(baselines, results, start, partial, positives=None):
    output = {
        "status": "partial" if partial else
                  ("complete" if len(results) == len(NS) * len(PATTERNS) * len(EPSILONS) else "incomplete"),
        "purpose": "Finite nonlinear proxy optimization plus sampled actual-median witnesses; numerical exploration only, not a proof.",
        "dimensions": list(NS), "patterns": list(PATTERNS), "epsilons": list(EPSILONS),
        "expected_parameter_count": len(NS) * len(PATTERNS) * len(EPSILONS),
        "theta_nodes": ANGLE_COUNT, "sample_powers": list(SAMPLE_POWERS),
        "sobol": {"dimension": 6, "scrambled": True, "base_seed": BASE_SEED,
                  "seed_formula": "base+n*10000+(0 sparse,4000 dense_energy)+int(epsilon)*100+replicate_index"},
        "frequency_assumption": "t=1 for sparse +/-1; t=1/sqrt(n) for dense_energy +/-1/sqrt(n), per parent clarification.",
        "mu": compact_float(MU), "p_mu": compact_float(P_MU),
        "baseline_medians": [{k: compact_float(v) if isinstance(v, float) else v
                               for k, v in rec.items()} for rec in baselines.values()],
        "baseline_model": "Finite Gamma(n,1)+mu*Binomial(n,p_mu) CDF solved independently by SciPy brentq; finite binomial-gamma sum.",
        "proxy_model": "max over cplus,cminus in [0,8] of the supplied geometric proxy, using SciPy minimize(-proxy) with the -epsilon Re phase term.",
        "proxy_search_domain": [0, 8],
        "negative_side_scope": "Negative c is legal in the broader model but not searched in this authorized [0,8] main batch. Every zero-boundary hit is recorded per parameter; absence of positive relative-to-baseline scores here is not an obstruction for negative controls.",
        "memory_diagnostics": {
            "max_measured_single_worker_peak_rss_bytes": max((x.get("worker_peak_rss_bytes", 0) for x in results), default=0),
            "parent_peak_rss_bytes": peak_rss_bytes(),
            "conservative_12_worker_plus_parent_bytes": (12 * max((x.get("worker_peak_rss_bytes", 0) for x in results), default=0) + peak_rss_bytes()),
            "limit_bytes": 3000000000,
        },
        "proxy_objective": "a M(cplus)+b M(cminus)+z mu-eps Re(exp(i theta) phi_cplus(t)^a phi_cminus(-t)^b)",
        "actual_model": "For each proxy-optimal node, sample S and W from grouped Gamma/Binomial variables with common scrambled-Sobol uniforms; compute sample median of S-epsilon*cos(theta+W).",
        "interpolation_scope": "64 proxy-optimal values are interpreted as periodic piecewise-linear cplus(theta),cminus(theta); the resulting witness family is Borel. Sample medians and angular quadrature are not rigorous bounds.",
        "actual_witness_warning": "Actual r_theta is the sample median of S-epsilon*cos(theta+W); it is not the proxy minimum value.",
        "sample_policy": "Gamma sums and Sobol base uniforms are generated once per parameter/replicate and reused for all 64 angles; binomial inverse transforms are recomputed for each node's optimized controls. No sample arrays are saved.",
        "processes": WORKERS, "blas_threads_per_worker": 1,
        "results": results, "positive_candidate_summaries": positives or [],
        "status_counts": {s: sum(1 for x in results if x.get("status") == s)
                          for s in ("ok", "error", "worker_error")},
        "wall_seconds": time.time() - start,
    }
    # Compact serialization keeps the complete 64-node optimizer and median
    # diagnostics below the requested storage budget.
    (HERE / "nonlinear_search.json").write_text(json.dumps(output, separators=(",", ":")) + "\n")
    print(json.dumps({"status": output["status"], "completed": len(results),
                      "status_counts": output["status_counts"],
                      "positive_candidate_count": len(output["positive_candidate_summaries"]),
                      "wall_seconds": output["wall_seconds"]}), flush=True)


if __name__ == "__main__":
    main()
