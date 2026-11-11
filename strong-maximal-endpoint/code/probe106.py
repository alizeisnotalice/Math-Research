#!/usr/bin/env python3
"""Finite deterministic CDF diagnostics on the expanded r106 left-tail window.

The Gamma/geometric/binomial positive-convolution kernels are imported from
the previously checked r094 implementation. Results are numerical probes,
not a proof of a uniform constant or a general maximal-operator bound.
"""
from __future__ import annotations

import concurrent.futures as cf
import json
import math
import multiprocessing as mp
import os
import sys
import time
from pathlib import Path

for _name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
              "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[_name] = "1"


import tilted_probe94 as kernel
from scipy.optimize import brentq


NS = (64, 128, 256, 512, 1024, 2048, 4096, 8192, 16384, 32768, 65536)
Z_SCALES = (0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0)
FOREGROUND_LABELS = ("1", "2", "4", "8", "n^0.25")
TAIL_TARGET = 1e-20
WORKERS = 12
RESOLUTION_RELATIVE_TAIL_MAX = 1e-6


def positive_root(fun, hi: float) -> float:
    while fun(hi) > 0:
        hi *= 2.0
        if hi > 1e4:
            raise ArithmeticError("positive saddlepoint root failed to bracket")
    return float(brentq(fun, 0.0, hi, xtol=5e-15,
                        rtol=4 * kernel.np.finfo(float).eps))


def one_case(spec: tuple[int, float, str]) -> dict:
    n, z_scale, foreground_label = spec
    started = time.monotonic()

    L = math.log(n)
    z = z_scale * math.sqrt(L)
    q = kernel.MU*n - kernel.SIGMA*math.sqrt(n)*z
    if q <= 0:
        return {"status": "outside_positive_support", "N": n,
                "z_scale": z_scale, "z": z, "q": q}

    def ref_mean_residual(s: float) -> float:
        mean, _, _ = kernel.mean_tilt(s)
        return n*mean-q

    sref = positive_root(ref_mean_residual, max(0.5, math.sqrt(L/n)))
    ref_saddle_residual = ref_mean_residual(sref)
    c0, c0_root_residual = kernel.action_c(sref)
    Fref, ref_tail = kernel.reference_cdf(n, q, c0, TAIL_TARGET)

    foreground = n**0.25 if foreground_label == "n^0.25" else float(foreground_label)
    bulk = (n-foreground)/(n-1.0)
    if not (0 < foreground < n and bulk > 0):
        raise ArithmeticError(f"invalid weights: N={n}, front={foreground}, bulk={bulk}")

    def weighted_mean_residual(s: float) -> float:
        front_mean, _, _ = kernel.mean_tilt(s*foreground)
        bulk_mean, _, _ = kernel.mean_tilt(s*bulk)
        return foreground*front_mean + (n-1.0)*bulk*bulk_mean-q

    sweighted = positive_root(weighted_mean_residual, max(0.5, sref))
    weighted_saddle_residual = weighted_mean_residual(sweighted)
    c1, c1_root_residual = kernel.action_c(sweighted*foreground)
    c2, c2_root_residual = kernel.action_c(sweighted*bulk)

    if foreground == 1.0:
        # This is exactly the iid vector, so the distributional identity gives
        # the reference CDF without adding another floating-point convolution.
        Fweighted = Fref
        weighted_tail = dict(ref_tail)
        weighted_tail["cdf_evaluation"] = "exact iid identity"
    else:
        Fweighted, weighted_tail = kernel.weighted_cdf(
            n, q, foreground, bulk, c1, c2, TAIL_TARGET)
        weighted_tail["cdf_evaluation"] = "r094 positive Gamma/Geometric/Binomial mixture"

    if not (0.0 < Fweighted < 1.0 and 0.0 < Fref < 1.0):
        raise ArithmeticError(f"CDF outside (0,1): N={n}, z={z_scale}, a={foreground}, "
                              f"Fw={Fweighted}, F*={Fref}")
    epsw = float(weighted_tail["combined_tail_bound"])
    epsr = float(ref_tail["combined_tail_bound"])
    tail_over_fw = epsw/Fweighted
    tail_over_fr = epsr/Fref

    ratio = Fweighted/Fref
    ratio_lo = Fweighted/(Fref+epsr)
    ratio_hi = (Fweighted+epsw)/Fref
    ratio_width_rel = (ratio_hi-ratio_lo)/ratio if ratio else math.inf
    resolved = max(tail_over_fw, tail_over_fr, ratio_width_rel) <= RESOLUTION_RELATIVE_TAIL_MAX
    log_ratio = math.log(Fweighted)-math.log(Fref)
    scaled_log_ratio = math.sqrt(n)*log_ratio
    c_required_diagnostic = -scaled_log_ratio
    thresholds = {str(C): math.exp(-C/math.sqrt(n)) for C in (1.0, 4.0, 16.0)}

    return {
        "status": "pass",
        "N": n,
        "logN": L,
        "z_scale": z_scale,
        "z": z,
        "q": q,
        "lambda_front": foreground,
        "lambda_front_label": foreground_label,
        "lambda_bulk": bulk,
        "reference_saddlepoint": sref,
        "reference_saddlepoint_residual": ref_saddle_residual,
        "reference_action_c": c0,
        "reference_action_root_residual": c0_root_residual,
        "weighted_saddlepoint": sweighted,
        "weighted_saddlepoint_residual": weighted_saddle_residual,
        "weighted_front_action_c": c1,
        "weighted_front_action_root_residual": c1_root_residual,
        "weighted_bulk_action_c": c2,
        "weighted_bulk_action_root_residual": c2_root_residual,
        "F_weighted_numeric": Fweighted,
        "F_reference_numeric": Fref,
        "F_weighted_over_F_reference_numeric": ratio,
        "ratio_envelope_from_numerical_truncation_budget": [ratio_lo, ratio_hi],
        "log_ratio_numeric": log_ratio,
        "sqrtN_log_ratio_numeric": scaled_log_ratio,
        "C_required_diagnostic_only": c_required_diagnostic,
        "exp_minus_C_over_sqrtN_diagnostics": {
            C: {"threshold": val, "ratio_at_least_threshold": ratio_lo >= val}
            for C, val in thresholds.items()
        },
        "weighted_absolute_tail_budget_numeric": epsw,
        "reference_absolute_tail_budget_numeric": epsr,
        "tail_budget_over_F_weighted": tail_over_fw,
        "tail_budget_over_F_reference": tail_over_fr,
        "ratio_interval_relative_width": ratio_width_rel,
        "resolution_status": "resolved_for_truncation_diagnostic" if resolved else "finite_resolution",
        "weighted_tail_components": weighted_tail,
        "reference_tail_components": ref_tail,
        "arithmetic_scope": "SciPy double-precision finite positive sums; tail budgets are numerical diagnostics, not directed intervals.",
        "seconds": time.monotonic()-started,
        "pid": os.getpid(),
    }


def main() -> None:
    out_dir = Path(sys.argv[1]); out_dir.mkdir(parents=True, exist_ok=True)
    specs = [(n, z, label) for n in NS for z in Z_SCALES for label in FOREGROUND_LABELS]
    started = time.monotonic()
    context = mp.get_context("spawn")
    with cf.ProcessPoolExecutor(max_workers=WORKERS, mp_context=context) as pool:
        rows = list(pool.map(one_case, specs, chunksize=1))
    pids = sorted({row.get("pid") for row in rows if row.get("pid") is not None})
    resolved_count = sum(row.get("resolution_status") == "resolved_for_truncation_diagnostic" for row in rows)
    result = {
        "schema": "r106-expanded-left-tail-positive-convolution-v1",
        "pass_all_computations": len(rows) == len(specs) and all(row.get("status") == "pass" for row in rows),
        "case_count": len(rows),
        "expected_case_count": len(specs),
        "resolved_for_truncation_diagnostic_count": resolved_count,
        "finite_resolution_count": len(rows)-resolved_count,
        "workers_configured": WORKERS,
        "worker_pids": pids,
        "actual_workers": len(pids),
        "n_values": list(NS),
        "z_scale_values": list(Z_SCALES),
        "foreground_values": list(FOREGROUND_LABELS),
        "duplicate_foreground_aliases": "At N=256, n^0.25=4; at N=4096, n^0.25=8. Both requested labels are retained, so these repeated parameter values are deliberate.",
        "tail_target_each_component": TAIL_TARGET,
        "resolution_relative_tail_threshold": RESOLUTION_RELATIVE_TAIL_MAX,
        "formula_q": "mu*N - sigma*sqrt(N)*z, z=z_scale*sqrt(log(N))",
        "formula_weights": "lambda_front=a; lambda_bulk=(N-a)/(N-1)",
        "action_choice": "positive saddlepoint separately matching q for weighted and iid reference sums",
        "comparison": "F_weighted/F_reference and diagnostics against exp(-C/sqrt(N)) for C=1,4,16; C is not fitted as a proof",
        "scope": "Finite deterministic numerical probe only; no uniform C proof or general maximal-operator bound.",
        "seconds": time.monotonic()-started,
        "rows": rows,
    }
    path = out_dir/"probe106-results.json"
    path.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({k: v for k, v in result.items() if k != "rows"}, indent=2))


if __name__ == "__main__":
    main()
