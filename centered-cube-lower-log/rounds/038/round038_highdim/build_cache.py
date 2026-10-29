"""One preallocated source-feature cache. Import never allocates the full cache.

build_cache() is intentionally called only by the parent's driver. Running this
file directly writes the small preflight JSON and does NOT build any features.
All numerical error guards retain source_tables' ordinary-FP qualification.
"""
from pathlib import Path
from fractions import Fraction
import hashlib
import json
import math
import resource
import sys
import time
import numpy as np
from source_tables import Source

HERE = Path(__file__).resolve().parent
SOURCE_SHA256 = "85c592ac80d54b376a43038af7395cfc757482b4851b996ffe4d72be2a2fa12d"
FEATURE_NAMES = ["log_total_nominal", "plus_probability_nominal",
                 "bernoulli_variance_nominal", "log1p_relative_class_variance_cap"]


def _source_hash():
    actual = hashlib.sha256((HERE/"source_tables.py").read_bytes()).hexdigest()
    if actual != SOURCE_SHA256:
        raise RuntimeError("Frozen source_tables.py hash changed.")
    return actual


def exact_grid(src):
    natural = {Fraction(10,9)+Fraction(j,src.K*src.v)
               for j in range(src.K*src.v//9+1)}
    old = {Fraction(80+j,72) for j in range(9)}
    grid = sorted(natural | old)
    assert len(grid) == 233
    assert grid[0] == Fraction(10,9) and grid[-1] == Fraction(11,9)
    assert all((length*src.Q).denominator == 1 for length in grid)
    return grid


def _rss():
    # Darwin ru_maxrss is bytes; Linux uses KiB. This is lifetime process peak,
    # not a cache-only incremental allocation or a live resident-set reading.
    value = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return {"raw": value, "unit": "bytes" if sys.platform=="darwin" else "KiB",
            "platform": sys.platform, "process_lifetime_peak": True}


def _metadata(src, grid):
    r, j, q = src.R, len(grid), src.Q
    feature_bytes = 4*r*j*q*np.dtype(np.float64).itemsize
    # norms plus six reported (R,J) numerical error arrays.
    small_array_bytes = 7*r*j*np.dtype(np.float64).itemsize
    return {"source_SHA256": _source_hash(),
            "cache_builder_SHA256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            "n": src.n, "K": src.K, "Q": src.Q, "R": src.R, "J": j,
            "v": src.v, "w": src.w, "delta": src.delta, "p": src.p,
            "pchi": src.pchi, "field": src.field,
            "positive_background_lower_bound": src.background_lower,
            "feature_shape": [4,r*j,q], "feature_dtype": "float64",
            "feature_names": FEATURE_NAMES,
            "flattened_target_order": "row=b*J+j; b outer, exact-grid j inner",
            "feature_array_bytes": feature_bytes,
            "feature_array_GiB": feature_bytes/(1024**3),
            "small_returned_array_bytes": small_array_bytes,
            "estimated_returned_array_bytes": feature_bytes+small_array_bytes,
            "transient_policy": "one evaluate dict and Q-vectors; no list/stack or R/J evaluator cache",
            "grid_first": str(grid[0]), "grid_last": str(grid[-1]),
            "natural_grid_count": src.K*src.v//9+1, "old_grid_count": 9,
            "grid_as_fractions": list(map(str,grid)),
            "normalizer_definition": "mean(exp(features[0,b*J+j,:]))",
            "error_qualification": "analytic midpoint bound plus separately labelled empirical FP guard; not interval certified",
            "python": sys.executable, "python_version": sys.version,
            "numpy": np.__version__}


def preflight():
    """Construct source metadata/grid only: no Source.evaluate or feature allocation."""
    started = time.perf_counter()
    src = Source(R=4,w=4,v=8)
    grid = exact_grid(src)
    meta = _metadata(src,grid)
    return {"scope": "metadata/grid preflight only; full cache NOT built",
            "source_evaluate_calls": 0, "full_feature_allocations": 0,
            **meta, "preflight_seconds": time.perf_counter()-started,
            "ru_maxrss": _rss()}


def build_cache():
    """Return src, exact Fraction grid, features, norms, errors (with errors['meta']).

    Own exactly one (4,R*J,Q) float64 allocation. Nothing is automatically
    written to disk; no samples or RNG calls are made. Caller owns its lifetime.
    """
    started = time.perf_counter()
    src = Source(R=4,w=4,v=8)
    grid = exact_grid(src)
    r,j,q = src.R,len(grid),src.Q
    meta = _metadata(src,grid)
    before = _rss()
    features = np.empty((4,r*j,q),dtype=np.float64)
    norms = np.empty((r,j),dtype=np.float64)
    errors = {name: np.empty((r,j),dtype=np.float64) for name in [
        "mean_error_abs", "analytic_mean_error_abs", "floating_point_guard_abs",
        "min_total_nom", "loginflate", "p_sum_error"]}
    max_normalizer_discrepancy = 0.
    for b in range(r):
        for index,length in enumerate(grid):
            row = b*j+index
            evaluated = src.evaluate(float(length),b)
            total = evaluated["mean_plus"]+evaluated["mean_minus"]
            if not np.all(np.isfinite(total)) or total.min()<=0:
                raise ArithmeticError("Nonpositive/nonfinite nominal class total.")
            np.log(total,out=features[0,row])
            np.divide(evaluated["mean_plus"],total,out=features[1,row])
            np.multiply(features[1,row],1-features[1,row],out=features[2,row])
            np.log1p(evaluated["relative_variance_class_max"],out=features[3,row])
            # Reconstruct exactly the nominal tilt used by the count proposal;
            # mean(total) can differ in its last bit from mean(exp(log(total))).
            norm = float(np.exp(features[0,row]).mean())
            norms[b,index] = norm
            e = float(np.max(evaluated["mean_error_abs"]))
            min_s = float(total.min())
            ratio = 2*e/min_s
            if not (0<=ratio<1 and 0<norm and math.isfinite(norm)):
                raise ArithmeticError("Invalid error fraction or tilt normalizer.")
            errors["mean_error_abs"][b,index] = e
            errors["analytic_mean_error_abs"][b,index] = float(np.max(
                evaluated["analytic_mean_error_abs"]))
            errors["floating_point_guard_abs"][b,index] = float(np.max(
                evaluated["floating_point_guard_abs"]))
            errors["min_total_nom"][b,index] = min_s
            errors["loginflate"][b,index] = src.n*(-math.log1p(-ratio))
            errors["p_sum_error"][b,index] = src.n*e/min_s
            max_normalizer_discrepancy = max(max_normalizer_discrepancy,
                abs(norm-float(evaluated["normalizer"])))
            del evaluated,total
    log_c_z = max(0.,src.n*float(np.log(norms).max()))
    log_c_e = float(errors["loginflate"].max())
    errors["C_Z"] = math.exp(log_c_z)
    errors["C_E"] = math.exp(log_c_e)
    errors["log_C_Z"] = log_c_z
    errors["log_C_E"] = log_c_e
    meta.update({"scope": "actual full feature cache build; no samples",
        "build_runtime_seconds": time.perf_counter()-started,
        "source_evaluate_calls": r*j,
        "actual_normalizer_max_abs_residual": float(np.abs(norms-1).max()),
        "actual_normalizer_min": float(norms.min()),
        "actual_normalizer_max": float(norms.max()),
        "source_vs_feature_normalizer_max_abs_discrepancy": max_normalizer_discrepancy,
        "actual_returned_array_bytes": features.nbytes+norms.nbytes+
             sum(v.nbytes for v in errors.values() if isinstance(v,np.ndarray)),
        "ru_maxrss_before_build": before, "ru_maxrss_after_build": _rss(),
        "cache_saved_to_disk": False, "random_samples": 0})
    errors["meta"] = meta
    return src,grid,features,norms,errors


if __name__=="__main__":
    target = HERE/"CACHE_PREFLIGHT.json"
    if target.exists():
        raise RuntimeError("Preflight already saved; refusing replacement.")
    result = preflight()
    target.write_text(json.dumps(result,separators=(",",":")))
    print(json.dumps({k:result[k] for k in ["scope","J","Q","feature_shape",
        "feature_array_bytes","estimated_returned_array_bytes","preflight_seconds",
        "ru_maxrss","cache_builder_SHA256"]},indent=2))
