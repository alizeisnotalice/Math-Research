"""Finite-dimensional, exactly normalized shared-phase mixture probe.

Run with system Python + NumPy. No Gaussian or independent-scale surrogate.
All outputs stay alongside this file.
"""
from pathlib import Path
import json
import math
import sys
import time
import numpy as np

OUT = Path(__file__).resolve().parent
SEED = 20261029031
K = 3
P = 0.5
DELTA = 1 / 192
ATOMS = np.array([1 / 12, 1 / 2, 11 / 12])
PHASES = np.arange(4) / 4
ALPHA = 729 / 29
B1 = 4 * ALPHA * (
    (math.sin(math.pi * 10 / 27) / math.pi) ** 4
    + 3 * (math.sin(math.pi / 9) / math.pi) ** 4
)
C1 = math.sqrt(2 * B1 / 3)
COS_AMP = P / K * math.pi * C1
LS = np.linspace(10 / 9, 11 / 9, 49)
ANCHOR_INDICES = [0, 24, 48]
ANCHORS = LS[ANCHOR_INDICES]
MODELS = ["product", "shared", "independent"]


def lse(z, axis=-1):
    zmax = np.max(z, axis=axis, keepdims=True)
    return np.squeeze(zmax + np.log(np.sum(np.exp(z - zmax), axis=axis,
                                         keepdims=True)), axis=axis)


def mu_primitive(t):
    """Primitive of periodic unit-mass equal-box measure; boxes inside [0,1]."""
    periods = np.floor(t)
    u = t - periods
    result = np.zeros_like(u)
    for a in ATOMS:
        result += np.clip((u - a + DELTA / 2) / DELTA, 0, 1)
    return periods + result / K


def average_from_shifted(x, length, cosine=None):
    """Exact centered interval average of base rho at centers x (any shape)."""
    mass = mu_primitive(x + length / 2) - mu_primitive(x - length / 2)
    if cosine is None:
        cosine = np.cos(2 * np.pi * x)
    return ((1 - P) + P * mass / length
            + (P / K) * C1 * np.sin(np.pi * length) / length * cosine)


def density(x):
    x = x % 1
    result = (1 - P) + COS_AMP * np.cos(2 * np.pi * x)
    for a in ATOMS:
        result += (P / (K * DELTA)) * (np.abs(x - a) <= DELTA / 2)
    return result


def sample_rho(rng, shape):
    count = math.prod(shape)
    atomic = rng.random(count) < P
    values = np.empty(count)
    na = int(atomic.sum())
    values[atomic] = (ATOMS[rng.integers(K, size=na)]
                      + DELTA * (rng.random(na) - 0.5))
    remaining = np.flatnonzero(~atomic)
    ratio = COS_AMP / (1 - P)
    assert 0 <= ratio < 1
    while remaining.size:
        candidate = rng.random(remaining.size)
        accept = rng.random(remaining.size) < (
            1 + ratio * np.cos(2 * np.pi * candidate)) / (1 + ratio)
        values[remaining[accept]] = candidate[accept]
        remaining = remaining[~accept]
    return values.reshape(shape)


def sample_proposal(rng, samples, dimension, proposal):
    x = rng.random((samples, dimension))
    if proposal == "uniform":
        return x
    if proposal != "half_uniform_half_product_tilts":
        raise ValueError(proposal)
    tilted = rng.random(samples) >= 0.5
    rows = np.flatnonzero(tilted)
    source_points = sample_rho(rng, (rows.size, dimension))
    length = ANCHORS[rng.integers(len(ANCHORS), size=rows.size)]
    phase = PHASES[rng.integers(len(PHASES), size=rows.size)]
    x[rows] = (source_points
               - (rng.random(source_points.shape) - 0.5) * length[:, None]
               + phase[:, None]) % 1
    return x


def evaluate_grid(x, lengths, proposal):
    shifted = (x[:, :, None] - PHASES[None, None, :]) % 1
    cosine = np.cos(2 * np.pi * shifted)
    maxima = np.full((x.shape[0], 3), -np.inf)
    proposal_components = []
    for length in lengths:
        y = average_from_shifted(shifted, length, cosine)
        assert np.min(y) > 0
        phase_log_products = np.log(y).sum(axis=1)
        values = np.column_stack([
            phase_log_products[:, 0],
            lse(phase_log_products) - math.log(len(PHASES)),
            np.log(np.mean(y, axis=2)).sum(axis=1),
        ])
        maxima = np.maximum(maxima, values)
        if np.any(np.isclose(length, ANCHORS, atol=1e-14, rtol=0)):
            proposal_components.append(phase_log_products)
    if proposal == "uniform":
        weights = np.ones(x.shape[0])
    else:
        assert len(proposal_components) == len(ANCHORS)
        comp = np.concatenate(proposal_components, axis=1)
        log_tilt_part = lse(comp) - math.log(comp.shape[1])
        log_q = np.logaddexp(math.log(0.5), math.log(0.5) + log_tilt_part)
        weights = np.exp(-log_q)
        assert np.all(weights <= 2 + 1e-12)
    return maxima, weights


def weak_scan(log_values, weights, min_ess):
    order = np.argsort(log_values)[::-1]
    log_values = log_values[order]
    weights = weights[order]
    sums = np.cumsum(weights)
    squares = np.cumsum(weights * weights)
    ess = sums * sums / squares
    values = np.exp(log_values) * sums / len(weights)
    eligible = ess >= min_ess
    if not np.any(eligible):
        return None
    index = int(np.argmax(np.where(eligible, values, -np.inf)))
    prob = sums[index] / len(weights)
    second = squares[index] / len(weights)
    se = math.sqrt(max(0, second - prob * prob) / len(weights))
    return {
        "coefficient": float(values[index]),
        "log_threshold": float(log_values[index]),
        "threshold": float(np.exp(log_values[index])),
        "tail_probability": float(prob),
        "probability_se_descriptive": float(se),
        "tail_ess": float(ess[index]),
        "raw_tail_count": index + 1,
        "min_tail_ess_filter": min_ess,
    }


def tail_at(log_values, weights, threshold):
    values = weights * (log_values > math.log(threshold))
    total = values.sum()
    sq = np.sum(values * values)
    return {"threshold": float(threshold), "coefficient": float(threshold * values.mean()),
            "probability": float(values.mean()),
            "probability_se_descriptive": float(values.std(ddof=1)
                                                / np.sqrt(len(values))),
            "tail_ess": float(total * total / sq) if sq else 0.0,
            "raw_tail_count": int(np.count_nonzero(values))}


def fixed_validation(logs, weights):
    """Fixed thresholds, and discovery/holdout split without more sampling."""
    fixed = [0.5, 0.75, 1.0, 1.5, 2.0, 4.0, 8.0, 16.0]
    cut = len(weights) // 2
    result = {"fixed_thresholds": {},
              "discovery_selected_threshold_then_holdout": {},
              "paired_shared_minus_independent": []}
    # All main-case fixed thresholds: 3 dimensions x 3 models x 8 thresholds.
    result["simultaneous_95pct_Hoeffding_probability_radius_fixed_table"] = (
        2 * math.sqrt(math.log(2 * 72 / 0.05) / (2 * len(weights))))
    for j, model in enumerate(MODELS):
        result["fixed_thresholds"][model] = [
            tail_at(logs[:, j], weights, t) for t in fixed]
        discovery = weak_scan(logs[:cut, j], weights[:cut], 64)
        threshold = discovery["threshold"]
        holdout = tail_at(logs[cut:, j], weights[cut:], threshold)
        # Conditional on discovery, the threshold is fixed for holdout.
        holdout["pointwise_95pct_Hoeffding_probability_radius"] = (
            2 * math.sqrt(math.log(40) / (2 * (len(weights) - cut))))
        result["discovery_selected_threshold_then_holdout"][model] = {
            "discovery": discovery, "holdout": holdout}
    for t in fixed:
        difference = weights * ((logs[:, 1] > math.log(t)).astype(float)
                                - (logs[:, 2] > math.log(t)).astype(float))
        result["paired_shared_minus_independent"].append({
            "threshold": t,
            "coefficient_difference": float(t * difference.mean()),
            "coefficient_difference_se_descriptive":
                float(t * difference.std(ddof=1) / math.sqrt(len(weights)))})
    return result


def summarize(logs, weights, dimension, proposal):
    data = {"n": dimension, "samples": len(weights), "proposal": proposal,
            "estimated_uniform_mass": float(weights.mean()),
            "mass_se": float(weights.std(ddof=1) / np.sqrt(len(weights))),
            "weight_min": float(weights.min()), "weight_max": float(weights.max()),
            "models": {}}
    for j, model in enumerate(MODELS):
        full = weak_scan(logs[:, j], weights, 1)
        supported = weak_scan(logs[:, j], weights, 64)
        block_stats = []
        for block in np.array_split(np.arange(len(weights)), 16):
            block_stats.append(weak_scan(logs[block, j], weights[block], 8)
                               ["coefficient"])
        data["models"][model] = {
            "empirical_unrestricted": full, "tail_ess_ge_64": supported,
            "block_coefficient_quartiles_descriptive":
                np.quantile(block_stats, [0.25, 0.5, 0.75]).tolist(),
            "observed_max_log_average": float(logs[:, j].max()),
            "mean_grid_max_average": float(np.mean(weights * np.exp(logs[:, j]))),
        }
    return data


def run_case(rng, dimension, samples, lengths, proposal, name):
    logs_parts, weights_parts = [], []
    start = time.time()
    for left in range(0, samples, 256):
        count = min(256, samples - left)
        x = sample_proposal(rng, count, dimension, proposal)
        logs, weights = evaluate_grid(x, lengths, proposal)
        logs_parts.append(logs)
        weights_parts.append(weights)
    logs = np.concatenate(logs_parts)
    weights = np.concatenate(weights_parts)
    np.savez_compressed(OUT / f"{name}_samples.npz",
                        log_max=logs, weights=weights, lengths=lengths)
    summary = summarize(logs, weights, dimension, proposal)
    summary["runtime_seconds"] = time.time() - start
    print(json.dumps({"case": name, "runtime": summary["runtime_seconds"],
                      "mass": summary["estimated_uniform_mass"],
                      "coefficients": {m: summary["models"][m]["tail_ess_ge_64"]
                                       ["coefficient"] for m in MODELS}},
                     ensure_ascii=False), flush=True)
    return summary


def direct_cube_check(rng):
    """Direct integration in a real 2D centered cube, against exact formula."""
    center = np.array([[0.18, 0.71]])
    length = 7 / 6
    exact_logs, _ = evaluate_grid(center, np.array([length]), "uniform")
    exact = np.exp(exact_logs[0])
    count = 262144
    u = center + length * (rng.random((count, 2)) - 0.5)
    r = density(u[:, :, None] - PHASES[None, None, :])
    values = np.column_stack([np.prod(r[:, :, 0], axis=1),
                              np.mean(np.prod(r, axis=1), axis=1),
                              np.prod(np.mean(r, axis=2), axis=1)])
    result = {}
    for j, model in enumerate(MODELS):
        mc = float(values[:, j].mean())
        se = float(values[:, j].std(ddof=1) / np.sqrt(count))
        result[model] = {"analytic_average": float(exact[j]),
                         "direct_cube_MC": mc, "MC_se": se,
                         "difference_in_se": (mc - exact[j]) / se}
    return {"center": center[0].tolist(), "L": length, "samples": count,
            "models": result}


def resolution_check(rng):
    x = rng.random((512, 64))
    coarse, _ = evaluate_grid(x, LS, "uniform")
    fine, _ = evaluate_grid(x, np.linspace(10 / 9, 11 / 9, 193), "uniform")
    gap = fine - coarse
    assert gap.min() >= -1e-10
    return {model: {"mean_log_gain_193_vs_49": float(gap[:, j].mean()),
                    "max_log_gain": float(gap[:, j].max()),
                    "p95_log_gain": float(np.quantile(gap[:, j], 0.95))}
            for j, model in enumerate(MODELS)}


def synchronization_check():
    """U is a diagnostic upper surrogate, never a true cube maximal value."""
    rng = np.random.default_rng(SEED + 1)
    result = {"seed": SEED + 1, "samples_per_dimension": 1024, "dimensions": {}}
    for n in [32, 64]:
        x = rng.random((1024, n))
        shifted = (x[:, :, None] - PHASES[None, None, :]) % 1
        cosine = np.cos(2 * np.pi * shifted)
        phase_maxima = np.full((len(x), len(PHASES)), -np.inf)
        shared_max = np.full(len(x), -np.inf)
        for length in LS:
            phase_products = np.log(average_from_shifted(
                shifted, length, cosine)).sum(axis=1)
            phase_maxima = np.maximum(phase_maxima, phase_products)
            shared_max = np.maximum(shared_max,
                                     lse(phase_products) - math.log(len(PHASES)))
        log_u = lse(phase_maxima) - math.log(len(PHASES))
        gap = log_u - shared_max
        assert gap.min() >= -1e-12
        assert gap.max() <= math.log(len(PHASES)) + 1e-12
        result["dimensions"][str(n)] = {
            "log_U_over_V_quantiles_0_50_90_99_100":
                np.quantile(gap, [0, .5, .9, .99, 1]).tolist(),
            "mean_log_U_over_V": float(gap.mean()),
            "fixed_thresholds": [
                {"threshold": t, "P_U_above": float(np.mean(log_u > math.log(t))),
                 "P_V_above": float(np.mean(shared_max > math.log(t))),
                 "P_U_above_V_not": float(np.mean((log_u > math.log(t))
                                                   & (shared_max <= math.log(t))))}
                for t in [0.75, 1, 1.5, 2, 4]],
        }
    return result


def main():
    assert (1 - P) - abs(COS_AMP) > 0
    assert np.all(ATOMS - DELTA / 2 > 0)
    assert np.all(ATOMS + DELTA / 2 < 1)
    derived = np.sort(((ATOMS[:, None] + PHASES[None, :]) % 1).ravel())
    assert np.allclose(derived, np.arange(12) / 12)
    rng = np.random.default_rng(SEED)
    meta = {"seed": SEED, "python_executable": sys.executable,
            "python_version": sys.version, "numpy_version": np.__version__,
            "K": K, "p": P, "delta": DELTA, "atoms": ATOMS.tolist(),
            "phases": PHASES.tolist(), "b1": B1, "c1": C1,
            "cos_density_amplitude": COS_AMP,
            "background_lower_bound": (1 - P) - abs(COS_AMP),
            "base_density_upper_bound": (1 - P) + P / (K * DELTA)
                                        + abs(COS_AMP),
            "scale_interval": [float(LS[0]), float(LS[-1])],
            "scale_grid_size": len(LS), "proposal_anchors": ANCHORS.tolist()}
    (OUT / "parameters.json").write_text(json.dumps(meta, indent=2))
    direct = direct_cube_check(rng)
    (OUT / "direct_cube_check.json").write_text(json.dumps(direct, indent=2))
    pilot = run_case(rng, 16, 2048, LS[::3], "uniform", "pilot_n16")
    (OUT / "pilot.json").write_text(json.dumps(pilot, indent=2))
    results = []
    for dimension in [8, 32, 64]:
        results.append(run_case(rng, dimension, 8192, LS,
                                "half_uniform_half_product_tilts",
                                f"main_n{dimension}"))
    (OUT / "results.json").write_text(json.dumps(results, indent=2))
    validations = {}
    for dimension in [8, 32, 64]:
        data = np.load(OUT / f"main_n{dimension}_samples.npz")
        validations[str(dimension)] = fixed_validation(data["log_max"],
                                                       data["weights"])
    (OUT / "fixed_threshold_validation.json").write_text(
        json.dumps(validations, indent=2))
    resolution = resolution_check(rng)
    (OUT / "resolution_check.json").write_text(json.dumps(resolution, indent=2))
    (OUT / "synchronization_check.json").write_text(
        json.dumps(synchronization_check(), indent=2))
    print("DONE", flush=True)


if __name__ == "__main__":
    main()
