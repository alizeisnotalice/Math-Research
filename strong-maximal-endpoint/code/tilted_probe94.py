#!/usr/bin/env python3
"""Deterministic tilted Gamma/Binomial CDF probes for r094.

This is a numerical diagnostic only. It uses finite positive mixture sums with
explicit numerical tail budgets; SciPy double precision values are not directed
intervals and do not certify a physical halo norm.
"""
from __future__ import annotations

import concurrent.futures as cf
import argparse
import hashlib
import json
import math
import multiprocessing as mp
import os
import platform
import sys
import time
from pathlib import Path

for _thread_env in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ[_thread_env] = "1"

import numpy as np
import scipy
from scipy.optimize import brentq
from scipy.special import gammainc, ndtr
from scipy.stats import binom

OUT_DIR = Path(sys.argv[1]).resolve() if len(sys.argv)>1 else Path("/tmp/r094-tilted-replay")
NS = (64, 128, 256, 512, 1024, 4096, 16384, 65536)
TOLS = (1e-12, 1e-15)
XIS = (-4.0, -2.0, -1.0, 0.0, 1.0, 2.0, 4.0)
GAMMA_CHUNK = 32
TAIL_CAP = 1e-14
WORKERS = 12
CASE_TIMEOUT_S = 120.0


def mu_root() -> float:
    return float(brentq(lambda x: x - 1.0 - math.exp(-x), 1.0, 2.0,
                        xtol=5e-15, rtol=4 * np.finfo(float).eps))


MU = mu_root()
SIGMA = math.sqrt(MU)
C_Q = (9.0 - MU) / 24.0


def action_residual(u: float, c: float) -> float:
    if u == 0.0:
        return c - 1.0 - math.exp(-c)
    return math.expm1(u * c) / u - 1.0 - math.exp(-c)


def action_c(u: float) -> tuple[float, float]:
    if u <= 0.0:
        return MU, action_residual(0.0, MU)
    hi = 1.0
    while action_residual(u, hi) <= 0.0:
        hi *= 2.0
        if hi > 128.0:
            raise ArithmeticError(f"action root bracket failed for u={u}")
    c = float(brentq(lambda x: action_residual(u, x), 0.0, hi,
                     xtol=5e-15, rtol=4 * np.finfo(float).eps))
    return c, action_residual(u, c)


def mean_tilt(u: float) -> tuple[float, float, float]:
    c, residual = action_c(u)
    m = (1.0 + c / (1.0 + math.exp(c))) / (1.0 + u)
    return m, c, residual


def z_root(n: int) -> tuple[float, float]:
    L = math.log(n)
    hi = max(3.0, math.sqrt(L) + 3.0)
    f = lambda z: float(ndtr(-z)) - SIGMA * z / math.sqrt(n)
    z = float(brentq(f, 0.0, hi, xtol=5e-15,
                     rtol=4 * np.finfo(float).eps))
    return z, f(z)


def scalar_data(n: int, q_offset: float = 0.0, q_shift_label: str = "base") -> dict:
    L = math.log(n)
    z, z_resid = z_root(n)
    tau = SIGMA * z / math.sqrt(n)
    d = n * tau
    q0 = MU * n - d
    q = q0 + C_Q * L + q_offset
    if not (q > 0.0 and 0.0 < tau < 1.0):
        raise ArithmeticError(f"invalid q/tau at N={n}: q={q}, tau={tau}")

    def ref_fun(s):
        m, _, _ = mean_tilt(s)
        return n * m - q

    s_hi = max(0.5, math.sqrt(L / n))
    while ref_fun(s_hi) > 0.0:
        s_hi *= 2.0
    s0 = float(brentq(ref_fun, 0.0, s_hi, xtol=5e-15,
                      rtol=4 * np.finfo(float).eps))
    ref_mean_resid = ref_fun(s0)
    c0, c0_resid = action_c(s0)
    return {
        "N": n, "L": L, "mu": MU, "sigma": SIGMA, "c_Q": C_Q,
        "z": z, "z_root_residual": z_resid, "tau": tau, "d": d,
        "q0": q0, "q": q, "q_offset": q_offset,
        "q_shift_label": q_shift_label, "s0": s0,
        "reference_mean_residual": ref_mean_resid,
        "c0": c0, "c0_root_residual": c0_resid,
    }


def root_cutoff_binomial(nn: int, p: float, target: float) -> dict:
    """Central binomial support with each omitted tail <= target (numerically)."""
    if nn == 0:
        return {"lo": 0, "hi": 0, "ks": np.array([0], dtype=np.int64),
                "pmf": np.array([1.0]), "pmfsum": 1.0,
                "left_tail": 0.0, "right_tail": 0.0,
                "target_each_tail": target}
    lo = int(binom.ppf(target, nn, p))
    hi = int(binom.isf(target, nn, p))
    lo = max(0, min(nn, lo))
    hi = max(lo, min(nn, hi))
    # Guard the inverse-CDF convention/rounding using evaluated CDF/SF values.
    while lo > 0 and float(binom.cdf(lo - 2, nn, p)) >= target:
        lo -= 1
    while lo < nn and float(binom.cdf(lo - 1, nn, p)) >= target:
        lo += 1
    while hi < nn and float(binom.sf(hi, nn, p)) > target:
        hi += 1
    while hi > 0 and float(binom.sf(hi - 1, nn, p)) <= target:
        hi -= 1
    if lo > hi:
        raise ArithmeticError(f"empty binomial support nn={nn}, p={p}, target={target}")
    ks = np.arange(lo, hi + 1, dtype=np.int64)
    pmf = np.asarray(binom.pmf(ks, nn, p), dtype=float)
    pmfsum = math.fsum(float(x) for x in pmf)
    left_tail = float(binom.cdf(lo - 1, nn, p)) if lo else 0.0
    right_tail = float(binom.sf(hi, nn, p)) if hi < nn else 0.0
    if left_tail > target * (1.0 + 1e-5) or right_tail > target * (1.0 + 1e-5):
        raise ArithmeticError(f"binomial tail budget exceeded: {left_tail}, {right_tail}, target={target}")
    return {"lo": lo, "hi": hi, "ks": ks, "pmf": pmf,
            "pmfsum": pmfsum, "left_tail": left_tail,
            "right_tail": right_tail, "target_each_tail": target}


def geometric_support(rho: float, target: float) -> dict:
    if rho >= 1.0 - 4 * np.finfo(float).eps:
        return {"k": np.array([0], dtype=np.int64), "pmf": np.array([1.0]),
                "pmfsum": 1.0, "tail_bound": 0.0, "Kmax": 0,
                "rho": 1.0}
    if not (0.0 < rho < 1.0):
        raise ArithmeticError(f"invalid geometric success probability {rho}")
    logtail = math.log1p(-rho)
    Kmax = max(0, int(math.ceil(math.log(target) / logtail)) - 1)
    kk = np.arange(Kmax + 1, dtype=np.int64)
    pmf = rho * np.exp(kk.astype(float) * logtail)
    pmfsum = math.fsum(float(x) for x in pmf)
    tail = math.exp((Kmax + 1) * logtail)
    if tail > target * (1.0 + 1e-12):
        raise ArithmeticError(f"geometric tail budget exceeded: {tail}>{target}")
    return {"k": kk, "pmf": pmf, "pmfsum": pmfsum,
            "tail_bound": tail, "Kmax": Kmax, "rho": rho}


def reference_cdf(n: int, q: float, c0: float, tail_target: float) -> tuple[float, dict]:
    p0 = 1.0 / (1.0 + math.exp(c0))
    support = root_cutoff_binomial(n, p0, tail_target)
    x = np.maximum(q - c0 * support["ks"].astype(float), 0.0)
    vals = gammainc(float(n), x)
    cdf = math.fsum(float(a * b) for a, b in zip(support["pmf"], vals))
    return cdf, {
        "binomial_n": n, "binomial_p": p0,
        "binomial_lo": support["lo"], "binomial_hi": support["hi"],
        "binomial_pmfsum": support["pmfsum"],
        "binomial_left_tail": support["left_tail"],
        "binomial_right_tail": support["right_tail"],
        "binomial_target_each_tail": support["target_each_tail"],
        "geometric_pmfsum": 1.0, "geometric_tail_bound": 0.0,
        "combined_tail_bound": support["left_tail"] + support["right_tail"],
    }


def weighted_cdf(n: int, q: float, a: float, b: float,
                 c1: float, c2: float, tail_target: float,
                 deadline: float | None = None) -> tuple[float, dict]:
    p1 = 1.0 / (1.0 + math.exp(c1))
    p2 = 1.0 / (1.0 + math.exp(c2))
    bsup = root_cutoff_binomial(n - 1, p2, tail_target)
    if abs(a - b) <= 8 * np.finfo(float).eps * max(a, b):
        gsup = {"k": np.array([0], dtype=np.int64), "pmf": np.array([1.0]),
                "pmfsum": 1.0, "tail_bound": 0.0, "Kmax": 0, "rho": 1.0}
        shapes = np.array([n], dtype=float)
    else:
        rho = b / a
        gsup = geometric_support(rho, tail_target)
        shapes = n + gsup["k"].astype(float)

    branch_cdfs = []
    for foreground_bit in (0, 1):
        shift0 = a * c1 * foreground_bit
        branch_total = 0.0
        for start in range(0, len(shapes), GAMMA_CHUNK):
            if deadline is not None and time.monotonic() > deadline:
                raise TimeoutError("case exceeded 120-second wall limit during Gamma chunks")
            shape_block = shapes[start:start + GAMMA_CHUNK]
            geom_block = gsup["pmf"][start:start + GAMMA_CHUNK]
            x = np.maximum(
                (q - shift0 - b * c2 * bsup["ks"].astype(float)) / b,
                0.0,
            )
            gamma_block = gammainc(shape_block[:, None], x[None, :])
            bin_average = gamma_block @ bsup["pmf"]
            branch_total += float(np.dot(geom_block, bin_average))
        branch_cdfs.append(branch_total)
    cdf = (1.0 - p1) * branch_cdfs[0] + p1 * branch_cdfs[1]
    return cdf, {
        "foreground_p": p1, "background_p": p2,
        "binomial_n": n - 1, "binomial_lo": bsup["lo"],
        "binomial_hi": bsup["hi"], "binomial_pmfsum": bsup["pmfsum"],
        "binomial_left_tail": bsup["left_tail"],
        "binomial_right_tail": bsup["right_tail"],
        "binomial_target_each_tail": bsup["target_each_tail"],
        "geometric_success_probability": gsup["rho"],
        "geometric_Kmax": gsup["Kmax"],
        "geometric_pmfsum": gsup["pmfsum"],
        "geometric_tail_bound": gsup["tail_bound"],
        "combined_tail_bound": gsup["tail_bound"] + bsup["left_tail"] + bsup["right_tail"],
    }


def candidate_job(spec: tuple) -> dict:
    if len(spec) == 3:
        n, a, xi = spec
        q_offset, q_shift_label = 0.0, "base"
    elif len(spec) == 5:
        n, a, xi, q_offset, q_shift_label = spec
    else:
        raise ValueError(f"unsupported case spec: {spec}")
    n, a, xi = int(n), float(a), float(xi)
    q_offset, q_shift_label = float(q_offset), str(q_shift_label)
    deadline = time.monotonic() + CASE_TIMEOUT_S
    scal = scalar_data(n, q_offset, q_shift_label)
    b = (n - a) / (n - 1.0)
    if not (a > 0.0 and b > 0.0 and a < n):
        raise ArithmeticError(f"invalid normalized weights N={n}, a={a}, b={b}")

    def joint_mean(s: float) -> float:
        m1, _, _ = mean_tilt(s * a)
        mb, _, _ = mean_tilt(s * b)
        return a * m1 + (n - 1.0) * b * mb

    if a == 1.0:
        s = scal["s0"]
    else:
        hi = max(0.5, scal["s0"])
        while joint_mean(hi) > scal["q"]:
            hi *= 2.0
            if hi > 1e4:
                raise ArithmeticError(f"joint tilt root bracket failed N={n},a={a}")
        s = float(brentq(lambda x: joint_mean(x) - scal["q"], 0.0, hi,
                         xtol=5e-15, rtol=4 * np.finfo(float).eps))
    joint_mean_resid = joint_mean(s) - scal["q"]
    c1, c1_resid = action_c(s * a)
    c2base, c2base_resid = action_c(s * b)
    offset = xi / math.sqrt(n * scal["L"])
    c2 = c2base + offset
    c2_resid = action_residual(s * b, c2base)
    if c2 <= 0.0:
        raise ArithmeticError(f"perturbed action c2 nonpositive at {spec}: {c2}")

    exact_identical = (a == 1.0 and xi == 0.0)
    tol_rows = []
    for tol in TOLS:
        tail_target = min(tol, TAIL_CAP)
        Fref, ref_status = reference_cdf(n, scal["q"], scal["c0"], tail_target)
        if exact_identical:
            # With a=b=1 and c1=c2=c0, B1+Bin(N-1,p0) is Bin(N,p0).
            Fw = Fref
            weighted_status = dict(ref_status)
            weighted_status["cdf_evaluation"] = "exact iid distributional identity; reference summation path"
        else:
            Fw, weighted_status = weighted_cdf(
                n, scal["q"], a, b, c1, c2, tail_target, deadline
            )
            weighted_status["cdf_evaluation"] = "positive Gamma/Geometric/Binomial mixture; shape chunk 32"
        if not (math.isfinite(Fw) and math.isfinite(Fref) and 0.0 < Fw < 1.0 and 0.0 < Fref < 1.0):
            raise ArithmeticError(f"CDF outside finite open unit interval: spec={spec}, Fw={Fw}, Fref={Fref}")
        if weighted_status["combined_tail_bound"] > 3e-14 * (1.0 + 1e-5):
            raise ArithmeticError(f"weighted combined tail budget exceeded: {weighted_status}")
        if ref_status["combined_tail_bound"] > 2e-14 * (1.0 + 1e-5):
            raise ArithmeticError(f"reference tail budget exceeded: {ref_status}")
        lr = math.log(Fw) - math.log(Fref)
        tol_rows.append({
            "tol": tol, "tail_target_each_component": tail_target,
            "F_weighted": Fw, "F_reference": Fref,
            "logratio": lr, "sqrtN_logratio": math.sqrt(n) * lr,
            "sqrtNL_logratio": math.sqrt(n * scal["L"]) * lr,
            "weighted_status": weighted_status,
            "reference_status": ref_status,
        })

    if exact_identical:
        if s != scal["s0"] or c1 != scal["c0"] or c2 != scal["c0"]:
            raise AssertionError(f"identity parameters did not match exactly: {spec}")
        if any(r["F_weighted"] != r["F_reference"] or r["logratio"] != 0.0 for r in tol_rows):
            raise AssertionError(f"a=1,xi=0 CDF identity failed: {spec}")

    return {
        **scal,
        "a": a, "b": b, "lambda1": a, "lambda_bulk": b,
        "xi": xi, "s": s,
        "joint_mean_residual": joint_mean_resid,
        "c1": c1, "c1_root_residual": c1_resid,
        "c2_base": c2base, "c2_base_root_residual": c2base_resid,
        "c2_perturbation": offset, "c2": c2,
        "c2_unperturbed_root_residual": c2_resid,
        "exact_iid_identity_case": exact_identical,
        "tolerance_logratio_difference": abs(tol_rows[0]["logratio"] - tol_rows[1]["logratio"]),
        "tol_rows": tol_rows,
        "pid": os.getpid(),
    }


def unique_a_values(n: int) -> list[float]:
    vals = {1.0, 2.0, 4.0, 8.0, float(n) ** 0.25, 2.0 * float(n) ** 0.25}
    return sorted(x for x in vals if x < n)


def all_specs() -> list[tuple[int, float, float]]:
    return [(n, a, xi) for n in NS for a in unique_a_values(n) for xi in XIS]


def validate_case(r: dict) -> None:
    for key in ("q", "s0", "s", "c0", "c1", "c2", "joint_mean_residual"):
        if not math.isfinite(float(r[key])):
            raise ArithmeticError(f"nonfinite case field {key}: {r}")
    if not (r["s"] > 0.0 and r["s0"] > 0.0 and len(r["tol_rows"]) == 2):
        raise ArithmeticError(f"invalid roots/tolerance rows: {r}")
    for tr in r["tol_rows"]:
        for k in ("F_weighted", "F_reference", "logratio", "sqrtN_logratio", "sqrtNL_logratio"):
            if not math.isfinite(float(tr[k])):
                raise ArithmeticError(f"nonfinite {k}: {r}")
        if not (0.0 < tr["F_weighted"] < 1.0 and 0.0 < tr["F_reference"] < 1.0):
            raise ArithmeticError(f"CDF domain failure: {r}")
        if tr["weighted_status"]["combined_tail_bound"] > 3e-14 * (1.0 + 1e-5):
            raise ArithmeticError(f"weighted tail status failed: {r}")


def file_sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def main() -> None:
    specs = all_specs()
    pilot_n = max(NS)
    pilot_specs = [s for s in specs if s[0] == pilot_n]
    if len(pilot_specs) != 42:
        raise AssertionError(f"expected 42 maximum-N pilot cases; got {len(pilot_specs)}")

    started = time.monotonic()
    pids: set[int] = set()
    completed: dict[tuple[int, float, float], dict] = {}
    ctx = mp.get_context("spawn")
    with cf.ProcessPoolExecutor(max_workers=WORKERS, mp_context=ctx) as pool:
        for i, result in enumerate(pool.map(candidate_job, pilot_specs, chunksize=1), start=1):
            validate_case(result)
            pids.add(int(result["pid"]))
            spec = (int(result["N"]), float(result["a"]), float(result["xi"]))
            completed[spec] = result
            if i % 12 == 0:
                print(json.dumps({"phase": "pilot", "completed": i,
                                  "total": len(pilot_specs), "workers": len(pids),
                                  "seconds": time.monotonic() - started}), flush=True)

        if len(completed) != len(pilot_specs) or len(pids) != WORKERS:
            raise RuntimeError(f"pilot incomplete: cases={len(completed)}, workers={len(pids)}")
        pilot_identity = [r for r in completed.values() if r["exact_iid_identity_case"]]
        if len(pilot_identity) != 1 or pilot_identity[0]["N"] != pilot_n:
            raise AssertionError("maximum-N pilot lacks its a=1,xi=0 identity case")
        print(json.dumps({"phase": "pilot_passed", "cases": len(completed),
                          "workers": len(pids),
                          "max_two_tolerance_logratio_difference": max(r["tolerance_logratio_difference"] for r in completed.values())}, indent=2), flush=True)

        rest = [s for s in specs if s not in completed]
        for i, result in enumerate(pool.map(candidate_job, rest, chunksize=1), start=1):
            validate_case(result)
            pids.add(int(result["pid"]))
            spec = (int(result["N"]), float(result["a"]), float(result["xi"]))
            if spec in completed:
                raise AssertionError(f"duplicate case {spec}")
            completed[spec] = result
            if i % 24 == 0:
                print(json.dumps({"phase": "full", "new_completed": i,
                                  "new_total": len(rest), "workers": len(pids),
                                  "seconds": time.monotonic() - started}), flush=True)
    elapsed = time.monotonic() - started

    rows = [completed[s] for s in specs]
    if len(rows) != len(specs) or len(pids) != WORKERS:
        raise RuntimeError(f"full batch incomplete: rows={len(rows)}, specs={len(specs)}, workers={len(pids)}")
    for n in NS:
        for a in unique_a_values(n):
            case = completed[(n, a, 0.0)]
            if a == 1.0 and not case["exact_iid_identity_case"]:
                raise AssertionError("a=1,xi=0 exact identity missing")

    OUTPUT_DIR = OUT_DIR
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    rows_path = OUTPUT_DIR / "rows.json"
    rows_path.write_text(json.dumps(rows, indent=2, allow_nan=False) + "\n")
    rel_diffs = [r["tolerance_logratio_difference"] for r in rows]
    tail_vals = [tr["weighted_status"]["combined_tail_bound"]
                 for r in rows for tr in r["tol_rows"]]
    metric_records = {
        key: [
            {"value": float(tr[key]), "N": r["N"], "a": r["a"], "b": r["b"],
             "xi": r["xi"], "tol": tr["tol"], "s": r["s"],
             "c1": r["c1"], "c2": r["c2"]}
            for r in rows for tr in r["tol_rows"]
        ]
        for key in ("logratio", "sqrtN_logratio", "sqrtNL_logratio")
    }

    def metric_extremes(records):
        return {"min": min(records, key=lambda r: r["value"]),
                "max": max(records, key=lambda r: r["value"])}

    sqrt_nl_records = metric_records["sqrtNL_logratio"]
    perturbation_groups = {
        "xi_zero": [x for x in sqrt_nl_records if x["xi"] == 0.0],
        "xi_nonzero": [x for x in sqrt_nl_records if x["xi"] != 0.0],
        "a_equals_one_xi_nonzero": [x for x in sqrt_nl_records if x["a"] == 1.0 and x["xi"] != 0.0],
        "a_greater_than_one": [x for x in sqrt_nl_records if x["a"] > 1.0],
    }
    summary = {
        "pass_all": True,
        "scope": "Deterministic double-precision tilted Gamma/Geometric/Binomial CDF comparison only; finite mixture tail budgets are numerical estimates, not directed intervals; no physical halo norm or onset dimension is certified.",
        "actual_workers": len(pids), "pids": sorted(pids),
        "target_precision": "SciPy float64",
        "scipy_version": scipy.__version__, "numpy_version": np.__version__,
        "python_version": platform.python_version(),
        "N_values": list(NS), "xi_values": list(XIS), "tol_values": list(TOLS),
        "chunk_gamma_shapes": GAMMA_CHUNK,
        "pilot_cases": len(pilot_specs), "pilot_N": pilot_n,
        "total_cases": len(rows), "expected_total_cases": len(specs),
        "cases_by_N": {str(n): len(unique_a_values(n)) * len(XIS) for n in NS},
        "exact_identity_cases": sum(bool(r["exact_iid_identity_case"]) for r in rows),
        "all_CDFs_finite_open_unit_interval": True,
        "all_candidate_component_tail_budgets_at_most_1e-14": max(
            max(tr["weighted_status"]["geometric_tail_bound"],
                tr["weighted_status"]["binomial_left_tail"],
                tr["weighted_status"]["binomial_right_tail"])
            for r in rows for tr in r["tol_rows"]
        ) <= 1e-14 * (1.0 + 1e-5),
        "all_candidate_combined_tail_budgets_at_most_3e-14": max(tail_vals) <= 3e-14 * (1.0 + 1e-5),
        "max_weighted_combined_tail_bound": max(tail_vals),
        "max_two_tolerance_logratio_difference": max(rel_diffs),
        "min_two_tolerance_logratio_difference": min(rel_diffs),
        "CDF_logratio_extremes": {k: metric_extremes(v) for k, v in metric_records.items()},
        "sqrtNL_logratio_extremes_by_perturbation": {
            k: metric_extremes(v) for k, v in perturbation_groups.items()
        },
        "max_abs_joint_mean_residual_over_N": max(abs(r["joint_mean_residual"]) / r["N"] for r in rows),
        "max_abs_reference_mean_residual_over_N": max(abs(r["reference_mean_residual"]) / r["N"] for r in rows),
        "wall_seconds": elapsed,
        "resource_scope": "Launched under math_guard with reserve 1 GiB and owned RSS budget 6.5 GiB; guard report is a sampled macOS RSS/available-memory estimate, not a hard OS reservation.",
        "source_sha256": {
            str(Path(__file__).resolve()): file_sha256(Path(__file__).resolve()),
            "rounds/r094/core94.md": file_sha256(Path(__file__).resolve().parent.parent / "rounds/r094/core94.md"),
            "rounds/r094/audit94.md": file_sha256(Path(__file__).resolve().parent.parent / "rounds/r094/audit94.md"),
        },
        "rows_sha256": file_sha256(rows_path),
    }
    summary_path = OUTPUT_DIR / "summary.json"
    summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True, allow_nan=False) + "\n")
    print(json.dumps({"phase": "complete", "summary": str(summary_path),
                      "rows": str(rows_path), "cases": len(rows),
                      "workers": len(pids), "wall_seconds": elapsed,
                      "max_tail": max(tail_vals), "max_tol_diff": max(rel_diffs)}, indent=2), flush=True)


if __name__ == "__main__":
    main()
