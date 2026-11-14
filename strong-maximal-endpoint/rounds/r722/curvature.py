#!/usr/bin/env python3
"""High precision finite-angle Hessian diagnostics for the uniform baseline.

All outputs are numerical, not interval-certified. Baseline CDF sums use the
finite Gamma-binomial mixture; angular transforms use the stated hypergeometric
closed forms. Two precision passes (40 and 70 dps) are retained for comparison.
"""
from __future__ import annotations

import itertools
import json
import math
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
NS = (4, 8, 16, 32, 64, 128)
AMPS = ("0.1", "0.3", "1", "3")
PRECISIONS = (40, 70)
WORKERS = 12
MAXTERMS = 20000
CASE_BUDGET_SECONDS = 300


def gamma_pdf(N, x):
    if x <= 0:
        return mp.mpf(0)
    return mp.exp((N - 1) * mp.log(x) - x - mp.loggamma(N))


def gamma_cdf(N, x):
    if x <= 0:
        return mp.mpf(0)
    return mp.gammainc(N, 0, x) / mp.gamma(N)


def binom_weight(N, j, p):
    return mp.binomial(N, j) * p**j * (1 - p)**(N - j)


def mixture_pdf(N, s, c, p):
    return mp.fsum(binom_weight(N, j, p) * gamma_pdf(N, s - j * c)
                   for j in range(N + 1))


def mixture_cdf(N, s, c, p):
    return mp.fsum(binom_weight(N, j, p) * gamma_cdf(N, s - j * c)
                   for j in range(N + 1))


def baseline(n, dps):
    mp.mp.dps = dps
    mp.mp.maxterms = MAXTERMS

    def p_of(c):
        return 1 / (1 + mp.exp(c))

    # Solve the two median equations simultaneously. The initial scale follows
    # the mean-one Bernoulli-exponential increment fixed point c=1+c p(c).
    c0 = mp.mpf("1.28")
    r0 = n * (1 + c0 * p_of(c0))

    def equations(r, c):
        p = p_of(c)
        return (mixture_cdf(n, r, c, p) - mp.mpf("0.5"),
                mixture_cdf(n - 1, r - c, c, p) - mp.mpf("0.5"))

    r, c = mp.findroot(equations, (r0, c0), tol=mp.mpf(10) ** (-(dps - 8)),
                       maxsteps=80, solver="mdnewton")
    p = p_of(c)
    residuals = equations(r, c)
    return {"n": n, "dps": dps, "r": r, "c": c, "p": p,
            "cdf_n_residual": residuals[0],
            "cdf_nminus1_shifted_residual": residuals[1]}


def poly_coeffs(a, b, z, phase):
    coeff = [mp.mpc(1)]
    for value, multiplicity in ((phase, a), (1 / phase, b), (mp.mpc(1), z)):
        for _ in range(multiplicity):
            updated = [mp.mpc(0)] * (len(coeff) + 1)
            for j, old in enumerate(coeff):
                updated[j] += old
                updated[j + 1] += old * value
            coeff = updated
    return coeff


def angular_chi(a, b, z, N, t, x, mode):
    if x <= 0:
        return mp.mpc(0)
    itx = 1j * t * x
    if mode == "dense":
        # No zero coordinates: the Dirichlet sum of the positive coordinates
        # is Beta(a,b); this includes the balanced 0F1 specialization.
        return mp.exp(-itx) * mp.hyp1f1(a, N, 2 * itx)
    if mode == "sparse" or mode == "sparse_remove_zero":
        return mp.hyper([1], [mp.mpf(N) / 2, mp.mpf(N + 1) / 2],
                        -(t * x) ** 2 / 4)
    if mode == "sparse_remove_plus":
        return mp.hyp1f1(1, N, -itx)
    raise ValueError(f"unknown angular transform mode {mode}")


def weighted_density(N, s, c, p, a, b, z, t, mode):
    phase = mp.exp(1j * t * c)
    coeff = poly_coeffs(a, b, z, phase)
    terms = []
    for j, ej in enumerate(coeff):
        x = s - j * c
        if x <= 0 or ej == 0:
            continue
        chi = angular_chi(a, b, z, N, t, x, mode)
        terms.append(p**j * (1 - p)**(N - j) * ej * chi * gamma_pdf(N, x))
    return mp.fsum(terms)


def fmt(x, dps):
    if isinstance(x, complex):
        x = mp.mpc(x)
    return mp.nstr(x, n=dps)


def evaluate_case(task):
    n, family, ampstr, dps, base_strings = task
    tstart = time.monotonic()
    mp.mp.dps = dps
    mp.mp.maxterms = MAXTERMS
    try:
        r = mp.mpf(base_strings["r"])
        c = mp.mpf(base_strings["c"])
        p = mp.mpf(base_strings["p"])
        amp = mp.mpf(ampstr)
        if family == "sparse":
            a, b, z, angle, mode = 1, 1, n - 2, amp, "sparse"
        elif family == "dense":
            a, b, z, angle, mode = n // 2, n // 2, 0, amp, "dense"
        elif family == "dense_energy":
            a, b, z, angle, mode = n // 2, n // 2, 0, amp / mp.sqrt(n), "dense"
        else:
            raise ValueError(f"unknown family {family}")

        f_n_r = mixture_pdf(n, r, c, p)
        f_nm1_r = mixture_pdf(n - 1, r, c, p)
        f_nm1_shift = mixture_pdf(n - 1, r - c, c, p)
        F_n2_shift = mixture_cdf(n - 2, r - 2 * c, c, p)
        if min(f_n_r, f_nm1_r, f_nm1_shift) <= 0:
            raise ArithmeticError("nonpositive baseline density")

        Ffull = weighted_density(n, r, c, p, a, b, z, angle, mode)
        m = Ffull / f_n_r
        if family in ("sparse",):
            F_removed_plus = weighted_density(
                n - 1, r, c, p, 0, 1, n - 2, angle, "sparse_remove_plus")
            F_removed_plus_shift = weighted_density(
                n - 1, r - c, c, p, 0, 1, n - 2, angle,
                "sparse_remove_plus")
            F_removed_zero_shift = weighted_density(
                n - 1, r - c, c, p, 1, 1, n - 3, angle,
                "sparse_remove_zero")
        else:
            F_removed_plus = weighted_density(
                n - 1, r, c, p, a - 1, b, z, angle, "dense")
            F_removed_plus_shift = weighted_density(
                n - 1, r - c, c, p, a - 1, b, z, angle, "dense")
            F_removed_zero_shift = mp.mpc(0)

        m_plus0 = F_removed_plus / f_nm1_r
        m_plusc = mp.exp(1j * angle * c) * F_removed_plus_shift / f_nm1_shift
        if z:
            m_zero = F_removed_zero_shift / f_nm1_shift
        else:
            m_zero = None

        h = p * f_nm1_shift / f_n_r
        offb = p**2 * (F_n2_shift - mp.mpf("0.5")) / f_n_r
        D_over_f = ((1 - p) * f_nm1_r / f_n_r *
                    (1 + abs(m)**2 - 2 * mp.re(mp.conj(m) * m_plus0)) +
                    h * (1 + abs(m)**2 - 2 * mp.re(mp.conj(m) * m_plusc)))
        gplus = h * (m - m_plusc)
        gzero = h * (m - m_zero) if m_zero is not None else mp.mpc(0)
        schur = ((2 * a * abs(gplus)**2 + z * abs(gzero)**2) / (h - offb) -
                 offb * abs(2 * a * mp.re(gplus) + z * gzero)**2 /
                 ((h - offb) * (h + (n - 1) * offb)))
        C_all = (schur - D_over_f) / 4
        C_fixed = -D_over_f / 4
        C_singleaxis = (abs(gplus)**2 / h - D_over_f) / 4
        # The Hessian's diagonal term is h (not D/f); D/f belongs to the
        # Schur-complement curvature formula above.
        H_transverse = h - offb
        H_longitudinal = h + (n - 1) * offb
        stationarity = {
            "baseline_cdf_n_minus_half": base_strings["cdf_n_residual"],
            "baseline_cdf_nminus1_shifted_minus_half": base_strings["cdf_nminus1_shifted_residual"],
            "symmetric_angle_imaginary_m_residual": fmt(mp.im(m), dps),
        }
        elapsed = time.monotonic() - tstart
        return {"n": n, "family": family, "amplitude": ampstr,
                "angle_t": fmt(angle, dps), "precision_dps": dps,
                "status": "ok" if elapsed <= CASE_BUDGET_SECONDS else "timeout_budget_exceeded",
                "elapsed_seconds": elapsed, "hypergeometric_maxterms": MAXTERMS,
                "counts_plus_minus_zero": [a, b, z],
                "baseline": {"r": fmt(r, dps), "c": fmt(c, dps), "p": fmt(p, dps),
                             "f_n_r": fmt(f_n_r, dps), "h": fmt(h, dps),
                             "off_diagonal": fmt(offb, dps)},
                "stationarity_residuals": stationarity,
                "moments": {"m": fmt(m, dps), "imag_m": fmt(mp.im(m), dps),
                            "m_plus0": fmt(m_plus0, dps), "m_plusc": fmt(m_plusc, dps),
                            "m_zero": fmt(m_zero, dps) if m_zero is not None else None},
                "D_over_f": fmt(D_over_f, dps), "Schur": fmt(schur, dps),
                "C_all": fmt(C_all, dps), "C_fixed": fmt(C_fixed, dps),
                "C_singleaxis": fmt(C_singleaxis, dps),
                "H_eigenvalues": {"transverse_multiplicity_n_minus_1": fmt(H_transverse, dps),
                                  "longitudinal_multiplicity_1": fmt(H_longitudinal, dps)},
                "H_eigen_signs": {"transverse": signof(H_transverse),
                                  "longitudinal": signof(H_longitudinal)},
                "curvature_signs": {"C_all": signof(C_all), "C_fixed": signof(C_fixed),
                                    "C_singleaxis": signof(C_singleaxis)}}
    except Exception as exc:
        return {"n": n, "family": family, "amplitude": ampstr,
                "precision_dps": dps, "status": "error", "elapsed_seconds": time.monotonic()-tstart,
                "error": f"{type(exc).__name__}: {exc}",
                "traceback": traceback.format_exc(limit=3)}


def signof(x):
    if x > 0:
        return "positive"
    if x < 0:
        return "negative"
    return "zero"


def main():
    start = time.time()
    baselines = {}
    baseline_records = []
    for n in NS:
        for dps in PRECISIONS:
            try:
                b = baseline(n, dps)
                encoded = {k: fmt(v, dps) if isinstance(v, mp.mpf) else v
                           for k, v in b.items()}
                baselines[(n, dps)] = {k: encoded[k] for k in (
                    "r", "c", "p", "cdf_n_residual", "cdf_nminus1_shifted_residual")}
                baseline_records.append({"n": n, "precision_dps": dps,
                                         "status": "ok", **encoded})
            except Exception as exc:
                baseline_records.append({"n": n, "precision_dps": dps,
                                         "status": "error", "error": f"{type(exc).__name__}: {exc}",
                                         "traceback": traceback.format_exc(limit=3)})

    tasks = []
    for n in NS:
        for family in ("sparse", "dense", "dense_energy"):
            for amp in AMPS:
                for dps in PRECISIONS:
                    b = baselines.get((n, dps))
                    if b is None:
                        tasks.append((n, family, amp, dps, None))
                    else:
                        tasks.append((n, family, amp, dps, b))
    results = []
    with ProcessPoolExecutor(max_workers=WORKERS) as pool:
        futures = {}
        for n, family, amp, dps, base in tasks:
            if base is None:
                results.append({"n": n, "family": family, "amplitude": amp,
                                "precision_dps": dps, "status": "baseline_error"})
            else:
                futures[pool.submit(evaluate_case, (n, family, amp, dps, base))] = (n, family, amp, dps)
        for fut in as_completed(futures):
            try:
                results.append(fut.result())
            except Exception as exc:
                n, family, amp, dps = futures[fut]
                results.append({"n": n, "family": family, "amplitude": amp,
                                "precision_dps": dps, "status": "worker_error",
                                "error": f"{type(exc).__name__}: {exc}"})
            if len(results) % 12 == 0:
                write_output(baseline_records, results, start, partial=True)

    results.sort(key=lambda x: (x.get("n", 0), x.get("family", ""),
                                AMPS.index(x.get("amplitude", "0.1")),
                                x.get("precision_dps", 0)))
    bykey = {(x.get("n"), x.get("family"), x.get("amplitude"), x.get("precision_dps")): x
             for x in results}
    diffs = []
    for n in NS:
        for family in ("sparse", "dense", "dense_energy"):
            for amp in AMPS:
                low = bykey.get((n, family, amp, 40))
                high = bykey.get((n, family, amp, 70))
                if low and high and low.get("status") == high.get("status") == "ok":
                    for field in ("C_all", "C_fixed", "C_singleaxis", "D_over_f", "Schur"):
                        try:
                            delta = abs(mp.mpf(low[field]) - mp.mpf(high[field]))
                            diffs.append({"n": n, "family": family, "amplitude": amp,
                                          "field": field, "abs_difference_40_vs_70": mp.nstr(delta, 12)})
                        except Exception:
                            pass
    write_output(baseline_records, results, start, partial=False, differences=diffs)


def write_output(baselines, results, start, partial, differences=None):
    output = {"status": "partial" if partial else
              ("complete" if len(results) == 72 * len(PRECISIONS) else "incomplete"),
              "purpose": "Finite-angle all-axis curvature diagnostics; numerical exploration only, not a rigorous certificate.",
              "dimensions": list(NS), "amplitudes": list(AMPS), "families": ["sparse", "dense", "dense_energy"],
              "case_count_expected_per_precision": 72, "precision_dps": list(PRECISIONS),
              "processes": WORKERS, "blas_threads_per_worker": 1,
              "hypergeometric_maxterms": MAXTERMS,
              "case_budget_seconds": CASE_BUDGET_SECONDS,
              "time_limit_semantics": "Each task reports elapsed time and flags elapsed>budget after completion; hypergeometric calls are additionally capped by maxterms. Exceptions and baseline failures are retained in JSON.",
              "baseline_equations": ["CDF_n(r,c)=1/2", "CDF_{n-1}(r-c,c)=1/2"],
              "cdf_model": "Gamma(N,1)+c*Binomial(N,p), p=1/(1+exp(c)); finite binomial-gamma mixture.",
              "numerical_scope": "mpmath hypergeometric and finite-mixture calculations; 40/70 dps agreement is a diagnostic, not interval certification.",
              "execution_correction": "Post-run audit found an output-only Hessian field error in the initial execution: H eigenvalues had been formed with D_over_f in place of h. The Schur and curvature formulas were unaffected. Corrected script now uses H_transverse=h-offb and H_longitudinal=h+(n-1)*offb; stored h/offb were used to recompute all 144 existing JSON H fields at their original dps without rerunning hypergeometric calculations.",
              "baseline_solutions": baselines, "results": results,
              "precision_differences": differences or [],
              "status_counts": {s: sum(1 for x in results if x.get("status") == s)
                                for s in ("ok", "timeout_budget_exceeded", "error", "worker_error", "baseline_error")},
              "wall_seconds": time.time() - start}
    (HERE / "curvature.json").write_text(json.dumps(output, indent=2) + "\n")
    print(json.dumps({"status": output["status"], "results": len(results),
                      "status_counts": output["status_counts"],
                      "wall_seconds": output["wall_seconds"]}, indent=2))


if __name__ == "__main__":
    main()
