#!/usr/bin/env python3
"""Full-rank original-source formulas and SMALL exact symbol guards.

No spatial MC runs here. Importing does not run anything or write files.
Use --guard ONLY after reviewing the saved registration.
"""
from fractions import Fraction as F
from pathlib import Path
import argparse
import hashlib
import json
import time

BASE = Path(__file__).resolve().parent
PREFIX = "full_rank_source_probe_20261007"


def atanh_log_interval(y, terms):
    """Exact log(y) bracket for 1 <= y <= 2, positive power-series tail."""
    y = F(y)
    assert 1 <= y <= 2
    z = (y - 1) / (y + 1)
    power = z
    value = F(0)
    for k in range(terms):
        value += 2 * power / (2 * k + 1)
        power *= z * z
    tail = 2 * power / ((2 * terms + 1) * (1 - z * z))
    return value, value + tail


def log_interval(y, terms):
    """Exact rational range reduction; no floating-point log/rounding."""
    y = F(y)
    assert y >= 1
    exponent = 0
    while y >= 2:
        y /= 2
        exponent += 1
    lo, hi = atanh_log_interval(y, terms)
    lo2, hi2 = atanh_log_interval(F(2), terms)
    return lo + exponent * lo2, hi + exponent * hi2


def b_interval(v, terms):
    v = F(v)
    if v == 0:
        return F(0), F(0)
    assert v > 0
    lo, hi = log_interval(1 + v, terms)
    return v / hi - 1, v / lo - 1


def p_squared_interval(r, scale, xi_squared, terms):
    """Original c=1 symbol P_r^scale squared, one marginal."""
    blo, bhi = b_interval(F(scale) ** 2 * F(xi_squared), terms)
    assert blo >= 0 and 0 < r < 1
    # 1-r B/(1+B) is DECREASING in B.
    return 1 - r * bhi / (1 + bhi), 1 - r * blo / (1 + blo)


def quotient_squared_interval(r, small, large, xi_squared, terms):
    alo, ahi = p_squared_interval(r, small, xi_squared, terms)
    blo, bhi = p_squared_interval(r, large, xi_squared, terms)
    assert alo > 0 and blo > 0
    return blo / ahi, bhi / alo


def rational_string(x):
    return f"{x.numerator}/{x.denominator}"


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def tensor_response(a_outer, a_layers, d_layers, amplitude, weights, kappa=1):
    """Original finite-r K_r nu tensor formula, arbitrary receiver batch.

    Caller provides 1D T_r tables already evaluated at all n coordinates.
    Last axis is coordinate. Prefix/suffix avoid division at zero.
    This FLOAT helper is not a certified continuum convolution or MC runner.
    It never clips negative values or adjusts thresholds.
    """
    import numpy as np
    a_outer = np.asarray(a_outer)
    n = a_outer.shape[-1]
    out = kappa * np.prod(a_outer, axis=-1)
    for weight, avec, dvec in zip(weights, a_layers, d_layers):
        avec, dvec = np.asarray(avec), np.asarray(dvec)
        assert avec.shape == dvec.shape and avec.shape[-1] == n
        ones = np.ones_like(avec[..., :1])
        prefix = np.concatenate((ones, np.cumprod(avec[..., :-1], axis=-1)), axis=-1)
        suffix = np.concatenate((np.cumprod(avec[..., :0:-1], axis=-1)[..., ::-1], ones), axis=-1)
        out = out + amplitude * weight * (n * np.prod(avec, axis=-1) - np.sum(dvec * prefix * suffix, axis=-1))
    return out


def count_statistics(counts, node_count):
    """Accepted-count identities; only algebra, no measure reinterpretation."""
    import numpy as np
    counts = np.asarray(counts)
    assert np.all((counts >= 0) & (counts <= node_count) & (counts == np.floor(counts)))
    mean = float(np.mean(counts))
    stop = float(np.mean(counts > 0))
    deficit = float(np.mean(np.maximum(counts - 1, 0)))
    pair = float(np.mean(counts * (counts - 1) / 2))
    return dict(weak_over_W=mean, Tstop_over_W=stop, D_over_W=deficit,
                pair_over_W=pair, identity_residual=mean-stop-deficit)


def run_guard():
    registration_path = BASE / f"{PREFIX}_registration.json"
    registration = json.loads(registration_path.read_text())
    result_path = BASE / f"{PREFIX}_guard_results.json"
    assert not result_path.exists(), "Never overwrite an old result."
    start = time.perf_counter()
    records = []
    corner_checks = []
    for plan in registration["executed_small_guard_rounds"]:
        n, terms = plan["n"], plan["log_series_terms"]
        for rstr in plan["r_values"]:
            r = F(rstr)
            for xstr in plan["xi_squared_values"]:
                x = F(xstr)
                lo, hi = quotient_squared_interval(r, F(1), F(2), x, terms)
                assert 0 < lo <= hi < 1
                reverse_lo, reverse_hi = 1 / hi, 1 / lo
                assert reverse_lo > 1
                records.append(dict(n=n, r=rstr, scale_small="1", scale_large="2",
                    xi_squared=xstr, terms=terms,
                    squared_ratio_interval=[rational_string(lo), rational_string(hi)],
                    squared_ratio_display=[float(lo), float(hi)],
                    reverse_squared_ratio_interval=[rational_string(reverse_lo), rational_string(reverse_hi)],
                    passed=True))
        for dstr in registration["source_formula_parameters"]["delta_layers"]:
            delta = F(dstr)
            L = F(101, 100) * F(16, 35) * delta
            assert 0 < L <= 1
            # Symmetry reduces all 2^n multiaffine vertices to zero-count classes.
            for zeros in range(n + 1):
                value = n * L - n if zeros == 0 else L if zeros == 1 else F(0)
                assert value <= L
                corner_checks.append(dict(n=n, delta=dstr, L=rational_string(L),
                    zeros=zeros, vertex_value=rational_string(value), passed=True))
    # Independent exact sanity checks on log series formula and source constants.
    assert atanh_log_interval(1, 8) == (F(0), F(0))
    assert F(8, 3)**6 > 201  # with pi>3 and e>8/3 gives coth(pi)<101/100.
    payload = dict(status="passed", scope="original c=1 marginal symbol and full-rank source corner algebra only",
        no_spatial_MC=True, no_weak_or_geom_estimate=True,
        registration_sha256=sha(registration_path), script_sha256=sha(Path(__file__)),
        symbol_checks=len(records), corner_checks=len(corner_checks), records=records,
        source_corner_records=corner_checks, elapsed_seconds=time.perf_counter()-start,
        rounding="Fraction throughout; displays are non-authoritative floating-point summaries",
        global_limit="Proved analytically in companion note; finite high-frequency nodes do not certify an infinite-frequency limit.",
        spatial_MC_status=registration["spatial_MC_status"])
    result_path.write_text(json.dumps(payload, indent=2) + "\n")
    print(json.dumps({k: payload[k] for k in ["status", "symbol_checks", "corner_checks", "elapsed_seconds"]}))
    print(str(result_path))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--guard", action="store_true")
    args = parser.parse_args()
    if args.guard:
        run_guard()
    else:
        parser.print_help()
