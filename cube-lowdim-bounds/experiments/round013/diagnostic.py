#!/usr/bin/env python3
"""Bounded mpmath diagnostic of the exterior lemma; no search or LP.

Prints actual run evidence as JSON. This is not an interval certificate.
Default mode preserves the first run's finite-precision u==a endpoint test;
this test can miss exact mathematical r=1 contact, so those values must not
be called a reliable boundary check. --boundary-only explicitly enforces
the mathematical r=1 endpoint and reproduces the bounded correction run.
At most 920 one-dimensional quadratures: 440 exterior-generator point
evaluations, each using at most two slices, plus 40 reference h integrals.
"""
import argparse
import json
import platform
import time

import mpmath as mp

mp.mp.dps = 60
quadratures = 0


def quad(f, intervals):
    global quadratures
    quadratures += 1
    if quadratures > 1000:
        raise RuntimeError("Diagnostic quadrature budget exceeded")
    return mp.quad(f, intervals)


def h(p, r):
    # Endpoint cancellation is exact when r=1.
    if r == 1:
        return 2 * quad(lambda y: (1+y)**p/(1-y)**2, [-1, 0]) + 2/(p-1)
    return 2 * quad(lambda y: (1-abs(y))**p/(r-y)**2, [-1, 0, 1])


def slice_integral(p, a, u, endpoint=False):
    if a <= 0:
        return mp.mpf(0)
    # Literal unscaled slice, independent of the h-based scaling identity.
    if endpoint or u == a:
        return 2 * quad(lambda x: (a+x)**p/(x-a)**2, [-a, 0]) + 2*a**(p-1)/(p-1)
    return 2 * quad(lambda x: (a-abs(x))**p/(x-u)**2, [-a, 0, a])


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--boundary-only", action="store_true",
                        help="Only recheck r=1 with an explicit exact endpoint flag")
    args = parser.parse_args()
    started = time.monotonic()
    powers = list(map(mp.mpf, ["2", "2.5", "3", "4"]))
    radii = list(map(mp.mpf, ["1", "1.00000001", "1.001", "1.1", "1.5", "1.999", "2", "2.001", "3", "5"]))
    if args.boundary_only:
        radii = [mp.mpf(1)]
    angles = list(map(mp.mpf, ["0", "0.0001", "0.01", "0.1", "0.25", "0.5", "0.75", "0.9", "0.99", "0.9999", "1"]))
    minimum = mp.inf
    record = None
    max_axis_error = mp.mpf(0)
    nonaxis_minimum = mp.inf
    point_count = 0
    for p in powers:
        for r in radii:
            reference = h(p, r)
            for q in angles:
                # q=u/r, so both axes, diagonal, and support contacts occur.
                u, v = q*r, (1-q)*r
                explicit_endpoint = args.boundary_only and r == 1
                actual = (slice_integral(p, 1-v, u, endpoint=explicit_endpoint)
                          + slice_integral(p, 1-u, v, endpoint=explicit_endpoint))
                gap = reference-actual
                point_count += 1
                if gap < minimum:
                    minimum = gap
                    record = {"p": str(p), "r": str(r), "q": str(q),
                              "h": mp.nstr(reference, 55), "Apsi": mp.nstr(actual, 55)}
                if q in (0, 1):
                    max_axis_error = max(max_axis_error, abs(gap))
                else:
                    nonaxis_minimum = min(nonaxis_minimum, gap)
    print(json.dumps({
        "status": "finite_mpmath_diagnostic_only",
        "python": platform.python_version(), "mpmath": mp.__version__,
        "decimal_precision": mp.mp.dps,
        "point_evaluations": point_count, "quadratures": quadratures,
        "budget": 1000, "minimum_h_minus_Apsi": mp.nstr(minimum, 55),
        "minimum_record": record,
        "minimum_nonaxis_gap": mp.nstr(nonaxis_minimum, 55),
        "maximum_axis_absolute_error": mp.nstr(max_axis_error, 55),
        "elapsed_seconds": time.monotonic()-started,
        "scope": ("44 fixed r=1 endpoint points" if args.boundary_only
                  else "440 fixed exterior points") + "; no random sampling, no LP, no new upper bound",
        "explicit_r1_endpoint_handling": args.boundary_only
    }, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
