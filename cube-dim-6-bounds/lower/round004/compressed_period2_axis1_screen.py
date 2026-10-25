#!/usr/bin/env python3
"""Reflection-compressed period-two, axis-one discovery screen.

This uses only a reflected witness subfamily in
H=[0,a1] x product_{j>1}[0,aj/2].  It is a floating discovery calculation,
not an exact certificate.  The generated boxes retain one common radius per
source block and are checked against every selected atom in that block.
"""
from functools import lru_cache
from itertools import product
from math import ceil, floor, prod
from pathlib import Path
import hashlib
import json
import time

TOL = 2e-12
RAW_CAP = 5000
WALL_SECONDS = 20 * 60
HERE = Path(__file__).resolve()
# This independent lower-bound clone intentionally contains only lower/.
# The frozen comparator is read (never written) from the integration clone.
MAIN_LEDGER = Path("/Users/zhengzhihao/Desktop/T/Math-Research-cube-dim-6/cube-dim-6-bounds/results/bounds_ledger.json")


def spacing():
    u, v = 2 ** (1 / 6), 3 ** (1 / 6)
    return (1 + u**4 / 2, (u**3 + u**4 * v**3) / 4,
            *((1 + u) / 2 for _ in range(4)))


def centered_interval(a, ell, radius):
    """Old R6 witness interval in [0,a/2], with its source endpoints."""
    start = -(ell // 2) if ell % 2 else 1 - ell // 2
    if ell % 2:
        lo, hi = 0.0, min(a / 2, radius - (ell - 1) * a / 2)
    else:
        lo, hi = max(0.0, ell * a / 2 - radius), a / 2
    return lo, hi, start


def first_axis_starts(ell, rest_mass, delta, a1):
    """All integer s whose block center interval can meet [0,a1].

    For a fixed s, I_s=[(s+ell-1)a1-r_s, s*a1+r_s].  We enumerate a
    proven finite superset using r_max, then impose I_s intersection with
    the actual parity-dependent radius.  This is the exact range condition,
    not a fixed translate grid.
    """
    largest = ell if ell % 2 == 0 else ell + delta
    rmax = (rest_mass * largest) ** (1 / 6) / 2
    left = floor(-rmax / a1) - 1
    right = ceil(1 + rmax / a1 - (ell - 1)) + 1
    for s in range(left, right + 1):
        selected = ell if ell % 2 == 0 else ell + (delta if s % 2 == 0 else -delta)
        radius = (rest_mass * selected) ** (1 / 6) / 2
        lo, hi = (s + ell - 1) * a1 - radius, s * a1 + radius
        if lo < a1 - TOL and hi > TOL:
            yield s, selected, radius, max(0.0, lo), min(a1, hi)


def boxes(delta, started):
    """Return raw and maximal boxes of the reflected witness subfamily."""
    a = spacing()
    raw = []
    for ell in product(range(1, 5), repeat=6):
        if time.monotonic() - started > WALL_SECONDS:
            raise TimeoutError("20-minute wall-time budget reached")
        rest_mass = prod(ell[1:])
        # The five unchanged directions use exactly the old centered witnesses.
        rest = []
        valid = True
        base_radius_for_check = None
        for s, selected, radius, lo1, hi1 in first_axis_starts(ell[0], rest_mass, delta, a[0]):
            coords = [lo1, hi1]
            valid = True
            for j in range(1, 6):
                lo, hi, start = centered_interval(a[j], ell[j], radius)
                if hi <= lo + TOL:
                    valid = False
                    break
                coords += [lo, hi]
                for x in (lo, hi):
                    for atom in (start * a[j], (start + ell[j] - 1) * a[j]):
                        if abs(x - atom) > radius + 4e-12:
                            raise AssertionError((ell, s, delta, j, x, atom, radius))
            if not valid or hi1 <= lo1 + TOL:
                continue
            for x in (lo1, hi1):
                for atom in (s * a[0], (s + ell[0] - 1) * a[0]):
                    if abs(x - atom) > radius + 4e-12:
                        raise AssertionError((ell, s, delta, 0, x, atom, radius))
            raw.append((ell, s, selected, tuple(coords)))
            if len(raw) > RAW_CAP:
                raise RuntimeError("raw-box cap exceeded")
    unique = {box: (ell, s, selected) for ell, s, selected, box in raw}
    kept = []
    for box, label in unique.items():
        if any(other != box and all(other[2*j] <= box[2*j] + TOL and
                                    box[2*j+1] <= other[2*j+1] + TOL
                                    for j in range(6)) for other in unique):
            continue
        kept.append((*label, box))
    return raw, kept


def box_union(B):
    @lru_cache(None)
    def rec(X):
        if not X:
            return 0.0
        if len(X[0]) == 2:
            intervals = sorted(X)
            lo, hi = intervals[0]
            answer = 0.0
            for c, d in intervals[1:]:
                if c > hi + TOL:
                    answer += hi - lo
                    lo, hi = c, d
                elif d > hi:
                    hi = d
            return answer + hi - lo
        answer = 0.0
        cuts = sorted({x for b in X for x in b[:2]})
        for lo, hi in zip(cuts, cuts[1:]):
            if hi <= lo + TOL:
                continue
            child = tuple(sorted({b[2:] for b in X if b[0] <= lo + TOL and hi <= b[1] + TOL}))
            answer += (hi - lo) * rec(child)
        return answer
    return rec(tuple(sorted(B)))


def evaluate(delta, started):
    raw, kept = boxes(delta, started)
    half_union = box_union(tuple(b for *_, b in kept))
    return {"delta": f"{delta:.17g}", "raw_boxes": len(raw), "kept_boxes": len(kept),
            "half_cell_union_float": half_union,
            "ratio_float": 32 * half_union,
            "ratio_formula": "32 * |U_H| (full periodic reflected mass is 2)"}


def main():
    started = time.monotonic()
    ledger_bytes = MAIN_LEDGER.read_bytes()
    out = {"status": "running", "scope": "numeric reflected-witness discovery only; no C6 claim",
           "family": "axis-1 period-two weights 1+delta,1-delta at breakpoint spacing",
           "domain": "H=[0,a1] x product_{j>1}[0,aj/2]",
           "normalization": "full reflected union=64*|U_H|; period source mass=2; ratio=32*|U_H|",
           "raw_cap": RAW_CAP, "wall_seconds": WALL_SECONDS,
           "ledger_sha256": hashlib.sha256(ledger_bytes).hexdigest(),
           "frozen_complete_upper": json.loads(ledger_bytes)["lower"]["exact_interval"][1]}
    try:
        cal = evaluate(0.0, started)
        # The actual certified breakpoint value, used only as a float screen.
        baseline = 2.06692382601088
        cal["baseline_float"] = baseline
        cal["absolute_error_from_existing_witness_union"] = abs(cal["ratio_float"] - baseline)
        out["delta_zero_calibration"] = cal
        if cal["absolute_error_from_existing_witness_union"] > 2e-9:
            out["status"] = "calibration_mismatch_stop"
            out["delta_1_20_started"] = False
        else:
            out["delta_1_20"] = evaluate(1 / 20, started)
            out["delta_1_20_started"] = True
            out["status"] = "completed_numeric_screen"
    except (RuntimeError, TimeoutError) as exc:
        out["status"] = "budget_stop"
        out["reason"] = str(exc)
        out.setdefault("delta_1_20_started", False)
    out["elapsed_seconds"] = time.monotonic() - started
    out["code_sha256"] = hashlib.sha256(HERE.read_bytes()).hexdigest()
    HERE.with_name("compressed_period2_axis1_results.json").write_text(json.dumps(out, indent=2) + "\n")
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
