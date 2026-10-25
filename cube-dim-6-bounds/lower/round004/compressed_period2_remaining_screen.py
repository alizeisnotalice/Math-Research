#!/usr/bin/env python3
"""Eight remaining reflected period-two witness screens (numeric only)."""
from functools import lru_cache
from fractions import Fraction
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
MAIN_LEDGER = Path("/Users/zhengzhihao/Desktop/T/Math-Research-cube-dim-6/cube-dim-6-bounds/results/bounds_ledger.json")


def base_spacing():
    u, v = 2 ** (1 / 6), 3 ** (1 / 6)
    return (1 + u**4 / 2, (u**3 + u**4 * v**3) / 4,
            *((1 + u) / 2 for _ in range(4)))


def reordered_spacing(target_zero_based):
    a = base_spacing()
    order = (target_zero_based,) + tuple(j for j in range(6) if j != target_zero_based)
    return tuple(a[j] for j in order), order


def centered_interval(a, ell, radius):
    start = -(ell // 2) if ell % 2 else 1 - ell // 2
    if ell % 2:
        return 0.0, min(a / 2, radius - (ell - 1) * a / 2), start
    return max(0.0, ell * a / 2 - radius), a / 2, start


def first_axis_starts(ell, rest_mass, delta, a1):
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


def boxes(target, delta, started):
    a, order = reordered_spacing(target)
    raw = []
    for ell in product(range(1, 5), repeat=6):
        if time.monotonic() - started > WALL_SECONDS:
            raise TimeoutError("20-minute wall-time budget reached")
        rest_mass = prod(ell[1:])
        for s, selected, radius, lo1, hi1 in first_axis_starts(ell[0], rest_mass, delta, a[0]):
            if hi1 <= lo1 + TOL:
                continue
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
                            raise AssertionError((target, ell, s, delta, j, x, atom, radius))
            if not valid:
                continue
            for x in (lo1, hi1):
                for atom in (s * a[0], (s + ell[0] - 1) * a[0]):
                    if abs(x - atom) > radius + 4e-12:
                        raise AssertionError((target, ell, s, delta, 0, x, atom, radius))
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
    return raw, kept, order


def box_union(B):
    @lru_cache(None)
    def rec(X):
        if not X:
            return 0.0
        if len(X[0]) == 2:
            ints = sorted(X)
            lo, hi = ints[0]
            total = 0.0
            for c, d in ints[1:]:
                if c > hi + TOL:
                    total += hi - lo
                    lo, hi = c, d
                elif d > hi:
                    hi = d
            return total + hi - lo
        total = 0.0
        cuts = sorted({x for b in X for x in b[:2]})
        for lo, hi in zip(cuts, cuts[1:]):
            if hi > lo + TOL:
                child = tuple(sorted({b[2:] for b in X if b[0] <= lo + TOL and hi <= b[1] + TOL}))
                total += (hi - lo) * rec(child)
        return total
    return rec(tuple(sorted(B)))


def one_case(target, delta):
    started = time.monotonic()
    raw, kept, order = boxes(target, delta, started)
    half = box_union(tuple(b for *_, b in kept))
    return {"period_axis_one_based": target + 1, "delta": str(delta),
            "coordinate_order_zero_based": list(order), "raw_boxes": len(raw),
            "kept_boxes": len(kept), "half_cell_union_float": half,
            "ratio_float": 32 * half,
            "ratio_formula": "32*|U_H|", "elapsed_seconds": time.monotonic() - started}


def main():
    ledger_bytes = MAIN_LEDGER.read_bytes()
    upper = Fraction(json.loads(ledger_bytes)["lower"]["exact_interval"][1])
    # Axis 1 / 1/20 was already evaluated in Sub-batch B.
    cases = [(0, Fraction(1, 10)), (0, Fraction(1, 5)),
             (1, Fraction(1, 20)), (1, Fraction(1, 10)), (1, Fraction(1, 5)),
             (2, Fraction(1, 20)), (2, Fraction(1, 10)), (2, Fraction(1, 5))]
    out = {"status": "running", "scope": "eight numeric witness screens only; no C6 claim",
           "normalization": "reflected period cell has mass 2 and union 64|U_H|, hence ratio 32|U_H|",
           "frozen_complete_upper": str(upper), "ledger_sha256": hashlib.sha256(ledger_bytes).hexdigest(),
           "raw_cap_each_case": RAW_CAP, "wall_seconds_each_case": WALL_SECONDS,
           "case_definition": "period axes 1,2,representative tail 3; axis1 delta=1/20 excluded as already run",
           "results": []}
    try:
        for target, d in cases:
            result = one_case(target, float(d))
            # float comparison is only a promotion trigger; exact work is required after it.
            result["margin_over_frozen_upper_float"] = result["ratio_float"] - float(upper)
            result["requires_exact_certification"] = result["margin_over_frozen_upper_float"] > 1e-7
            out["results"].append(result)
        out["status"] = "completed_numeric_screen"
    except (RuntimeError, TimeoutError) as exc:
        out["status"] = "budget_stop"
        out["reason"] = str(exc)
    out["best"] = max(out["results"], key=lambda x: x["ratio_float"], default=None)
    out["any_requires_exact_certification"] = any(x["requires_exact_certification"] for x in out["results"])
    out["code_sha256"] = hashlib.sha256(HERE.read_bytes()).hexdigest()
    HERE.with_name("compressed_period2_remaining_results.json").write_text(json.dumps(out, indent=2) + "\n")
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
