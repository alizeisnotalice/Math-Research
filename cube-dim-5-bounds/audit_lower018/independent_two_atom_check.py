#!/usr/bin/env python3
"""Registered exact enclosure audit; no imports from the candidate's programs."""
from fractions import Fraction as F
from itertools import combinations, product
from pathlib import Path
import hashlib
import json
import math
import platform
import sys
import time

HERE = Path(__file__).resolve().parent


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def qstr(x):
    return str(x.numerator) + "/" + str(x.denominator)


def fifth_root_enclosure(m, denominator):
    """Use integer comparisons only; detect exact roots at the chosen scale."""
    assert m >= 0
    lo, hi = 0, denominator
    target = m.numerator * denominator**5
    while hi**5 * m.denominator < target:
        hi *= 2
    while lo + 1 < hi:
        mid = (lo + hi) // 2
        if mid**5 * m.denominator <= target:
            lo = mid
        else:
            hi = mid
    if lo**5 * m.denominator == target:
        return F(lo, denominator), F(lo, denominator)
    if hi**5 * m.denominator == target:
        return F(hi, denominator), F(hi, denominator)
    assert lo**5 * m.denominator < target < hi**5 * m.denominator
    return F(lo, denominator), F(hi, denominator)


def coordinate_enclosure(spec, denominator):
    if isinstance(spec, str):
        x = F(spec)
        return x, x
    a, b = fifth_root_enclosure(F(spec["fifth_root"]), denominator)
    c = F(spec.get("factor", "1"))
    return (c*a, c*b) if c >= 0 else (c*b, c*a)


def receiver_boxes(case, denominator):
    weights = list(map(F, case["weights"]))
    alpha = F(case["alpha"])
    points = [[coordinate_enclosure(x, denominator) for x in p]
              for p in case["points"]]
    inner, outer, roots = [], [], []
    for k in (1, 2):
        for subset in combinations(range(2), k):
            mass = sum((weights[i] for i in subset), F(0))
            if mass == 0:
                continue
            slo, shi = fifth_root_enclosure(mass/alpha, denominator)
            roots.append({"subset": list(subset), "side_lower": qstr(slo),
                          "side_upper": qstr(shi), "mass_over_alpha": qstr(mass/alpha)})
            bi, bo = [], []
            for axis in range(5):
                lows = [points[i][axis][0] for i in subset]
                highs = [points[i][axis][1] for i in subset]
                bi.append((max(highs)-slo/2, min(lows)+slo/2))
                bo.append((max(lows)-shi/2, min(highs)+shi/2))
            if all(a < b for a, b in bi):
                inner.append(tuple(bi))
            if all(a < b for a, b in bo):
                outer.append(tuple(bo))
    return inner, outer, roots


def cell_union_volume(boxes, cell_limit, deadline):
    """Disjoint coordinate cells, with exact midpoint membership decisions."""
    if not boxes:
        return F(0), 0
    cuts = [sorted({z for box in boxes for z in box[j]}) for j in range(5)]
    intervals = [[(a, b) for a, b in zip(cut, cut[1:]) if a < b]
                 for cut in cuts]
    total, cells = F(0), 0
    for cell in product(*intervals):
        cells += 1
        assert cells <= cell_limit, "registered cell budget exceeded"
        if cells % 64 == 0:
            assert time.monotonic() < deadline, "registered time budget exceeded"
        midpoint = [(a+b)/2 for a, b in cell]
        if any(all(a < x < b for (a, b), x in zip(box, midpoint)) for box in boxes):
            total += math.prod(b-a for a, b in cell)
    return total, cells


def inclusion_exclusion_volume(boxes):
    """Small secondary check, implemented without coordinate subdivision."""
    total = F(0)
    for k in range(1, len(boxes)+1):
        for chosen in combinations(boxes, k):
            volume = F(1)
            for axis in range(5):
                a = max(box[axis][0] for box in chosen)
                b = min(box[axis][1] for box in chosen)
                volume *= max(F(0), b-a)
            total += volume if k % 2 else -volume
    return total


def main():
    registration_path = HERE / "registration.json"
    registration = json.loads(registration_path.read_text())
    assert digest(__file__) == registration["script_sha256"]
    for source in registration["locked_sources"]:
        assert digest(source["path"]) == source["sha256"], "source changed after registration"
    baseline_path = registration["standing_C_path"]
    baseline = json.loads(Path(baseline_path).read_text())
    clo, chi = map(F, baseline["ratio_interval"])
    assert F(3, 2) < clo <= chi
    denominator = 2**registration["root_denominator_power"]
    started = time.monotonic()
    outputs = []
    for case in registration["cases"]:
        ts = time.monotonic()
        deadline = ts + registration["per_case_seconds"]
        inner, outer, roots = receiver_boxes(case, denominator)
        vlow, nlow = cell_union_volume(inner, registration["per_union_cell_limit"], deadline)
        vhigh, nhigh = cell_union_volume(outer, registration["per_union_cell_limit"], deadline)
        assert vlow == inclusion_exclusion_volume(inner)
        assert vhigh == inclusion_exclusion_volume(outer)
        scale = F(case["alpha"])/sum(map(F, case["weights"]))
        rlow, rhigh = vlow*scale, vhigh*scale
        assert 0 <= rlow <= rhigh <= F(3, 2)
        if case.get("exact_expected_ratio"):
            expected = F(case["exact_expected_ratio"])
            assert rlow <= expected <= rhigh
        outputs.append({"id": case["id"], "ratio_interval": [qstr(rlow), qstr(rhigh)],
                        "display_only": [float(rlow), float(rhigh)],
                        "inner_outer_boxes": [len(inner), len(outer)],
                        "cell_counts": [nlow, nhigh],
                        "root_certificates": roots,
                        "two_union_methods_agree": True,
                        "exact_upper_at_most_three_halves": True,
                        "seconds": time.monotonic()-ts})
    by_id = {r["id"]: r for r in outputs}
    assert (by_id["small_mixed_unequal"]["ratio_interval"] ==
            by_id["mass_level_rescaling"]["ratio_interval"])
    result = {"status": "completed_registered_exact_degeneracy_audit",
              "scope": "Two nonnegative atoms in R^5 only; general all-scale theorem is proved separately.",
              "python": sys.executable, "python_version": platform.python_version(),
              "arithmetic": "Python standard-library Fraction and integer fifth-power comparisons",
              "registration_sha256": digest(registration_path),
              "script_sha256": digest(__file__),
              "locked_sources": registration["locked_sources"],
              "standing_C_interval": baseline["ratio_interval"],
              "three_halves_strictly_below_full_C_lower_endpoint": True,
              "old_C_baseline_rerun": False,
              "case_count": len(outputs), "cases": outputs,
              "mass_level_rescaling_identical": True,
              "elapsed_seconds": time.monotonic()-started}
    path = HERE / "results.json"
    path.write_text(json.dumps(result, indent=2)+"\n")
    print(json.dumps({"status": result["status"], "case_count": len(outputs),
                      "max_ratio_upper_display_only": max(r["display_only"][1] for r in outputs),
                      "elapsed_seconds": result["elapsed_seconds"],
                      "results_sha256": digest(path)}))


if __name__ == "__main__":
    main()
