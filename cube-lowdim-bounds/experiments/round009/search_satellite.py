"""Round 9: one freely placed positive satellite per Lin 2h x V period cell.

For every candidate this script constructs the 3x3 period patch, enumerates
every coordinate-rectangular closed atom subset by its x/y endpoint pair,
and uses all atoms of that rectangle for a centred-square witness.  A witness
with mass M has side sqrt(M), hence uses one common radius in both coordinates.
Only these finite witnesses are counted, so omitted patch or larger-scale
witnesses can lower a score but never create a false lower bound.

Numerical discovery only; a value above Phi would be a candidate requiring a
separate rational/interval certificate.
"""

from __future__ import annotations

import functools
import itertools
import json
import math
import platform
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent
SEED = 2026102909
BUDGET = 500
PHI = 1.68550999335552518466528


def lin_parameters() -> tuple[float, float, float, float, float]:
    u = (2.0 + math.sqrt(22.0)) / 3.0
    h = (1.0 + u) / 2.0
    v = h + 1.0
    w = u * u - 1.0
    return h, v, w, 2.0 * h, v


def split_periodic(lo: float, hi: float, period: float) -> list[tuple[float, float]]:
    length = hi - lo
    if length <= 0.0:
        return []
    if length >= period:
        return [(0.0, period)]
    lo %= period
    if lo + length <= period:
        return [(lo, lo + length)]
    return [(lo, period), (0.0, lo + length - period)]


def prune(boxes: list[tuple[float, float, float, float]]) -> list[tuple[float, float, float, float]]:
    unique = list(dict.fromkeys(boxes))
    if not unique:
        return []
    a = np.asarray(unique)
    keep: list[tuple[float, float, float, float]] = []
    for i, box in enumerate(a):
        contains = np.all(a[:, ::2] <= box[::2], axis=1) & np.all(a[:, 1::2] >= box[1::2], axis=1)
        contains[i] = False
        if not bool(np.any(contains)):
            keep.append(tuple(float(z) for z in box))
    return keep


def union_area(boxes: list[tuple[float, float, float, float]]) -> float:
    @functools.lru_cache(maxsize=None)
    def rec(data: tuple[tuple[float, ...], ...]) -> float:
        if not data:
            return 0.0
        if len(data[0]) == 2:
            intervals = sorted(data)
            lo, hi = intervals[0]
            total = 0.0
            for a, b in intervals[1:]:
                if a > hi:
                    total += hi - lo
                    lo, hi = a, b
                else:
                    hi = max(hi, b)
            return total + hi - lo
        cuts = sorted({z for box in data for z in box[:2]})
        total = 0.0
        for lo, hi in zip(cuts, cuts[1:]):
            mid = (lo + hi) / 2.0
            rest = tuple(sorted(set(box[2:] for box in data if box[0] < mid < box[1])))
            total += (hi - lo) * rec(rest)
        return total
    return rec(tuple(sorted(boxes)))


def patch_points(x_fraction: float, y_fraction: float, satellite_mass: float) -> tuple[np.ndarray, np.ndarray, float, float, float]:
    """Build the requested 3x3 period patch and return complete cell mass."""
    h, v, w, px, py = lin_parameters()
    seed_points = [(0.0, 0.0, 1.0), (h, 0.0, w)]
    if satellite_mass > 0.0:
        seed_points.append((x_fraction * px, y_fraction * py, satellite_mass))
    points: list[tuple[float, float]] = []
    weights: list[float] = []
    for kx in (-1, 0, 1):
        for ky in (-1, 0, 1):
            for x, y, mass in seed_points:
                points.append((x + kx * px, y + ky * py))
                weights.append(mass)
    return np.asarray(points), np.asarray(weights), px, py, 1.0 + w + satellite_mass


def witness_boxes(x_fraction: float, y_fraction: float, satellite_mass: float) -> tuple[list[tuple[float, float, float, float]], float]:
    points, weights, px, py, period_mass = patch_points(x_fraction, y_fraction, satellite_mass)
    xs = sorted(set(float(x) for x in points[:, 0]))
    ys = sorted(set(float(y) for y in points[:, 1]))
    raw: list[tuple[float, float, float, float]] = []
    # Every x/y endpoint choice produces the full coordinate-rectangular
    # closure in this patch, not a selected subset with hidden omitted atoms.
    for left, right in itertools.combinations_with_replacement(xs, 2):
        xmask = (points[:, 0] >= left) & (points[:, 0] <= right)
        for bottom, top in itertools.combinations_with_replacement(ys, 2):
            mask = xmask & (points[:, 1] >= bottom) & (points[:, 1] <= top)
            if not bool(np.any(mask)):
                continue
            selected = points[mask]
            side = math.sqrt(float(weights[mask].sum()))
            lo_x, hi_x = float(selected[:, 0].max() - side / 2.0), float(selected[:, 0].min() + side / 2.0)
            lo_y, hi_y = float(selected[:, 1].max() - side / 2.0), float(selected[:, 1].min() + side / 2.0)
            if hi_x <= lo_x or hi_y <= lo_y:
                continue
            for xpart, ypart in itertools.product(split_periodic(lo_x, hi_x, px), split_periodic(lo_y, hi_y, py)):
                raw.append((xpart[0], xpart[1], ypart[0], ypart[1]))
    return prune(raw), period_mass


def score(x_fraction: float, y_fraction: float, satellite_mass: float) -> tuple[float, list[tuple[float, float, float, float]], float]:
    boxes, mass = witness_boxes(x_fraction, y_fraction, satellite_mass)
    return union_area(boxes) / mass, boxes, mass


def independent_edge_cells(boxes: list[tuple[float, float, float, float]]) -> float:
    """Independent endpoint-cell marking, not the recursive discovery sweep."""
    xs = sorted({z for box in boxes for z in box[:2]})
    ys = sorted({z for box in boxes for z in box[2:]})
    total = 0.0
    for left, right in zip(xs, xs[1:]):
        for bottom, top in zip(ys, ys[1:]):
            x, y = (left + right) / 2.0, (bottom + top) / 2.0
            if any(a < x < b and c < y < d for a, b, c, d in boxes):
                total += (right - left) * (top - bottom)
    return total


def schedule() -> list[tuple[str, float, float, float]]:
    rng = np.random.default_rng(SEED)
    out: list[tuple[str, float, float, float]] = [
        ("zero_mass_control", 0.0, 0.0, 0.0),
        ("near_zero_control", 0.371, 0.619, 1e-12),
    ]
    # Structured small masses probe each quarter and cell midlines before
    # broad random samples; all actual satellite candidates respect m>=1e-4.
    locations = [(x, y) for x in (0.0, 0.25, 0.5, 0.75) for y in (0.0, 0.25, 0.5, 0.75)]
    for m in (1e-4, 1e-3, 1e-2):
        for x, y in locations:
            out.append(("structured_small_mass", x, y, m))
    while len(out) < BUDGET:
        x, y = rng.uniform(0.0, 1.0, size=2)
        m = math.exp(rng.uniform(math.log(1e-4), math.log(1.0)))
        out.append(("joint_loguniform", float(x), float(y), float(m)))
    assert len(out) == BUDGET
    return out


def main() -> None:
    started = time.monotonic()
    records: list[dict] = []
    best_noncontrol = None
    for index, (source, x, y, mass) in enumerate(schedule(), 1):
        value, boxes, period_mass = score(x, y, mass)
        row = {"evaluation": index, "source": source, "x_fraction": x, "y_fraction": y,
               "satellite_mass": mass, "score": value, "box_count": len(boxes), "period_mass": period_mass}
        records.append(row)
        if mass >= 1e-4 and (best_noncontrol is None or value > best_noncontrol["score"]):
            best_noncontrol = row
            print(f"eval={index:3d} best_noncontrol={value:.15f} m={mass:.6g}", flush=True)
    zero, near_zero = records[:2]
    assert best_noncontrol is not None
    # Keep only audit-relevant data: controls, strongest satellite, and the
    # strongest structured and random satellite rows.
    structured_best = max((r for r in records if r["source"] == "structured_small_mass"), key=lambda r: r["score"])
    random_best = max((r for r in records if r["source"] == "joint_loguniform"), key=lambda r: r["score"])
    checks = []
    for label, row in (("zero_mass_control", zero), ("best_noncontrol", best_noncontrol)):
        _, boxes, period_mass = score(row["x_fraction"], row["y_fraction"], row["satellite_mass"])
        swept = union_area(boxes)
        edged = independent_edge_cells(boxes)
        checks.append({"label": label, "box_count": len(boxes), "period_mass": period_mass,
                       "sweep_area": swept, "edge_cell_area": edged,
                       "absolute_area_difference": abs(swept - edged),
                       "edge_cell_score": edged / period_mass,
                       "recorded_score": row["score"],
                       "absolute_score_difference": abs(edged / period_mass - row["score"])})
    result = {
        "round": 9, "family": "Lin_2h_by_V_two_atoms_plus_one_free_periodic_satellite",
        "status": "numerical_discovery_only", "seed": SEED, "evaluation_budget": BUDGET,
        "evaluations_completed": len(records), "patch_periods": [3, 3],
        "enumeration": "all x/y endpoint coordinate-rectangular closures within finite 3x3 patch",
        "controls": {"zero_mass": zero, "near_zero": near_zero, "phi": PHI,
                     "zero_mass_minus_phi": zero["score"] - PHI, "near_zero_minus_phi": near_zero["score"] - PHI},
        "best_noncontrol": best_noncontrol,
        "best_minus_phi": best_noncontrol["score"] - PHI,
        "strict_improvement_candidate": bool(best_noncontrol["score"] > PHI + 1e-10),
        "best_structured": structured_best, "best_random": random_best,
        "independent_edge_cell_checks": checks,
        "elapsed_seconds": time.monotonic() - started,
        "environment": {"python": sys.version, "platform": platform.platform(), "numpy": np.__version__},
    }
    (ROOT / "search_result.json").write_text(json.dumps(result, indent=2, ensure_ascii=False) + "\n")
    print(json.dumps(result, indent=2, ensure_ascii=False), flush=True)


if __name__ == "__main__":
    main()
