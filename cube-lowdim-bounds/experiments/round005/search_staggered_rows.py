"""Round 5 discovery search: a staggered-row periodic atomic motif.

This is deliberately different from the earlier 2x1/2x2 product-grid
searches.  The fundamental cell has four alternating columns and two rows;
the second row is shifted by ``delta`` in the horizontal direction.  Its
weights may also have a distinct multiplier.  For each rectangular finite
index block we use its coordinate bounding box to make one centred-square
witness.  Thus every reported box is a valid finite witness; we make no
claim that this list exhausts the superlevel set.

The search is numerical discovery only.  It evaluates exactly BUDGET motifs
after the baseline control, with a fixed NumPy seed.  It records the complete
best candidate and enough parameters to rerun the experiment.
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
SEED = 2026102505
BUDGET = 5_000
KMAX = 6
LMAX = 4
PHI = 1.68550999335552518466528


def split_periodic(lo: float, hi: float, period: float) -> list[tuple[float, float]]:
    """Intersect [lo, hi] with a fundamental period, after quotienting."""
    length = hi - lo
    if length <= 0.0:
        return []
    if length >= period:
        return [(0.0, period)]
    lo = lo % period
    if lo + length <= period:
        return [(lo, lo + length)]
    return [(lo, period), (0.0, lo + length - period)]


def union_area(boxes: list[tuple[float, float, float, float]]) -> float:
    """Exact-for-the-supplied-floats rectangle-union scan, independent of search."""
    if not boxes:
        return 0.0
    xs = sorted({z for b in boxes for z in b[:2]})
    area = 0.0
    for left, right in zip(xs, xs[1:]):
        if right <= left:
            continue
        mid = (left + right) / 2.0
        intervals = sorted((b[2], b[3]) for b in boxes if b[0] < mid < b[1])
        if not intervals:
            continue
        covered = 0.0
        lo, hi = intervals[0]
        for c, d in intervals[1:]:
            if c > hi:
                covered += hi - lo
                lo, hi = c, d
            else:
                hi = max(hi, d)
        area += (right - left) * (covered + hi - lo)
    return area


def make_points(p: np.ndarray) -> tuple[np.ndarray, np.ndarray, float, float, float]:
    """Return nearby periodic points (coordinate, mass), periods, and cell mass."""
    xg = p[:4]
    yg = p[4:6]
    delta_fraction, heavy, row_multiplier = p[6:]
    px, py = float(sum(xg)), float(sum(yg))
    xbase = np.r_[0.0, np.cumsum(xg[:-1])]
    ybase = np.r_[0.0, yg[0]]
    column_weight = np.array([1.0, heavy, 1.0, heavy])
    pts: list[tuple[float, float]] = []
    masses: list[float] = []
    # The blocks below can start within the fundamental 4x2 index cell and
    # extend at most KMAX/LMAX.  This neighbourhood contains all their atoms.
    for j in range(-1, LMAX + 2):
        shift_x = (j % 2) * delta_fraction * px
        row_scale = row_multiplier if j % 2 else 1.0
        y = (j // 2) * py + ybase[j % 2]
        for i in range(-1, KMAX + 3):
            x = (i // 4) * px + xbase[i % 4] + shift_x
            pts.append((x, y))
            masses.append(float(column_weight[i % 4] * row_scale))
    return np.asarray(pts), np.asarray(masses), px, py, float(2.0 * (2.0 + 2.0 * heavy) * (1.0 + row_multiplier) / 2.0)


def witness_boxes(p: np.ndarray) -> tuple[list[tuple[float, float, float, float]], float]:
    """Generate all finite rectangular-block witnesses used in this round."""
    if (not np.all(np.isfinite(p))) or np.min(p[:6]) <= 0 or p[7] <= 0 or p[8] <= 0:
        return [], 1.0
    points, masses, px, py, cell_mass = make_points(p)
    # Index into the explicit point array: j maps to (j+1) row block and i to
    # (i+1) column block.  The supplied margin is enough for all KMAX/LMAX.
    def ix(i: int, j: int) -> int:
        return (j + 1) * (KMAX + 4) + (i + 1)

    raw: list[tuple[float, float, float, float]] = []
    for start_j in range(2):
        for start_i in range(4):
            for ly in range(1, LMAX + 1):
                for lx in range(1, KMAX + 1):
                    ids = [ix(i, j) for j in range(start_j, start_j + ly)
                           for i in range(start_i, start_i + lx)]
                    selected = points[ids]
                    side = math.sqrt(float(np.sum(masses[ids])))
                    xmin, ymin = np.min(selected, axis=0)
                    xmax, ymax = np.max(selected, axis=0)
                    lo_x, hi_x = xmax - side / 2.0, xmin + side / 2.0
                    lo_y, hi_y = ymax - side / 2.0, ymin + side / 2.0
                    if hi_x <= lo_x or hi_y <= lo_y:
                        continue
                    for xx, yy in itertools.product(split_periodic(lo_x, hi_x, px),
                                                     split_periodic(lo_y, hi_y, py)):
                        raw.append((xx[0], xx[1], yy[0], yy[1]))
    # Duplicates and contained rectangles cannot alter the union.
    unique = list(dict.fromkeys(raw))
    keep: list[tuple[float, float, float, float]] = []
    for i, b in enumerate(unique):
        contained = False
        for j, a in enumerate(unique):
            if i == j:
                continue
            if a[0] <= b[0] and a[1] >= b[1] and a[2] <= b[2] and a[3] >= b[3]:
                contained = True
                break
        if not contained:
            keep.append(b)
    return keep, cell_mass


def score(p: np.ndarray) -> tuple[float, int]:
    boxes, mass = witness_boxes(p)
    return union_area(boxes) / mass, len(boxes)


def decode(z: np.ndarray) -> np.ndarray:
    """Map a bounded latent vector to positive spacings and motif parameters."""
    xg = np.exp(z[:4])
    yg = np.exp(z[4:6])
    # z[6] is already a shift fraction.  Heavy and row multiplier are log scale.
    return np.r_[xg, yg, z[6], np.exp(z[7]), np.exp(z[8])]


def baseline() -> np.ndarray:
    u = (2.0 + math.sqrt(22.0)) / 3.0
    h = (1.0 + u) / 2.0
    v = h + 1.0
    w = u * u - 1.0
    return np.array([h, h, h, h, v, v, 0.0, w, 1.0])


def run() -> dict:
    start = time.monotonic()
    rng = np.random.default_rng(SEED)
    base = baseline()
    base_score, base_boxes = score(base)
    calls = 0
    best = {"score": base_score, "parameters": base.tolist(), "boxes": base_boxes, "source": "baseline_control"}
    all_top: list[tuple[float, np.ndarray, int, str]] = []

    # Budget allocation: 3000 broad draws + 2000 local perturbations around
    # the most promising staggered motifs, so this is a finite exploration.
    lower = np.array([math.log(0.35)] * 4 + [math.log(0.45)] * 2 + [0.0, math.log(0.45), math.log(0.35)])
    upper = np.array([math.log(3.8)] * 4 + [math.log(4.8)] * 2 + [0.999, math.log(8.5), math.log(3.2)])

    def evaluate(z: np.ndarray, source: str) -> None:
        nonlocal calls, best
        p = decode(np.clip(z, lower, upper))
        value, nboxes = score(p)
        calls += 1
        all_top.append((value, p, nboxes, source))
        if value > best["score"]:
            best = {"score": value, "parameters": p.tolist(), "boxes": nboxes, "source": source}
            print(f"{calls:5d} best={value:.15f} boxes={nboxes} via {source}", flush=True)

    for _ in range(3_000):
        evaluate(rng.uniform(lower, upper), "broad_uniform")
    all_top.sort(key=lambda item: item[0], reverse=True)
    seeds = all_top[:40]
    for k in range(2_000):
        _, p0, _, _ = seeds[k % len(seeds)]
        z0 = np.r_[np.log(p0[:6]), p0[6], np.log(p0[7:])]
        # Temperature cycles modestly to sample facets of the union objective.
        sigma = 0.11 if k < 1_200 else 0.045
        evaluate(z0 + rng.normal(0.0, sigma, size=9), "local_from_top40")
    assert calls == BUDGET
    top = sorted(all_top, key=lambda item: item[0], reverse=True)[:20]
    outcome = {
        "round": 5,
        "family": "4x2_staggered_rows_rectangular_block_witnesses",
        "status": "numerical_discovery_only",
        "seed": SEED,
        "budget_evaluations": BUDGET,
        "elapsed_seconds": time.monotonic() - start,
        "witness_block_bounds": {"horizontal_columns": [1, KMAX], "vertical_rows": [1, LMAX]},
        "baseline_control": {"score": base_score, "phi": PHI, "difference": base_score - PHI, "parameters": base.tolist(), "boxes": base_boxes},
        "best": best,
        "best_minus_phi": best["score"] - PHI,
        "top20": [{"score": v, "parameters": p.tolist(), "boxes": b, "source": s} for v, p, b, s in top],
        "environment": {"python": sys.version, "platform": platform.platform(), "numpy": np.__version__},
    }
    (ROOT / "search_result.json").write_text(json.dumps(outcome, indent=2, ensure_ascii=False) + "\n")
    print(json.dumps(outcome, indent=2, ensure_ascii=False), flush=True)
    return outcome


if __name__ == "__main__":
    run()
