"""Round 6: four independent parity weights on a 2x2x1^(d-2) lattice.

The lattice point indexed by n in Z^d has coordinate n_j*a_j and weight
W[n_0 mod 2, n_1 mod 2].  The fundamental period is
(2*a_0, 2*a_1, a_2, ...), with mass sum(W), so no reflection multiplier is
used.  A rectangular atom block gives a valid same-radius cube witness: its
centre set is the coordinatewise intersection of radius sqrt[d](mass)/2
intervals.  Wrapped pieces are unioned inside the actual period cell.

This is bounded numerical discovery, not a certificate.  Every run evaluates
at most 200 parameter choices per dimension, with fixed seeds.
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
LIM = 4
MAX_EVALS = 200
SEEDS = {3: 2026102603, 4: 2026102604}
BASELINE = {3: 1.7939332277947704, 4: 1.9017426111009557}


def base_spacings(d: int) -> np.ndarray:
    u = 2.0 ** (1.0 / d)
    v = 3.0 ** (1.0 / d)
    return np.array([1.0 + u ** (d - 2) / 2.0, (2.0 + u ** (d - 1) * v) / 4.0]
                    + [(1.0 + u) / 2.0] * (d - 2), dtype=float)


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


def prune(boxes: list[tuple[float, ...]]) -> list[tuple[float, ...]]:
    unique = list(dict.fromkeys(boxes))
    if not unique:
        return []
    arr = np.asarray(unique)
    keep: list[tuple[float, ...]] = []
    for i, b in enumerate(arr):
        contains = np.all(arr[:, ::2] <= b[::2], axis=1) & np.all(arr[:, 1::2] >= b[1::2], axis=1)
        contains[i] = False
        if not bool(np.any(contains)):
            keep.append(tuple(float(x) for x in b))
    return keep


def union_volume(boxes: list[tuple[float, ...]]) -> float:
    """Recursive sweep used only by the discovery objective."""
    @functools.lru_cache(maxsize=None)
    def rec(data: tuple[tuple[float, ...], ...]) -> float:
        if not data:
            return 0.0
        if len(data[0]) == 2:
            intervals = sorted(data)
            left, right = intervals[0]
            total = 0.0
            for a, b in intervals[1:]:
                if a > right:
                    total += right - left
                    left, right = a, b
                else:
                    right = max(right, b)
            return total + right - left
        cuts = sorted({x for b in data for x in b[:2]})
        total = 0.0
        for left, right in zip(cuts, cuts[1:]):
            # With all endpoints as cuts, using the interior is equivalent and
            # avoids accidental boundary-only contributions.
            mid = (left + right) / 2.0
            rest = tuple(sorted(set(b[2:] for b in data if b[0] < mid < b[1])))
            total += (right - left) * rec(rest)
        return total
    return rec(tuple(sorted(boxes)))


def block_mass(weights: np.ndarray, starts: tuple[int, ...], lengths: tuple[int, ...]) -> float:
    first = [sum(weights[(starts[0] + i) % 2, (starts[1] + j) % 2]
                       for i in range(lengths[0]) for j in range(lengths[1]))]
    for ell in lengths[2:]:
        first.append(float(ell))
    return math.prod(first)


def boxes_for(spacings: np.ndarray, weights: np.ndarray) -> list[tuple[float, ...]]:
    d = len(spacings)
    periods = np.r_[2.0 * spacings[:2], spacings[2:]]
    raw: list[tuple[float, ...]] = []
    # Starts 0 and 1 are the two genuine residue classes in each weighted
    # direction.  Every other dimension has period one lattice interval.
    for starts01 in itertools.product(range(2), repeat=2):
        for lengths in itertools.product(range(1, LIM + 1), repeat=d):
            starts = starts01 + (0,) * (d - 2)
            mass = block_mass(weights, starts, lengths)
            side = mass ** (1.0 / d)
            pieces: list[list[tuple[float, float]]] = []
            for j, (s, ell, gap, period) in enumerate(zip(starts, lengths, spacings, periods)):
                low = (s + ell - 1) * gap - side / 2.0
                high = s * gap + side / 2.0
                if high <= low:
                    pieces = []
                    break
                pieces.append(split_periodic(low, high, period))
            if pieces:
                for product_piece in itertools.product(*pieces):
                    raw.append(tuple(x for pair in product_piece for x in pair))
    return prune(raw)


def score(spacings: np.ndarray, weights: np.ndarray) -> tuple[float, list[tuple[float, ...]]]:
    if np.min(spacings) <= 0.0 or np.min(weights) <= 0.0:
        return -math.inf, []
    boxes = boxes_for(spacings, weights)
    return union_volume(boxes) / float(np.sum(weights)), boxes


def independent_grid_volume(boxes: list[tuple[float, ...]]) -> float:
    """Endpoint-cell marking, independent from the recursive sweep above."""
    d = len(boxes[0]) // 2
    edges = [sorted({b[2 * j + q] for b in boxes for q in (0, 1)}) for j in range(d)]
    shape = tuple(len(edge) - 1 for edge in edges)
    covered = np.zeros(shape, dtype=bool)
    index = [{x: i for i, x in enumerate(edge)} for edge in edges]
    for box in boxes:
        slices = tuple(slice(index[j][box[2 * j]], index[j][box[2 * j + 1]]) for j in range(d))
        covered[slices] = True
    cell = np.ones(shape, dtype=float)
    for j, edge in enumerate(edges):
        widths = np.diff(np.asarray(edge, dtype=float))
        reshape = [1] * d
        reshape[j] = len(widths)
        cell *= widths.reshape(reshape)
    return float(cell[covered].sum())


def normalize_weights(log_weights: np.ndarray) -> np.ndarray:
    z = log_weights - float(np.mean(log_weights))
    return np.exp(z).reshape(2, 2)


def parameter_schedule(d: int) -> list[tuple[str, np.ndarray, np.ndarray]]:
    """One baseline, structured independent changes, then fixed random probes."""
    rng = np.random.default_rng(SEEDS[d])
    base = base_spacings(d)
    out: list[tuple[str, np.ndarray, np.ndarray]] = [("uniform_control", base, np.ones((2, 2)))]
    # Vary each of the four weights alone (renormalized by geometric mean),
    # including directions outside the old two-weight checkerboard family.
    for q in range(4):
        for t in (-0.45, -0.18, -0.06, 0.06, 0.18, 0.45):
            z = np.zeros(4)
            z[q] = t
            out.append((f"single_weight_{q}_{t:+.2f}", base, normalize_weights(z)))
    # Independent row/column/checkerboard modes are retained as controls.
    modes = [np.array([1.0, -1.0, 1.0, -1.0]), np.array([1.0, 1.0, -1.0, -1.0]),
             np.array([1.0, -1.0, -1.0, 1.0])]
    for m, mode in enumerate(modes):
        for t in (-0.30, -0.12, 0.12, 0.30):
            out.append((f"weight_mode_{m}_{t:+.2f}", base, normalize_weights(t * mode)))
    for j in range(d):
        for t in (-0.12, -0.05, 0.05, 0.12):
            space = base.copy()
            space[j] *= math.exp(t)
            out.append((f"single_spacing_{j}_{t:+.2f}", space, np.ones((2, 2))))
    # The remaining evaluations are true multi-parameter perturbations.  The
    # variance is deliberately bounded; this is a local finite study, not a
    # statement about the global four-weight family.
    while len(out) < MAX_EVALS:
        log_space = rng.normal(0.0, 0.075, size=d)
        log_weight = rng.normal(0.0, 0.23, size=4)
        out.append(("joint_gaussian", base * np.exp(log_space), normalize_weights(log_weight)))
    assert len(out) == MAX_EVALS
    return out


def run_dimension(d: int) -> dict:
    started = time.monotonic()
    candidates: list[dict] = []
    best = None
    for number, (source, spacings, weights) in enumerate(parameter_schedule(d), 1):
        value, boxes = score(spacings, weights)
        item = {"evaluation": number, "source": source, "score": value,
                "spacings": spacings.tolist(), "weights": weights.tolist(), "box_count": len(boxes)}
        candidates.append(item)
        if best is None or value > best["score"]:
            best = item
            print(f"d={d} eval={number:3d} best={value:.15f} source={source} boxes={len(boxes)}", flush=True)
    assert best is not None
    base = candidates[0]
    # Store a small, deterministic subset sufficient for review and rerunning;
    # it includes the best noncontrol structured and joint probes.
    structured = sorted((c for c in candidates if c["source"] != "joint_gaussian" and c["evaluation"] != 1), key=lambda c: c["score"], reverse=True)[:1]
    top_joint = sorted((c for c in candidates if c["source"] == "joint_gaussian"), key=lambda c: c["score"], reverse=True)[:1]
    outcome = {
        "dimension": d,
        "family": "2x2x1^(d-2)_periodic_lattice_four_independent_positive_weights",
        "seed": SEEDS[d], "evaluation_budget": MAX_EVALS, "evaluations_completed": len(candidates),
        "witness_lengths_each_coordinate": [1, LIM],
        "normalization": "geometric_mean(weights)=1; score uses full period union volume / sum(weights)",
        "baseline": base, "baseline_reference": BASELINE[d],
        "baseline_difference": base["score"] - BASELINE[d],
        "best": best, "best_minus_baseline": best["score"] - BASELINE[d],
        "strict_improvement_candidate": bool(best["score"] > BASELINE[d] + 1e-10),
        "structured_candidates": structured, "top_joint_candidates": top_joint,
        "elapsed_seconds": time.monotonic() - started,
    }
    return outcome


def main() -> None:
    result = {"round": 6, "status": "numerical_discovery_only", "dimensions": [run_dimension(3), run_dimension(4)],
              "environment": {"python": sys.version, "platform": platform.platform(), "numpy": np.__version__}}
    (ROOT / "search_result.json").write_text(json.dumps(result, indent=2, ensure_ascii=False) + "\n")
    print(json.dumps(result, indent=2, ensure_ascii=False), flush=True)


if __name__ == "__main__":
    main()
