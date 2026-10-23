"""Independent edge-cell union check for the saved round-5 candidates.

This deliberately does not call ``union_area``.  It partitions a period cell
at every witness edge and classifies each open elementary cell by one interior
point.  Boundaries have zero area, so this gives the same Lebesgue union area
for the finite witness list if the discovery scan is correct.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

import search_staggered_rows as search

ROOT = Path(__file__).resolve().parent


def edge_cell_area(boxes: list[tuple[float, float, float, float]]) -> float:
    xs = sorted({x for b in boxes for x in b[:2]})
    ys = sorted({y for b in boxes for y in b[2:]})
    total = 0.0
    for left, right in zip(xs, xs[1:]):
        for bottom, top in zip(ys, ys[1:]):
            probe_x, probe_y = (left + right) / 2.0, (bottom + top) / 2.0
            if any(a < probe_x < b and c < probe_y < d for a, b, c, d in boxes):
                total += (right - left) * (top - bottom)
    return total


def check(label: str, record: dict) -> dict:
    boxes, mass = search.witness_boxes(np.array(record["parameters"], dtype=float))
    scan_area = search.union_area(boxes)
    grid_area = edge_cell_area(boxes)
    return {
        "label": label,
        "box_count": len(boxes),
        "mass_per_period": mass,
        "scan_area": scan_area,
        "edge_cell_area": grid_area,
        "absolute_area_difference": abs(scan_area - grid_area),
        "score_from_edge_cells": grid_area / mass,
        "recorded_score": record["score"],
        "absolute_score_difference": abs(grid_area / mass - record["score"]),
    }


def main() -> None:
    saved = json.loads((ROOT / "search_result.json").read_text())
    # ``best`` is the baseline when no new parameter wins.  Audit the strongest
    # genuinely non-control entry rather than checking that control twice.
    best_noncontrol = max(saved["top20"], key=lambda row: row["score"])
    out = {
        "method": "finite-edge-cell decomposition independent of the search scan",
        "checks": [check("baseline_control", saved["baseline_control"]),
                   check("best_noncontrol", best_noncontrol)],
    }
    for row in out["checks"]:
        assert row["absolute_area_difference"] < 1e-10
        assert row["absolute_score_difference"] < 1e-11
    (ROOT / "verification.json").write_text(json.dumps(out, indent=2, ensure_ascii=False) + "\n")
    print(json.dumps(out, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
