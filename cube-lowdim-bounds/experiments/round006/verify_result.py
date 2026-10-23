"""Separate endpoint-cell audit for the round-6 saved baseline and winner."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

import search_four_weights as search

ROOT = Path(__file__).resolve().parent


def check(dimension: int, label: str, record: dict) -> dict:
    boxes = search.boxes_for(np.asarray(record["spacings"], dtype=float),
                             np.asarray(record["weights"], dtype=float))
    full_mass = float(np.sum(np.asarray(record["weights"], dtype=float)))
    # Deliberately only use endpoint-cell marking here, not the recursive
    # sweep used by the search objective.
    volume = search.independent_grid_volume(boxes)
    return {"dimension": dimension, "label": label, "box_count": len(boxes),
            "full_period_mass": full_mass, "endpoint_cell_volume": volume,
            "endpoint_cell_score": volume / full_mass, "recorded_score": record["score"],
            "absolute_score_difference": abs(volume / full_mass - record["score"])}


def main() -> None:
    result = json.loads((ROOT / "search_result.json").read_text())
    rows = []
    for entry in result["dimensions"]:
        candidates = entry["structured_candidates"] + entry["top_joint_candidates"]
        noncontrol = max((row for row in candidates if row["evaluation"] != 1),
                         key=lambda row: row["score"])
        rows.append(check(entry["dimension"], "baseline", entry["baseline"]))
        rows.append(check(entry["dimension"], "best_noncontrol", noncontrol))
    for row in rows:
        assert row["absolute_score_difference"] < 1e-10
    output = {"method": "endpoint-cell marking separate from recursive sweep", "checks": rows}
    (ROOT / "verification.json").write_text(json.dumps(output, indent=2, ensure_ascii=False) + "\n")
    print(json.dumps(output, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
