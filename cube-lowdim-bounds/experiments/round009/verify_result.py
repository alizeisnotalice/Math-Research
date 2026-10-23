"""Independent endpoint-cell check for the saved round-9 control and winner."""

from __future__ import annotations

import json
from pathlib import Path

import search_satellite as search

ROOT = Path(__file__).resolve().parent


def check(label: str, row: dict) -> dict:
    boxes, period_mass = search.witness_boxes(row["x_fraction"], row["y_fraction"], row["satellite_mass"])
    # This deliberately uses endpoint-cell marking only.  It does not call the
    # recursive union sweep used by the discovery objective.
    area = search.independent_edge_cells(boxes)
    return {"label": label, "box_count": len(boxes), "period_mass": period_mass,
            "edge_cell_area": area, "edge_cell_score": area / period_mass,
            "recorded_score": row["score"],
            "absolute_score_difference": abs(area / period_mass - row["score"])}


def main() -> None:
    saved = json.loads((ROOT / "search_result.json").read_text())
    rows = [check("zero_mass_control", saved["controls"]["zero_mass"]),
            check("best_noncontrol", saved["best_noncontrol"])]
    for row in rows:
        assert row["absolute_score_difference"] < 1e-10
    result = {"method": "endpoint-cell marking independent of recursive search sweep", "checks": rows}
    (ROOT / "verification.json").write_text(json.dumps(result, indent=2, ensure_ascii=False) + "\n")
    print(json.dumps(result, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
