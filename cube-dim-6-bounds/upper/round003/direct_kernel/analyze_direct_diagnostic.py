#!/usr/bin/env python3
"""Summarize the fixed 32-point direct-kernel diagnostic without certification claims."""
from __future__ import annotations

import hashlib
import json
from fractions import Fraction
from itertools import combinations
from pathlib import Path

import mpmath as mp

HERE = Path(__file__).resolve().parent
RESULT = HERE / "direct_kernel_generic_results.json"


def dist2_s2(numerators, denominator):
    x = tuple(Fraction(v, denominator) for v in numerators)
    best = None
    # On this positive face x0=1, the value 1 is one of the two absolute
    # levels.  Its index must be assigned to that level in a nearest S2 point.
    for size in range(1, 7):
        for tail in combinations(range(1, 6), size - 1):
            ones = {0, *tail}
            rest = [i for i in range(6) if i not in ones]
            if rest:
                a = sum((x[i] for i in rest), Fraction()) / len(rest)
                squared = sum(((x[i] - 1) ** 2 for i in ones), Fraction()) + sum(((x[i] - a) ** 2 for i in rest), Fraction())
            else:
                squared = sum(((x[i] - 1) ** 2 for i in ones), Fraction())
            best = squared if best is None or squared < best else best
    return best


def main():
    data = json.loads(RESULT.read_text())
    rows = data["points"]
    assert len(rows) == 32 and all(row["status"] == "GENERIC_OK" and row["zero_node_terms"] == 0 and row["repeated_node_pairs"] == 0 for row in rows)
    mp.mp.dps = 80
    values = [mp.mpf(row["G_direct_long_double"]) for row in rows]
    q = mp.mpf(6) / 5
    mean_q = sum(value**q for value in values) / len(values)
    index = min(range(len(values)), key=values.__getitem__)
    minimum = values[index]
    cost = mp.power(6, q) * mean_q / minimum**q
    distances = [dist2_s2(row["numerators"], row["denominator"]) for row in rows]
    output = {
        "status": "numerical_route_selection_only",
        "points": len(rows),
        "q": "6/5",
        "sample_mean_G_to_q": mp.nstr(mean_q, 40),
        "sample_min_G": mp.nstr(minimum, 40),
        "sample_min_index": index,
        "sample_min_numerators": rows[index]["numerators"],
        "formal_sample_cost": mp.nstr(cost, 40),
        "max_cancellation_ratio": mp.nstr(max(mp.mpf(row["cancellation_ratio"]) for row in rows), 30),
        "total_elapsed_seconds": mp.nstr(sum(mp.mpf(row["elapsed_seconds"]) for row in rows), 30),
        "max_point_elapsed_seconds": mp.nstr(max(mp.mpf(row["elapsed_seconds"]) for row in rows), 30),
        "nearest_S2_distance_squared_exact_per_point": [str(value) for value in distances],
        "nearest_S2_distance_min": mp.nstr(mp.sqrt(min(distances)), 30),
        "nearest_S2_distance_at_minpoint": mp.nstr(mp.sqrt(distances[index]), 30),
        "bias": "Fixed 32-point Latin-hypercube sample excludes face boundaries and supplies neither integration error nor a global minimum certificate; proximity to S2 is reported only as a diagnostic.",
        "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "input_sha256": hashlib.sha256(RESULT.read_bytes()).hexdigest(),
    }
    (HERE / "direct_kernel_diagnostic_summary.json").write_text(json.dumps(output, indent=2) + "\n")
    print(json.dumps(output, indent=2))


if __name__ == "__main__": main()
