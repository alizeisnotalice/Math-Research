#!/usr/bin/env python3
"""Exact finite check of the positive-face traces of all 120 body planes."""
from __future__ import annotations

import hashlib
import itertools
import json
from pathlib import Path


def main() -> None:
    directions = [(1,) + signs for signs in itertools.product((-1, 1), repeat=4)]
    full_segments: set[tuple[int, ...]] = set()
    point_only: set[tuple[int, ...]] = set()
    segment_pair_count = point_pair_count = empty_pair_count = 0
    for v, w in itertools.combinations(directions, 2):
        J = tuple(j + 2 for j in range(4) if v[j + 1] != w[j + 1])
        equal = tuple(j for j in range(4) if v[j + 1] == w[j + 1])
        # x=a v+(1-a)w has x_j=v_j on equal coordinates.  A negative
        # equal coordinate cannot meet [0,1]^4.
        if any(v[j + 1] < 0 for j in equal):
            empty_pair_count += 1
            continue
        signs_on_J = {v[j + 1] for j in range(4) if j + 2 in J}
        if len(signs_on_J) == 1:
            full_segments.add(J)
            segment_pair_count += 1
        else:
            # Then a-(1-a)=0 is forced, giving precisely the s=0
            # endpoint of the J-segment.
            point_only.add(J)
            point_pair_count += 1
    wanted = {tuple(j + 2 for j in range(4) if mask & (1 << j)) for mask in range(1, 16)}
    assert full_segments == wanted
    assert point_only <= full_segments
    assert segment_pair_count + point_pair_count + empty_pair_count == 120
    result = {
        "status": "PASS_exact_120_pair_positive_face_trace",
        "body_directions": len(directions),
        "unordered_pairs": 120,
        "pair_counts": {"full_segment": segment_pair_count, "point_only": point_pair_count, "empty": empty_pair_count},
        "distinct_full_segments": len(full_segments),
        "point_only_traces_are_segment_endpoints": True,
        "all_nonempty_coordinate_subsets_obtained": True,
        "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    }
    target = Path(__file__).with_name("007_trace_exact_check_results.json")
    target.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
