#!/usr/bin/env python3
"""Make the registered, bounded coordinate search around generic point 14.

The body-frame is invariant under coordinate permutations, so the seed is
first put in the sorted positive chamber solely to describe the search.  The
direct evaluator still sums the complete signed-permutation group for every
output point; no symmetry orbit is deleted.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
PARENT = HERE.parent
D = 1_000_000_007
SEED_INDEX = 14
SEED_TAIL = (96_968_636, 467_291_939, 987_118_103, 892_738_545, 807_464_540)
STEP = D // 1024
LEVELS = tuple(range(1, 13))


def main() -> None:
    sorted_tail = tuple(sorted(SEED_TAIL, reverse=True))
    rows = [(D, *sorted_tail)]
    for coordinate in range(5):
        for sign in (-1, 1):
            for level in LEVELS:
                candidate = list(sorted_tail)
                candidate[coordinate] += sign * level * STEP
                assert 0 < candidate[coordinate] < D
                # The maximum total displacement is smaller than every gap
                # relevant to a chamber crossing for this fixed seed.
                assert candidate == sorted(candidate, reverse=True)
                rows.append((D, *candidate))
    assert len(rows) == 121 and len(set(rows)) == 121
    out = HERE / "local_points.txt"
    out.write_text(f"{D} {len(rows)}\n" + "\n".join(" ".join(map(str, row)) for row in rows) + "\n")
    metadata = {
        "status": "registered_local_coordinate_search_input",
        "scope": "121 deterministic sorted-chamber points; not a global-minimum search",
        "seed_index_in_32_point_batch": SEED_INDEX,
        "seed_tail_unsorted": list(SEED_TAIL),
        "seed_tail_sorted": list(sorted_tail),
        "denominator": D,
        "step_numerator": STEP,
        "levels_per_coordinate_and_sign": list(LEVELS),
        "points": len(rows),
        "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "points_sha256": hashlib.sha256(out.read_bytes()).hexdigest(),
        "parent_generic_results_sha256": hashlib.sha256((PARENT / "direct_kernel_generic_results.json").read_bytes()).hexdigest(),
    }
    (HERE / "input_metadata.json").write_text(json.dumps(metadata, indent=2) + "\n")
    print(json.dumps(metadata, indent=2))


if __name__ == "__main__":
    main()
