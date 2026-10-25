#!/usr/bin/env python3
"""Independent finite audit of the delta=(0,3/128,0,0,0) L5 witness.

It imports only the earlier audit's exact Q(u,v) scalar class/root isolator;
it imports neither a producer field constructor nor either producer union
engine.  Boxes and their union are reconstructed below using a fresh
coordinate-segment bitmask sweep.
"""
from __future__ import annotations

from fractions import Fraction as Q
from functools import lru_cache
from importlib.util import module_from_spec, spec_from_file_location
from itertools import product
from math import prod
from pathlib import Path
import hashlib
import json
import time


HERE = Path(__file__).resolve().parent
CORE_PATH = HERE / "002_independent_lower_audit.py"
spec = spec_from_file_location("audit001_scalar_core", CORE_PATH)
assert spec and spec.loader
core = module_from_spec(spec)
spec.loader.exec_module(core)
Alg, u, v, root = core.Alg, core.u, core.v, core.fifth_root_of_2_3_smooth

TARGET = (Q(0), Q(3, 128), Q(0), Q(0), Q(0))
OLD_CERT = Path("/Users/zhengzhihao/Desktop/T/Math-Research-cube5-lower-continue/cube-dim-5-bounds/lower/rectangular_certificate.json")
PRODUCER_RESULT = Path("/Users/zhengzhihao/Desktop/T/Math-Research-cube5-lower-continue/cube-dim-5-bounds/lower_round003/local_perturbations_results.json")
NODE_LIMIT = 200_000


def build_retained_boxes():
    spacing = [
        1 + u**3 / 2,
        (2 + u**4 * v) / 4 + Q(3, 128),
        (1 + u) / 2,
        (1 + u) / 2,
        (1 + u) / 2,
    ]
    raw = []
    for lengths in product(range(1, 5), repeat=5):
        radius = root(prod(lengths)) / 2
        coordinates = []
        for step, length in zip(spacing, lengths):
            if length % 2:
                low, high = Alg(0), min(step / 2, radius - (length - 1) * step / 2)
            else:
                low, high = max(Alg(0), length * step / 2 - radius), step / 2
            if high <= low:
                break
            coordinates.extend((low, high))
        else:
            raw.append((lengths, tuple(coordinates)))
    # The map to a dict deliberately deduplicates identical geometry before the
    # containment test, instead of relying on the producer's list order.
    unique = {box: lengths for lengths, box in raw}
    retained = []
    for box, lengths in unique.items():
        contained = any(
            other != box and all(other[2*j] <= box[2*j] and box[2*j+1] <= other[2*j+1] for j in range(5))
            for other in unique
        )
        if not contained:
            retained.append((lengths, box))
    return tuple(sorted(retained))


def bitmask_union(boxes, deadline):
    """Exact 5D union by coordinate slab masks; no producer union import."""
    per_axis = []
    for axis in range(5):
        cuts = sorted({box[2*axis + side] for box in boxes for side in (0, 1)})
        slabs = []
        for left, right in zip(cuts, cuts[1:]):
            mask = 0
            for index, box in enumerate(boxes):
                if box[2*axis] <= left and right <= box[2*axis + 1]:
                    mask |= 1 << index
            if mask:
                slabs.append((right - left, mask))
        per_axis.append(tuple(slabs))
    nodes = 0

    @lru_cache(maxsize=None)
    def sweep(axis, mask):
        nonlocal nodes
        nodes += 1
        if nodes > NODE_LIMIT:
            raise RuntimeError("independent node cap reached")
        if time.monotonic() > deadline:
            raise TimeoutError("independent 120-second cap reached")
        if axis == 5:
            return Alg(1)
        answer = Alg(0)
        for width, slab_mask in per_axis[axis]:
            overlap = mask & slab_mask
            if overlap:
                answer += width * sweep(axis + 1, overlap)
        return answer

    all_boxes = (1 << len(boxes)) - 1
    return sweep(0, all_boxes), [len(slabs) for slabs in per_axis], nodes


def selected_old_upper():
    payload = json.loads(OLD_CERT.read_text(encoding="utf-8"))
    for row in payload["candidates"]:
        if row["delta1"] == "0" and row["delta2"] == "1/64":
            return Q(row["ratio_interval"][1])
    raise KeyError("locked L5-001 certificate row not found")


def main():
    started = time.monotonic()
    retained = build_retained_boxes()
    boxes = tuple(box for _, box in retained)
    volume, segment_counts, nodes = bitmask_union(boxes, started + 120)
    ratio = 32 * volume
    lower, upper = ratio.interval()
    old_upper = selected_old_upper()
    threshold = Q(19937, 10000)
    producer = json.loads(PRODUCER_RESULT.read_text(encoding="utf-8"))
    row = next(item for item in producer["candidates"] if item["delta"] == [str(x) for x in TARGET])
    assert [str(x) for x in ratio.coeff] == [row["coefficients"].get(f"{i//5},{i%5}", "0") for i in range(25)]
    assert [list(lengths) for lengths, _ in retained] == row["retained_witness_counts"]
    assert lower > old_upper
    assert lower > threshold
    result = {
        "status": "PASS_independent_new_finite_witness_reconstruction",
        "candidate_delta": [str(x) for x in TARGET],
        "shared_component": "only audit/002_independent_lower_audit.py's Q(u,v) scalar class and root-isolation endpoints; no producer field or union engine was imported",
        "retained_witness_count": len(retained),
        "retained_witness_counts": [list(lengths) for lengths, _ in retained],
        "coordinate_segment_counts": segment_counts,
        "bitmask_sweep_nodes": nodes,
        "ratio_interval": [str(lower), str(upper)],
        "ratio_coefficients": {f"{i//5},{i%5}": str(c) for i, c in enumerate(ratio.coeff) if c},
        "old_locked_certificate_upper": str(old_upper),
        "strictly_above_old_upper": True,
        "rational_threshold": str(threshold),
        "strictly_above_19937_over_10000": True,
        "elapsed_seconds": time.monotonic() - started,
        "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "shared_scalar_core_sha256": hashlib.sha256(CORE_PATH.read_bytes()).hexdigest(),
        "old_certificate_sha256": hashlib.sha256(OLD_CERT.read_bytes()).hexdigest(),
        "producer_result_sha256": hashlib.sha256(PRODUCER_RESULT.read_bytes()).hexdigest(),
    }
    output = HERE / "010_independent_new_lower_results.json"
    output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
