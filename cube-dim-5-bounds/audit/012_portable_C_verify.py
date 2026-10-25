#!/usr/bin/env python3
"""Minimal portable verifier for the round-003 C finite lower witness.

Place this file beside 002_independent_lower_audit.py in
cube-dim-5-bounds/audit/.  By default, certificates are read relative to the
enclosing cube-dim-5-bounds project; --project-root permits one explicit
project-root override for an external replay.  No producer field/union module
is imported.
"""
from __future__ import annotations

from argparse import ArgumentParser
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
SCALAR = HERE / "002_independent_lower_audit.py"
spec = spec_from_file_location("portable_audit_scalar", SCALAR)
assert spec and spec.loader
scalar = module_from_spec(spec)
spec.loader.exec_module(scalar)
Alg, u, v, root = scalar.Alg, scalar.u, scalar.v, scalar.fifth_root_of_2_3_smooth
LIMIT = 200_000


def retained_C_boxes():
    spacings = [1 + u**3 / 2] + [(1 + u) / 2] * 4
    raw = []
    for lengths in product(range(1, 5), repeat=5):
        radius = root(prod(lengths)) / 2
        endpoints = []
        for step, length in zip(spacings, lengths):
            if length % 2:
                low, high = Alg(0), min(step / 2, radius - (length - 1) * step / 2)
            else:
                low, high = max(Alg(0), length * step / 2 - radius), step / 2
            if high <= low:
                break
            endpoints += [low, high]
        else:
            raw.append((lengths, tuple(endpoints)))
    unique = {box: lengths for lengths, box in raw}
    return tuple(sorted((lengths, box) for box, lengths in unique.items() if not any(
        other != box and all(other[2*j] <= box[2*j] and box[2*j+1] <= other[2*j+1] for j in range(5))
        for other in unique)))


def bitmask_union(boxes, deadline):
    """Fresh exact slab-mask union, included here for standalone replay."""
    axes = []
    for axis in range(5):
        cuts = sorted({box[2*axis + side] for box in boxes for side in (0, 1)})
        slabs = []
        for low, high in zip(cuts, cuts[1:]):
            mask = 0
            for index, box in enumerate(boxes):
                if box[2*axis] <= low and high <= box[2*axis + 1]:
                    mask |= 1 << index
            if mask:
                slabs.append((high - low, mask))
        axes.append(tuple(slabs))
    nodes = 0

    @lru_cache(maxsize=None)
    def sweep(axis, mask):
        nonlocal nodes
        nodes += 1
        if nodes > LIMIT:
            raise RuntimeError("node cap reached")
        if time.monotonic() > deadline:
            raise TimeoutError("120-second cap reached")
        if axis == 5:
            return Alg(1)
        answer = Alg(0)
        for width, slab_mask in axes[axis]:
            overlap = mask & slab_mask
            if overlap:
                answer += width * sweep(axis + 1, overlap)
        return answer

    return sweep(0, (1 << len(boxes)) - 1), [len(axis) for axis in axes], nodes


def old_001_upper(payload):
    for candidate in payload["candidates"]:
        if candidate["delta1"] == "0" and candidate["delta2"] == "1/64":
            return Q(candidate["ratio_interval"][1])
    raise KeyError("round-001 locked candidate")


def main():
    parser = ArgumentParser()
    parser.add_argument("--project-root", type=Path, default=HERE.parent,
                        help="directory containing lower/ and lower_round003/")
    args = parser.parse_args()
    project = args.project_root.resolve()
    round3 = project / "lower_round003"
    result_path, lock_c_path, lock_b_path = (round3 / "a2_equals_a3_results.json",
                                              round3 / "lock_c_a2_equals_a3.json",
                                              round3 / "lock_b_1_32.json")
    old_path = project / "lower" / "rectangular_certificate.json"
    started = time.monotonic()
    producer, lock_c, lock_b, old = (json.loads(path.read_text()) for path in
                                     (result_path, lock_c_path, lock_b_path, old_path))
    assert lock_c["spacing"] == ["1+u^3/2"] + ["(1+u)/2"] * 4
    assert lock_c["delta2_expression"] == "(2*u-u^4*v)/4"
    retained = retained_C_boxes()
    union, segments, nodes = bitmask_union(tuple(box for _, box in retained), started + 120)
    ratio = 32 * union
    low, high = ratio.interval()
    coefficient_list = [producer["coefficients"].get(f"{i//5},{i%5}", "0") for i in range(25)]
    assert [str(c) for c in ratio.coeff] == coefficient_list
    assert len(retained) == producer["retained_witness_count"] == lock_c["retained_witness_count"] == 42
    lock_low, lock_high = map(Q, lock_c["exact_interval"])
    b_upper = Q(lock_b["exact_interval"][1])
    old_upper = old_001_upper(old)
    threshold = Q(19943, 10000)
    assert lock_low <= low <= high <= lock_high
    assert low > b_upper and low > old_upper and low > threshold
    output = {
        "status": "PASS_portable_C_finite_replay",
        "project_relative_inputs": ["lower_round003/a2_equals_a3_results.json", "lower_round003/lock_c_a2_equals_a3.json", "lower_round003/lock_b_1_32.json", "lower/rectangular_certificate.json"],
        "spacing": lock_c["spacing"],
        "delta2_expression": lock_c["delta2_expression"],
        "enumerated_length_vectors": 1024,
        "retained_witness_counts": [list(x) for x, _ in retained],
        "coordinate_segment_counts": segments,
        "bitmask_sweep_nodes": nodes,
        "ratio_interval": [str(low), str(high)],
        "audit_interval_within_C_lock": True,
        "strictly_above_B_upper": True,
        "strictly_above_001_upper": True,
        "strictly_above_19943_over_10000": True,
        "elapsed_seconds": time.monotonic() - started,
        "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "shared_audit_scalar_sha256": hashlib.sha256(SCALAR.read_bytes()).hexdigest(),
        "C_result_sha256": hashlib.sha256(result_path.read_bytes()).hexdigest(),
        "C_lock_sha256": hashlib.sha256(lock_c_path.read_bytes()).hexdigest(),
        "B_lock_sha256": hashlib.sha256(lock_b_path.read_bytes()).hexdigest(),
        "old_001_sha256": hashlib.sha256(old_path.read_bytes()).hexdigest(),
        "independence_limit": "shares only audit/002 scalar/root code; all C geometry and union code are in this file; no producer field/union import",
    }
    (HERE / "012_portable_C_verify_results.json").write_text(json.dumps(output, indent=2) + "\n")
    print(json.dumps(output, indent=2))


if __name__ == "__main__":
    main()
