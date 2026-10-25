#!/usr/bin/env python3
"""Finite audit of the exact orbit compression used for one anchor."""
from __future__ import annotations

import hashlib
import importlib.util
import itertools
import json
from collections import Counter
from fractions import Fraction
from pathlib import Path

HERE = Path(__file__).resolve().parent
MODULE_PATH = HERE / "truncated_anchor_batch.py"
spec = importlib.util.spec_from_file_location("batch", MODULE_PATH)
batch = importlib.util.module_from_spec(spec)
assert spec.loader is not None
spec.loader.exec_module(batch)


def all_group_images(point):
    images = Counter()
    for permutation in itertools.permutations(range(6)):
        moved = tuple(point[index] for index in permutation)
        for signs in itertools.product((-1, 1), repeat=6):
            images[tuple(Fraction(sign) * coordinate for sign, coordinate in zip(signs, moved, strict=True))] += 1
    return images


def main():
    batch.ctx.prec = 384
    point = batch.anchor(1, Fraction(1, 4))
    images = all_group_images(point)
    orbit = batch.signed_permutation_orbit(point)
    assert set(images) == set(orbit)
    multiplicities = set(images.values())
    assert len(multiplicities) == 1
    stabilizer = multiplicities.pop()
    assert len(images) * stabilizer == batch.GROUP_SIZE
    coefficients, orbit_size, compressed_stabilizer = batch.aggregate_nodes(point)
    assert orbit_size == len(images) and compressed_stabilizer == stabilizer
    assert sum(coefficients.values(), Fraction()) == batch.FULL_EXP_INTEGRAL
    rows = []
    for epsilon in (Fraction(1, 1000), Fraction(1, 2000)):
        unmerged = batch.arb(0)
        for image in orbit:
            for rays, rates, alpha in batch.BASE_CONES:
                nodes = tuple(sorted(sum((coordinate * ray_coordinate for coordinate, ray_coordinate in zip(image, ray, strict=True)), Fraction()) / rate for ray, rate in zip(rays, rates, strict=True)))
                unmerged += batch.arbq(stabilizer * alpha) * batch.hermite_dd(nodes, epsilon)
        unmerged /= (2 * batch.arb.pi())**6
        merged, metadata = batch.h_epsilon(point, epsilon)
        difference = merged - unmerged
        assert 0 in difference
        rows.append({"epsilon": str(epsilon), "unmerged_orbit_terms": len(orbit) * len(batch.BASE_CONES), "merged_exact_node_tuples": metadata["merged_exact_node_tuples"], "difference": difference.str(70, more=True)})
    output = {
        "status": "PASS",
        "anchor": [str(item) for item in point],
        "all_group_elements": batch.GROUP_SIZE,
        "distinct_orbit_points": len(images),
        "uniform_stabilizer": stabilizer,
        "coefficient_sum_equals_full_peak_integral": str(batch.FULL_EXP_INTEGRAL),
        "epsilon_checks": rows,
        "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "batch_source_sha256": hashlib.sha256(MODULE_PATH.read_bytes()).hexdigest(),
        "scope": "checks finite group orbit and exact node aggregation for one anchor; not an independent Green-potential proof",
    }
    (HERE / "verify_anchor_aggregation_results.json").write_text(json.dumps(output, indent=2) + "\n")
    print(json.dumps(output, indent=2))


if __name__ == "__main__":
    main()
