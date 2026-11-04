#!/usr/bin/env python3
"""Regression controls for contiguous atomic blocks across floating scales.

This test intentionally retains a local copy of the former prefix-difference
method only to demonstrate its failure.  Production code must use
``multiscale_search.block_halo`` instead.
"""
import json
import math

from multiscale_search import block_halo, merge


def old_prefix_halo(atoms):
    """The retired algorithm, used only to expose loss of a micro block."""
    atoms = sorted(atoms)
    x = [point for point, _ in atoms]
    prefix = [0.0]
    for _, mass in atoms:
        prefix.append(prefix[-1] + mass)
    intervals = []
    for left in range(len(atoms)):
        for right in range(left, len(atoms)):
            mass = prefix[right + 1] - prefix[left]
            a, b = x[right] - mass / 2, x[left] + mass / 2
            if a < b:
                intervals.append((a, b))
    return merge(intervals)


def all_subsets_halo(atoms):
    """Independent exhaustive reference for the small regression inputs."""
    intervals = []
    for mask in range(1, 1 << len(atoms)):
        selected = [atoms[index] for index in range(len(atoms)) if mask >> index & 1]
        mass = math.fsum(mass for _, mass in selected)
        a = max(point for point, _ in selected) - mass / 2
        b = min(point for point, _ in selected) + mass / 2
        if a < b:
            intervals.append((a, b))
    return merge(intervals)


def has_interval_covering(intervals, target):
    return any(left <= target[0] and target[1] <= right for left, right in intervals)


def main():
    # The outer mass is deliberately much larger than the central-scale mass.
    # Prefix accumulation rounds the latter away before subtraction; direct
    # fsum on [left,right] retains the physical micro block.
    cross_scale = [
        (-1.0, 1.0e-5),
        (2.5e-189, 1.5e-194),
        (5.0e-95, 5.0e-95),
        (1.0, 6.0e-6),
    ]
    expected_micro = (2.5e-95, 7.5e-95)
    retired = old_prefix_halo(cross_scale)
    direct = block_halo(cross_scale)
    exhaustive = all_subsets_halo(cross_scale)
    assert direct == exhaustive
    assert direct != retired
    assert has_interval_covering(direct, expected_micro)
    assert not has_interval_covering(retired, expected_micro)

    # A fixed two-atom configuration must remain unchanged by the repair.
    fixed_two = [(-1.0, 1.0 / 4.0), (1.0, 15.0 / 8.0)]
    assert block_halo(fixed_two) == all_subsets_halo(fixed_two)

    print(json.dumps({
        "status": "PASS",
        "cross_scale_direct_intervals": len(direct),
        "cross_scale_retired_intervals": len(retired),
        "fixed_two_intervals": len(block_halo(fixed_two)),
        "finite_precision_boundary": (
            "math.fsum prevents prefix-subtraction loss, but IEEE-754 endpoints "
            "still round when an interval width is below local ulp scale"
        ),
    }, sort_keys=True))


if __name__ == "__main__":
    main()
