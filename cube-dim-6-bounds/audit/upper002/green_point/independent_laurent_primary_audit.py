#!/usr/bin/env python3
"""Bounded independent four-jet audit of one Laurent Green-point certificate.

Only the locked fan JSON and the producer's JSON *claim* are read.  No
producer function, divided-difference routine, or node aggregation is
imported.  The finite result is a pointwise value and hence only an upper
bound on the cube infimum g; it is not a cost or C6 certificate.
"""
from __future__ import annotations

import hashlib
import itertools
import json
import math
import platform
from collections import Counter, defaultdict
from fractions import Fraction
from pathlib import Path

from flint import arb, ctx
import flint


HERE = Path(__file__).resolve().parent
SOURCE = HERE.parents[2] / "upper/round002"
FAN = SOURCE / "fan" / "exact_cells_results.json"
AUTHOR = SOURCE / "green_point" / "laurent_points_results.json"
AUTHOR_CODE = SOURCE / "green_point" / "laurent_points.py"
ORDER = (-5, -3, -1, 0)  # coefficients of eps^-5, eps^-3, eps^-1, pi
ZERO = (Fraction(),) * 4
GROUP_SIZE = 2**6 * math.factorial(6)


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def Q(s: str) -> Fraction:
    return Fraction(s)


def vadd(x, y):
    return tuple(a + b for a, b in zip(x, y, strict=True))


def vsub(x, y):
    return tuple(a - b for a, b in zip(x, y, strict=True))


def vscale(x, a: Fraction):
    return tuple(a * b for b in x)


def node_limit_jet(node: Fraction, derivative_order: int):
    """Four exact limiting Laurent coefficients of F^(m)(node)/m!.

    The constant coordinate is the coefficient of pi, not the real constant.
    The nonzero-node omitted remainder is O(epsilon), and is handled only
    after the finite complete sum; it is not asserted cone by cone at epsilon=0.
    """
    if node:
        if derivative_order == 0:
            return (Fraction(), Fraction(), Fraction(), Fraction(1 if node > 0 else -1, 2))
        return ZERO
    if derivative_order == 1:
        return (Fraction(), Fraction(), Fraction(1), Fraction())
    if derivative_order == 3:
        return (Fraction(), Fraction(-1, 3), Fraction(), Fraction())
    if derivative_order == 5:
        return (Fraction(1, 5), Fraction(), Fraction(), Fraction())
    return ZERO


def fourjet_dd5(nodes: tuple[Fraction, ...]):
    """Independent confluent Newton recurrence, over the four Laurent jets."""
    assert len(nodes) == 6 and tuple(sorted(nodes)) == nodes
    table = [[ZERO for _ in range(6)] for _ in range(6)]
    for i, node in enumerate(nodes):
        table[i][0] = node_limit_jet(node, 0)
    for depth in range(1, 6):
        for i in range(depth, 6):
            if nodes[i] == nodes[i - depth]:
                table[i][depth] = node_limit_jet(nodes[i], depth)
            else:
                table[i][depth] = vscale(vsub(table[i][depth - 1], table[i - 1][depth - 1]),
                                         Fraction(1, nodes[i] - nodes[i - depth]))
    return table[5][5]


def load_fan():
    doc = json.loads(FAN.read_text())
    cones = []
    for cell in doc["cells"]:
        rays = tuple(tuple(Q(v) for v in ray) for ray in cell["rays_a0_equals_one"])
        rates = tuple(Q(v) for v in cell["rates"])
        coeff = Q(cell["determinant"]) / math.prod(rates)
        assert len(rays) == len(rates) == 6
        assert all(rate > 0 for rate in rates) and coeff > 0
        cones.append((rays, rates, coeff))
    full = Q(doc["full_exp_minus_symbol_integral"])
    assert len(cones) == 21
    assert GROUP_SIZE * sum((coefficient for _, _, coefficient in cones), Fraction()) == full
    return cones, full


CONES, FULL = load_fan()


def images_with_multiplicity(point):
    result = Counter()
    for p in itertools.permutations(range(6)):
        moved = tuple(point[i] for i in p)
        for signs in itertools.product((-1, 1), repeat=6):
            result[tuple(sign * x for sign, x in zip(signs, moved, strict=True))] += 1
    assert sum(result.values()) == GROUP_SIZE
    assert len(set(result.values())) == 1
    return result


def nodes_for(point, rays, rates):
    return tuple(sorted(sum((x * r for x, r in zip(point, ray, strict=True)), Fraction()) / q
                        for ray, q in zip(rays, rates, strict=True)))


def audit_primary():
    point = (Fraction(1),) + (Fraction(3, 4),) * 4 + (Fraction(),)
    images = images_with_multiplicity(point)
    stabilizer = next(iter(images.values()))
    aggregate = defaultdict(Fraction)
    for image, mult in images.items():
        for rays, rates, alpha in CONES:
            aggregate[nodes_for(image, rays, rates)] += mult * alpha
    assert sum(aggregate.values(), Fraction()) == FULL
    total = ZERO
    for nodes, coefficient in aggregate.items():
        total = vadd(total, vscale(fourjet_dd5(nodes), coefficient))
    e5, e3, e1, pi_coeff = total
    assert e5 == e3 == e1 == 0
    q_over_pi5 = pi_coeff / 64
    # Check the author record only after reconstructing the value.
    author = json.loads(AUTHOR.read_text())["rows"][0]
    assert author["point"] == [str(x) for x in point]
    assert Q(author["coefficients"]["pi"]) == pi_coeff
    assert author["coefficients"]["epsilon^-1"] == author["coefficients"]["epsilon^-3"] == author["coefficients"]["epsilon^-5"] == "0"
    claimed_q = author["G_exact"].removeprefix("(").removesuffix(")/pi^5")
    assert Q(claimed_q) == q_over_pi5
    ctx.prec = 512
    value = (arb(q_over_pi5.numerator) / q_over_pi5.denominator) / arb.pi()**5
    claimed_ball = arb(author["G_arb_384"])
    difference = value - claimed_ball
    assert 0 in difference
    return {
        "point": [str(x) for x in point], "orbit_size": len(images), "stabilizer": stabilizer,
        "all_group_actions": sum(images.values()), "direct_preaggregation_terms": len(images) * len(CONES),
        "merged_exact_node_tuples": len(aggregate),
        "zero_node_tuples": sum(Fraction() in nodes for nodes in aggregate),
        "laurent_coefficients": {"epsilon^-5": str(e5), "epsilon^-3": str(e3), "epsilon^-1": str(e1), "pi": str(pi_coeff)},
        "G_exact": f"({q_over_pi5})/pi^5", "G_arb_512": value.str(90, more=True),
        "G_arb_512_minus_author_arb_384": difference.str(90, more=True),
        "overlaps_author_ball": bool(0 in difference),
    }


def main():
    row = audit_primary()
    payload = {
        "status": "PASS_primary_exact_laurent_point_independent_audit",
        "scope": "one non-S2 point only; proves finite exact G(point) conditional on the locked finite-cone identity, hence g <= G(point), but proves neither a cost bound nor C6",
        "fourjet_order": list(ORDER),
        "finite_remainder_reason": "For each fixed nonzero node, the discarded part of F_e and each derivative jet is O(epsilon); after exact principal-part cancellation this is a finite sum, so it tends to zero. No cone-wise deletion of a zero-node principal part is performed.",
        "inputs": {"fan_sha256": sha(FAN), "producer_source_sha256": sha(AUTHOR_CODE), "producer_results_sha256": sha(AUTHOR)},
        "software": {"python": platform.python_version(), "python_flint": getattr(flint, "__version__", "unknown")},
        "precision_bits_for_display_enclosure": 512,
        "row": row,
        "not_certified": ["a cube minimum equality", "a lower bound for g", "the face integral I", "the model cost", "any C6 upper or lower bound"],
        "source_sha256": sha(Path(__file__)),
    }
    out = HERE / "independent_laurent_primary_audit_results.json"
    out.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(json.dumps(payload, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
