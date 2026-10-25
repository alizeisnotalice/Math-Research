#!/usr/bin/env python3
"""Independent, bounded audit of two six-dimensional H_epsilon anchors.

This deliberately does not import the producer module and does not reuse its
node aggregation.  It reads only the locked 21-cell fan JSON, reconstructs the
signed-permutation action, then compares (i) a direct orbit sum and (ii) an
independently aggregated exact-node sum.  It is not a cover, continuity, g,
face-integral, or C6 certification.
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

import sympy as sp
from flint import arb, ctx
import flint


HERE = Path(__file__).resolve().parent
SOURCE = HERE.parents[2] / "upper/round002"
FAN = SOURCE / "fan" / "exact_cells_results.json"
AUTHOR_RESULTS = SOURCE / "anchors" / "truncated_anchor_results.json"
AUTHOR_BATCH = SOURCE / "anchors" / "truncated_anchor_batch.py"
GROUP_SIZE = 2**6 * math.factorial(6)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def Q(text: str) -> Fraction:
    return Fraction(text)


def aq(value: Fraction) -> arb:
    """Create an Arb ball from an exact rational, at the current precision."""
    return arb(value.numerator) / value.denominator


def load_cones():
    payload = json.loads(FAN.read_text())
    cones = []
    for cell in payload["cells"]:
        rays = tuple(tuple(Q(v) for v in ray) for ray in cell["rays_a0_equals_one"])
        rates = tuple(Q(v) for v in cell["rates"])
        alpha = Q(cell["determinant"]) / math.prod(rates)
        assert len(rays) == len(rates) == 6
        assert all(rate > 0 for rate in rates)
        assert alpha > 0
        cones.append((rays, rates, alpha))
    assert len(cones) == 21
    # The JSON calls this the full signed-permutation value; the 21 supplied
    # cells are one fundamental section, hence their coefficient sum is 1/|W|.
    full_integral = Q(payload["full_exp_minus_symbol_integral"])
    assert GROUP_SIZE * sum((alpha for _, _, alpha in cones), Fraction()) == full_integral
    return cones, full_integral


CONES, FULL_EXP_INTEGRAL = load_cones()


def all_action_images(point: tuple[Fraction, ...]) -> Counter[tuple[Fraction, ...]]:
    """Exact multiset of all 46080 signed-permutation images."""
    result: Counter[tuple[Fraction, ...]] = Counter()
    for perm in itertools.permutations(range(6)):
        moved = tuple(point[i] for i in perm)
        for signs in itertools.product((-1, 1), repeat=6):
            result[tuple(sign * coordinate for sign, coordinate in zip(signs, moved, strict=True))] += 1
    assert sum(result.values()) == GROUP_SIZE
    return result


def nodes_for(point: tuple[Fraction, ...], rays, rates) -> tuple[Fraction, ...]:
    return tuple(sorted(sum((x * r for x, r in zip(point, ray, strict=True)), Fraction()) / q
                        for ray, q in zip(rays, rates, strict=True)))


def derivative(node: Fraction, epsilon: Fraction, order: int) -> arb:
    x, e = aq(node), aq(epsilon)
    if order == 0:
        return (x / e).atan()
    d = e * e + x * x
    if order == 1:
        return e / d
    if order == 2:
        return -2 * e * x / d**2
    if order == 3:
        return 2 * e * (3 * x * x - e * e) / d**3
    if order == 4:
        return 24 * e * x * (e * e - x * x) / d**4
    if order == 5:
        return 24 * e * (e**4 - 10 * e * e * x * x + 5 * x**4) / d**5
    raise ValueError(order)


def hermite_dd(nodes: tuple[Fraction, ...], epsilon: Fraction) -> arb:
    """Confluent Newton DD.  Sorting makes equal nodes a contiguous block."""
    assert tuple(sorted(nodes)) == nodes
    table = [[arb(0) for _ in range(6)] for _ in range(6)]
    for i, node in enumerate(nodes):
        table[i][0] = derivative(node, epsilon, 0)
    for order in range(1, 6):
        for i in range(order, 6):
            if nodes[i] == nodes[i - order]:
                table[i][order] = derivative(nodes[i], epsilon, order) / math.factorial(order)
            else:
                table[i][order] = (table[i][order - 1] - table[i - 1][order - 1]) / aq(nodes[i] - nodes[i - order])
    return table[5][5]


def direct_and_aggregated(point: tuple[Fraction, ...], epsilon: Fraction):
    images = all_action_images(point)
    multiplicities = set(images.values())
    assert len(multiplicities) == 1
    stabilizer = next(iter(multiplicities))
    assert len(images) * stabilizer == GROUP_SIZE

    # Direct finite sum: no producer aggregation reused.
    direct = arb(0)
    aggregate: dict[tuple[Fraction, ...], Fraction] = defaultdict(Fraction)
    for image, multiplicity in images.items():
        for rays, rates, alpha in CONES:
            nodes = nodes_for(image, rays, rates)
            coefficient = multiplicity * alpha
            direct += aq(coefficient) * hermite_dd(nodes, epsilon)
            aggregate[nodes] += coefficient
    assert sum(aggregate.values(), Fraction()) == FULL_EXP_INTEGRAL
    merged = arb(0)
    for nodes, coefficient in aggregate.items():
        merged += aq(coefficient) * hermite_dd(nodes, epsilon)
    normalizer = (2 * arb.pi())**6
    direct /= normalizer
    merged /= normalizer
    difference = direct - merged
    assert 0 in difference
    return direct, merged, difference, {
        "orbit_size": len(images),
        "stabilizer": stabilizer,
        "full_action_count": sum(images.values()),
        "direct_terms": len(images) * len(CONES),
        "merged_exact_node_tuples": len(aggregate),
    }


def symbolic_derivative_check() -> bool:
    """Checks the five hard-coded F_e derivatives algebraically, not numerically."""
    x, e = sp.symbols("x e", positive=True, real=True)
    expected = {
        1: e / (e**2 + x**2),
        2: -2 * e * x / (e**2 + x**2)**2,
        3: 2 * e * (3*x**2 - e**2) / (e**2 + x**2)**3,
        4: 24 * e * x * (e**2 - x**2) / (e**2 + x**2)**4,
        5: 24 * e * (e**4 - 10*e**2*x**2 + 5*x**4) / (e**2 + x**2)**5,
    }
    return all(sp.simplify(sp.diff(sp.atan(x/e), x, order) - formula) == 0
               for order, formula in expected.items())


def author_ball(k: int, a: str, epsilon: str) -> arb:
    payload = json.loads(AUTHOR_RESULTS.read_text())
    for row in payload["anchors"]:
        if row["k"] == k and row["a"] == a:
            for value in row["values"]:
                if value["epsilon"] == epsilon:
                    return arb(value["interval"])
    raise KeyError((k, a, epsilon))


def render(value: arb) -> str:
    return value.str(90, more=True)


def main() -> None:
    ctx.prec = 512
    assert symbolic_derivative_check()
    # Case A independently reproduces a non-pilot author point.
    point_a = (Fraction(1),) + (Fraction(3, 4),) * 4 + (Fraction(0),)
    eps_a = Fraction(1, 1000)
    direct_a, merged_a, diff_a, meta_a = direct_and_aggregated(point_a, eps_a)
    declared_a = author_ball(4, "3/4", "1/1000")
    declared_difference = direct_a - declared_a
    assert 0 in declared_difference
    # Case B intentionally uses a fresh epsilon, hence no author result is read.
    point_b = (Fraction(1), Fraction(1, 2), Fraction(1, 2), Fraction(0), Fraction(0), Fraction(0))
    eps_b = Fraction(1, 1500)
    direct_b, merged_b, diff_b, meta_b = direct_and_aggregated(point_b, eps_b)
    payload = {
        "status": "independent_bounded_anchor_audit_pass",
        "scope": "Two finite H_epsilon values only; neither global H/g nor face integral nor C6 is certified.",
        "formula_scope": "Assumes the separately stated finite-cone identity H_epsilon=(2*pi)^-6 sum alpha DD5 atan(b/epsilon); this audit checks its finite exact-node realization.",
        "precision_bits": ctx.prec,
        "software": {"python": platform.python_version(), "python_flint": getattr(flint, "__version__", "unknown"), "sympy": sp.__version__},
        "inputs": {"fan": str(FAN), "fan_sha256": sha256(FAN), "author_batch_sha256": sha256(AUTHOR_BATCH), "author_results_sha256": sha256(AUTHOR_RESULTS)},
        "symbolic_F_e_derivatives_orders_1_to_5": True,
        "reality_conjugation_check": "All b_j, alpha and multiplicities are exact real rationals. The evaluator uses the real DD5[atan(b/epsilon)] identity, which is the real part of the complex cone integral in the stated kernel interface; no unpaired complex term is numerically discarded.",
        "epsilon_check": "Both epsilons are positive exact Fractions. The formula is evaluated for integral from epsilon to infinity, not substituted for the epsilon=0 Green potential.",
        "case_A_reproduced_author_nonpilot_anchor": {
            "point": [str(x) for x in point_a], "epsilon": str(eps_a),
            "direct_interval": render(direct_a), "aggregated_interval": render(merged_a),
            "direct_minus_aggregated": render(diff_a), "contains_zero": bool(0 in diff_a),
            "author_declared_interval": render(declared_a),
            "direct_minus_author_declared": render(declared_difference), "overlaps_author_interval": bool(0 in declared_difference), **meta_a,
        },
        "case_B_fresh_epsilon": {
            "point": [str(x) for x in point_b], "epsilon": str(eps_b),
            "direct_interval": render(direct_b), "aggregated_interval": render(merged_b),
            "direct_minus_aggregated": render(diff_b), "contains_zero": bool(0 in diff_b), **meta_b,
        },
        "interval_rounding_check": "All input nodes/coefficients are Fraction through equality, sorting, orbit multiplicity and aggregation; arithmetic after exact-to-Arb conversion is Arb ball arithmetic at 512 bits. The reported values are enclosures, not displayed floating-point values.",
        "not_checked": ["validity of the 21-cell fan as a full Fourier cone cover", "global convergence as epsilon tends to zero", "continuous face minimum", "g0", "I0", "a weak-(1,1) C6 bound"],
        "source_sha256": sha256(Path(__file__)),
    }
    output = HERE / "independent_anchor_audit_results.json"
    output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(json.dumps(payload, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
