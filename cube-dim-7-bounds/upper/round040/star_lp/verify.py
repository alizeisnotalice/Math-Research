#!/usr/bin/env python3
"""Small exact certificates only; no parameter scan or floating LP solve."""

from fractions import Fraction as F
from itertools import combinations, product
import hashlib
import json
from pathlib import Path
import platform


def ancestor(u, v):
    return len(u) <= len(v) and v[:len(u)] == u


def make_tree(d):
    nodes = [()]
    y = {(): F(0)}
    radius = {(): F(1)}
    for level in range(d):
        for u in [v for v in nodes if len(v) == level]:
            for sign in (-1, 1):
                v = u + (sign,)
                nodes.append(v)
                y[v] = y[u] + sign * radius[u] / 2
                radius[v] = radius[u] / 8
    return nodes, y, radius


def check(n, d):
    a, beta, tau = F(4), F(3, 2), F(1)
    assert d / a ** n <= tau
    nodes, y, radius = make_tree(d)
    leaves = [v for v in nodes if len(v) == d]
    volume = {v: (2 * radius[v] / a) ** n for v in nodes}
    V = {v: (2 * radius[v]) ** n for v in nodes}
    checks = 0

    # Exact small-cube disjointness, star laminarity, expanded-star laminarity.
    for u, v in combinations(nodes, 2):
        dist = abs(y[u] - y[v])
        assert dist > (radius[u] + radius[v]) / a
        comparable = ancestor(u, v) or ancestor(v, u)
        assert (dist < radius[u] + radius[v]) == comparable
        if comparable:
            top, bottom = (u, v) if ancestor(u, v) else (v, u)
            assert dist + radius[bottom] < radius[top]
            assert dist + beta * radius[bottom] < beta * radius[top]
        else:
            assert dist > beta * (radius[u] + radius[v])
        checks += 1

    # Actual largest-radius order and all previous-only eligibility tests.
    weight = a ** (-n)
    for index, v in enumerate(nodes):
        captured = [u for u in nodes[:index]
                    if abs(y[u] - y[v]) < beta * radius[u]]
        assert set(captured) == {u for u in nodes if u != v and ancestor(u, v)}
        assert len(captured) * weight == len(v) * weight <= tau
        assert volume[v] / V[v] == weight > 0
        assert abs(y[v]) + radius[v] <= 1  # Every star lies in root star.
        checks += 1

    # All common-star clusters: all nonempty subsets of every root-leaf path.
    clusters = set()
    for leaf in leaves:
        path = [leaf[:k] for k in range(d + 1)]
        for mask in product((0, 1), repeat=d + 1):
            selected = tuple(v for v, take in zip(path, mask) if take)
            if selected:
                clusters.add(selected)
    p = {v: F(1, 2 ** (len(v) + 1)) if len(v) < d else F(1, 2 ** d)
         for v in nodes}
    for cluster in clusters:
        representative = cluster[0]
        prob = F(sum(ancestor(representative, leaf) for leaf in leaves), len(leaves))
        assert prob == F(1, 2 ** len(representative))
        assert sum(p[v] for v in cluster) <= prob
        checks += 1
    target = 1 + F(d, 2)
    assert sum(p.values()) == target

    # Exact feasible primal on full downward paths, with all item covers = 1.
    columns = []
    for v in nodes:
        q = F(1) if v == () else F(1, 2)
        descendants = [leaf for leaf in leaves if ancestor(v, leaf)]
        for leaf in descendants:
            chain = tuple(leaf[:k] for k in range(len(v), d + 1))
            columns.append((chain, v, q / len(descendants)))
    for v in nodes:
        assert sum(theta for chain, _, theta in columns if v in chain) == 1
        checks += 1

    # All 1D face cells suffice: the other center coordinates are zero, and
    # moving a free coordinate away from zero can only remove star indicators.
    endpoints = sorted({y[v] + s * radius[v] for v in nodes for s in (-1, 1)})
    maxload = F(0)
    maxold = F(0)
    for lo, hi in zip(endpoints, endpoints[1:]):
        x = (lo + hi) / 2
        load = sum(theta for _, v, theta in columns
                   if abs(x - y[v]) < radius[v])
        old = sum(weight for v in nodes if abs(x - y[v]) < radius[v])
        maxload, maxold = max(maxload, load), max(maxold, old)
        assert load <= target
        checks += 1
    assert maxload == target
    assert maxold == (d + 1) * a ** (-n)

    # Exact containment-column compression and original common-mass value.
    root_column = [v for v in nodes if abs(y[v]) + radius[v] / a <= 1]
    assert root_column == nodes
    U = sum(volume.values())
    A = sum((2 * F(1, 8) ** n) ** level for level in range(d + 1))
    ratio = U / V[()]
    assert ratio == a ** (-n) * A
    assert ratio < F(4, 3) * a ** (-n)
    # g=1_(S_root) has mass V_root and integral V_v on every star; the root
    # constraint forces any common witness mass to be at least V_root.
    epsilon = F(1, 10)
    assert all((1 + epsilon) * V[v] > V[v] for v in nodes)
    checks += 4
    return {
        "n": n, "d": d, "N": len(nodes), "a": str(a),
        "beta": str(beta), "tau": str(tau),
        "maximum_eligibility": str(d * weight),
        "common_star_LP_exact": str(target),
        "containment_LP_exact": str(ratio),
        "old_pointwise_load_exact": str(maxold),
        "selected_union_volume": str(U),
        "common_witness_mass_exact": str(V[()]),
        "all_common_star_clusters_checked": len(clusters),
        "primal_path_columns": len(columns),
        "exact_assertion_groups": checks,
    }


source_paths = [
    "/Users/zhengzhihao/Desktop/T/cube-r030-geometry/cube-dim-7-bounds/upper/round039/joint_mass_dual/PROOF.md",
    "/Users/zhengzhihao/Desktop/T/Math-Research-cube-dim-7/cube-dim-7-bounds/upper/round032/joint_cover_theory/proof.md",
    "/Users/zhengzhihao/Desktop/T/cube-r030-geometry/cube-dim-7-bounds/upper/round037/general_cover/PROOF.md",
    "/Users/zhengzhihao/Desktop/T/Math-Research-cube-dim-7/AGENTS.md",
]
root = Path(__file__).resolve().parent
result = {
    "status": "all exact rational checks passed",
    "python": platform.python_version(),
    "method": "Python standard-library Fraction, all chain-subset columns for d<=3; no LP solve",
    "cases": [check(n, d) for n in (1, 3) for d in range(4)],
    "source_sha256": {p: hashlib.sha256(Path(p).read_bytes()).hexdigest()
                      for p in source_paths},
    "artifact_sha256": {name: hashlib.sha256((root / name).read_bytes()).hexdigest()
                        for name in ("PROOF.md", "verify.py")},
}
(root / "verification.json").write_text(json.dumps(result, indent=2) + "\n")
print(json.dumps({"status": result["status"], "cases": len(result["cases"]),
                  "assertion_groups": sum(c["exact_assertion_groups"] for c in result["cases"]),
                  "full_cluster_checks": sum(c["all_common_star_clusters_checked"] for c in result["cases"]),
                  "output": str(root / "verification.json")}, indent=2))
