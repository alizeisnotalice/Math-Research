#!/usr/bin/env python3
"""Strict rational necessary-span pruning for the R6-003 spacing family."""
from fractions import Fraction as Q
from itertools import product
from pathlib import Path
import hashlib
import json

N = 6
D = 10**35
HERE = Path(__file__).resolve().parent


def floor_root(x, n):
    lo, hi = 0, 1
    while hi**n <= x:
        hi *= 2
    while hi-lo > 1:
        mid = (lo+hi)//2
        if mid**n <= x:
            lo = mid
        else:
            hi = mid
    return lo


def sixth_interval(integer):
    q = floor_root(integer*D**6, 6)
    return Q(q, D), Q(q+1, D)


def spacing_intervals():
    u0, u1 = sixth_interval(2)
    v0, v1 = sixth_interval(3)
    a1 = (1+u0**4/2, 1+u1**4/2)
    # R6-003 breakpoint spacing, all monomials positive.
    a2 = ((u0**3+u0**4*v0**3)/4, (u1**3+u1**4*v1**3)/4)
    tail = ((1+u0)/2, (1+u1)/2)
    return (a1, a2, tail, tail, tail, tail)


def classify(m, a):
    stats = {"total": 0, "coordinate_span_excluded": 0,
             "product_only_excluded": 0, "not_excluded": 0,
             "certified_all_slab_spans_strict": 0,
             "new_shell_total": 0, "new_shell_not_excluded": 0,
             "all_lengths_at_least_2_total": 0,
             "all_lengths_at_least_2_not_excluded": 0,
             "new_shell_survivors": [],
             "examples": {"span": [], "product": [], "survive": []}}
    for ell in product(range(1, m+1), repeat=N):
        stats["total"] += 1
        shell = max(ell) == m
        if shell:
            stats["new_shell_total"] += 1
        all2 = min(ell) >= 2
        if all2:
            stats["all_lengths_at_least_2_total"] += 1
        P = 1
        for e in ell:
            P *= e
        root_lo, root_hi = sixth_interval(P)
        # Necessary: (ell_j-1)*a_j < P^(1/6) for every j.
        failing = [j for j in range(N) if (ell[j]-1)*a[j][0] >= root_hi]
        if failing:
            stats["coordinate_span_excluded"] += 1
            if len(stats["examples"]["span"]) < 4:
                stats["examples"]["span"].append({"ell": list(ell), "coordinate": failing[0]+1})
            continue
        # Multiplication of all six strict span inequalities gives
        # prod(a_j)*prod((ell_j-1)/ell_j) < 1.  Use its lower interval end
        # to rigorously reject a vector not already rejected above.
        geom_lo = Q(1)
        for e, (alo, _) in zip(ell, a):
            geom_lo *= alo * Q(e-1, e)
        if geom_lo >= 1:
            stats["product_only_excluded"] += 1
            if len(stats["examples"]["product"]) < 4:
                stats["examples"]["product"].append({"ell": list(ell)})
            continue
        stats["not_excluded"] += 1
        if shell:
            stats["new_shell_not_excluded"] += 1
            stats["new_shell_survivors"].append(list(ell))
        if all2:
            stats["all_lengths_at_least_2_not_excluded"] += 1
        if all((e-1)*ahi < root_lo for e, (_, ahi) in zip(ell, a)):
            stats["certified_all_slab_spans_strict"] += 1
        if len(stats["examples"]["survive"]) < 8:
            stats["examples"]["survive"].append({"ell": list(ell), "P": P})
    return stats


def main():
    a = spacing_intervals()
    output = {
        "status": "strict_necessary_pruning_only",
        "dimension": 6, "root_interval_denominator": D,
        "spacing": [{"lower": str(lo), "upper": str(hi)} for lo, hi in a],
        "necessary_conditions": [
            "(ell_j-1)*a_j < (prod ell)^(1/6) for each j",
            "prod_j a_j * prod_j ((ell_j-1)/ell_j) < 1",
        ],
        "m5": classify(5, a), "m6": classify(6, a),
        "scope_warning": "not_excluded means only that these necessary tests did not rule out the length vector; it is not a nonempty-box, union, containment, or C6 claim.",
        "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    }
    (HERE/'span_prune_results.json').write_text(json.dumps(output, indent=2)+'\n')
    print(json.dumps(output, indent=2))


if __name__ == '__main__':
    main()
