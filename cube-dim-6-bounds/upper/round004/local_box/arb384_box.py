#!/usr/bin/env python3
"""Arb384 complete-group H_{1/100} values for the five local-box centers.

The Counter reduction is exact: it groups identical transformed coordinate
tuples, retaining their number of occurrences in S_6 x {+/-}^6.
"""
from __future__ import annotations

from collections import Counter
from functools import cmp_to_key
import hashlib
import itertools
import json
import sys
from pathlib import Path
from time import perf_counter

from flint import acb, arb, ctx

HERE = Path(__file__).resolve().parent
FAN = HERE / "../../round003/direct_kernel/fan_integer_input.txt"


def parse_fan():
    tokens = iter(FAN.read_text().split())
    n = int(next(tokens)); cones = []
    for _ in range(n):
        det = int(next(tokens)); rays = []
        for _ in range(6):
            q = int(next(tokens)); rays.append((q, tuple(int(next(tokens)) for _ in range(6))))
        cones.append((det, rays))
    return cones


def cmp_node(a, b): return (a[0]*b[1] > b[0]*a[1]) - (a[0]*b[1] < b[0]*a[1])
def equal(a, b): return a[0]*b[1] == b[0]*a[1]
def aval(a): return arb(a[0]) / a[1]


def jet(x, epsilon, order):
    if order == 0: return (x/epsilon).atan()
    ii = acb(0, 1)
    plus, minus = acb(epsilon, x), acb(epsilon, -x)
    value = (ii**order/plus**order - (-ii)**order/minus**order) / (2*order*ii)
    return ((1 if order & 1 else -1) * value).real


def dd5(nodes, epsilon):
    z = sorted(nodes, key=cmp_to_key(cmp_node)); vals = [aval(v) for v in z]
    table = [[arb(0) for _ in range(6)] for _ in range(6)]
    for i in range(6): table[i][0] = jet(vals[i], epsilon, 0)
    for order in range(1, 6):
        for i in range(order, 6):
            table[i][order] = jet(vals[i], epsilon, order) if equal(z[i],z[i-order]) else (table[i][order-1]-table[i-1][order-1])/(vals[i]-vals[i-order])
    return table[5][5]


def orbit_counts(x):
    out = Counter()
    for perm in itertools.permutations(range(6)):
        for mask in range(64): out[tuple((-1 if (mask>>i)&1 else 1)*x[perm[i]] for i in range(6))] += 1
    assert sum(out.values()) == 46_080
    return out


def main():
    requested = [int(v) for v in sys.argv[1:]] or [0]
    raw = iter((HERE / "box_centers.txt").read_text().split()); den, n = int(next(raw)), int(next(raw))
    points = [tuple(int(next(raw)) for _ in range(6)) for _ in range(n)]
    assert all(0 <= k < n for k in requested)
    cones = parse_fan(); ctx.prec = 384; epsilon = arb(1)/100
    result_path = HERE / "arb384_results.json"
    prior = {}
    if result_path.exists():
        old = json.loads(result_path.read_text())
        prior = {row["index"]: row for row in old.get("outputs", [])}
    outputs=[]
    for index in requested:
        started = perf_counter(); orbit = orbit_counts(points[index]); total=arb(0); terms=0; equal_events=0
        for y, multiplicity in orbit.items():
            for det, rays in cones:
                nodes=[]; alpha=arb(det)
                for q, ray in rays:
                    nodes.append((sum(y[i]*ray[i] for i in range(6)), den*q)); alpha/=q
                ordered=sorted(nodes,key=cmp_to_key(cmp_node))
                equal_events += multiplicity * sum(equal(ordered[i-1],ordered[i]) for i in range(1,6))
                total += multiplicity * alpha * dd5(nodes,epsilon); terms += multiplicity
        result=total/(2*arb.pi())**6
        outputs.append({
            "status":"PASS_ARB384_COMPLETE_GROUP_H_EPSILON", "index":index,
            "point_numerators":list(points[index]), "denominator":den, "epsilon":"1/100",
            "raw_group_images":46080, "unique_group_images":len(orbit), "multiplicity_sum":sum(orbit.values()),
            "cone_terms_with_multiplicity":terms, "equal_adjacent_events_with_multiplicity":equal_events,
            "H_epsilon_arb_384":result.str(95,more=True), "elapsed_seconds":perf_counter()-started,
        })
    prior.update({row["index"]: row for row in outputs})
    payload={"status":"PASS_ARB384_LOCAL_BOX_CENTER_BATCH", "requested_indices":requested, "outputs":[prior[k] for k in sorted(prior)],
      "source_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
      "centers_sha256":hashlib.sha256((HERE/"box_centers.txt").read_bytes()).hexdigest(),
      "fan_sha256":hashlib.sha256(FAN.read_bytes()).hexdigest()}
    result_path.write_text(json.dumps(payload,indent=2)+"\n")
    print(json.dumps(payload,indent=2))


if __name__ == "__main__": main()
