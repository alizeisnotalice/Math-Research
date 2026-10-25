#!/usr/bin/env python3
"""Portable audit of lower_round007 using only archived project-relative inputs.

Run with the project directory containing lower/, lower_round003/ and
lower_round007/.  It deliberately reads final exactC files only: each embeds
the six exact union certificates needed for the polynomial check.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import sys
from fractions import Fraction as Q
from itertools import combinations, product
from math import comb, prod
from pathlib import Path


REQUIRED = [
    "lower/certify_rectangular.py",
    "lower_round003/a2_equals_a3_results.json",
    "lower_round007/common_tail_event_intervals.py",
] + [f"lower_round007/exactC_event_interval_{i}_results.json" for i in range(8)]


def decode(F, obj):
    a = [Q(0)] * 25
    for key, value in obj.items():
        i, j = map(int, key.split(",")); a[5 * i + j] = Q(value)
    return F(a=a)


def encode(x):
    return {f"{i // 5},{i % 5}": str(v) for i, v in enumerate(x.a) if v}


def lift(F, q): return F(q.numerator) / q.denominator


def fifth_root(n, u, v):
    a = b = 0
    while n % 2 == 0: n //= 2; a += 1
    while n % 3 == 0: n //= 3; b += 1
    assert n == 1
    return u ** a * v ** b


def raw_forms(F, u, v):
    a1, tail = 1 + u ** 3 / 2, (1 + u) / 2
    result = {(F(0), Q(0)), (a1 / 2, Q(0))}
    for L in product(range(1, 5), repeat=5):
        r = fifth_root(prod(L), u, v) / 2
        e = L[0]
        result.add((r - (e - 1) * a1 / 2 if e % 2 else e * a1 / 2 - r, Q(0)))
        for e in L[1:]:
            result.add((tail / 2, Q(1, 2)))
            result.add((r - (e - 1) * tail / 2 if e % 2 else e * tail / 2 - r,
                        Q(-(e - 1), 2) if e % 2 else Q(e, 2)))
    return result


def all_events(F, u, v):
    forms = raw_forms(F, u, v)
    low, high = lift(F, Q(-1, 16)), lift(F, Q(1, 16))
    out = {low.a: low, high.a: high, F(0).a: F(0)}
    for (b, m), (c, n) in combinations(forms, 2):
        if m != n:
            z = (c - b) / (m - n)
            if low < z and z < high: out[z.a] = z
    return forms, tuple(sorted(out.values()))


def interpolate(F, xs, ys):
    ans = [F(0)] * 5
    for i, (x, y) in enumerate(zip(xs[:5], ys[:5])):
        p, den = [F(1)], Q(1)
        for j, z in enumerate(xs[:5]):
            if i == j: continue
            nxt = [F(0)] * (len(p) + 1)
            for k, c in enumerate(p): nxt[k] += c * (-z); nxt[k + 1] += c
            p, den = nxt, den * (x - z)
        for k, c in enumerate(p): ans[k] += y * c / den
    return ans


def value(c, x):
    y = c[-1]
    for q in reversed(c[:-1]): y = y * x + q
    return y


def bernstein4(c, left, right):
    h, power = right - left, [c[0] * 0 for _ in range(5)]
    for i, z in enumerate(c):
        for k in range(i + 1): power[k] += z * Q(comb(i, k)) * left ** (i-k) * h ** k
    return [sum((power[i] * Q(comb(k, i), comb(4, i)) for i in range(k+1)), power[0] * 0)
            for k in range(5)]


def sha(path): return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    p = argparse.ArgumentParser()
    p.add_argument("project_root", type=Path)
    p.add_argument("--output", type=Path, required=True)
    args = p.parse_args(); root = args.project_root.resolve()
    missing = [x for x in REQUIRED if not (root / x).is_file()]
    if missing: raise FileNotFoundError("missing archived inputs: " + ", ".join(missing))
    sys.path.insert(0, str(root / "lower"))
    from certify_rectangular import field
    F, u, v, _, _, _ = field()
    forms, ev = all_events(F, u, v)
    left, right, zero = [x for x in ev if x < F(0)], [x for x in ev if F(0) < x], F(0)
    selected = [("left", a, b) for a, b in zip(left[-4:] + [zero], left[-3:] + [zero])]
    selected += [("right", a, b) for a, b in zip([zero] + right[:4], right[:4])]
    assert len(selected) == 8
    C = decode(F, json.loads((root / "lower_round003/a2_equals_a3_results.json").read_text())["coefficients"])
    rows = []
    for i, (side, lo, hi) in enumerate(selected):
        obj = json.loads((root / f"lower_round007/exactC_event_interval_{i}_results.json").read_text())
        row = obj["row"]
        assert obj["status"] == "finalized_event_interval_exact_C"
        assert row["side"] == side
        assert decode(F, row["left_delta_coefficients"]) == lo
        assert decode(F, row["right_delta_coefficients"]) == hi
        assert not any(lo < z < hi for z in ev)
        xs = [Q(x) for x in row["sample_deltas"]]
        ys = [decode(F, item["ratio_coefficients"]) for item in row["sample_union_certificates"]]
        assert len(xs) == len(set(xs)) == len(ys) == 6
        assert all(lo < lift(F, x) < hi for x in xs)
        poly = interpolate(F, xs, ys)
        assert value(poly, xs[5]) == ys[5]
        assert [encode(x) for x in poly] == row["ratio_polynomial_coefficients"]
        d = list(poly); d[0] -= C
        b = bernstein4(d, lo, hi)
        assert all(x <= F(0) for x in b)
        if lo == zero or hi == zero: assert value(poly, Q(0)) == C
        rows.append({"index": i, "side": side, "no_event_in_open_interval": True,
                     "six_exact_certificate_values": True, "sixth_value_validates_polynomial": True,
                     "bernstein_P_minus_C_nonpositive": True,
                     "delta0_equality": (lo != zero and hi != zero) or value(poly, Q(0)) == C})
    hashes = {x: sha(root / x) for x in REQUIRED}
    result = {"status": "PASS_portable_lower007_event_audit", "raw_endpoint_form_count": len(forms),
              "exact_event_count": len(ev), "selected_intervals": rows, "input_sha256": hashes,
              "structural_degree_proof": "All raw affine endpoint crossings are excluded from each open interval. Hence endpoint ordering, clamp branch and slab coverage mask are fixed. The first coordinate is constant and exactly four tail factors are affine, so degree<=4; sample fitting validates but does not establish this fact.",
              "shared_arithmetic_limitation": "Uses archived lower/certify_rectangular.py for Q(u,v) arithmetic and isolators. Does not import the producer event/interpolation/Bernstein code and does not access separate sample JSON files or rerun unions.",
              "source_sha256": sha(Path(__file__))}
    raw = json.dumps(result, indent=2, sort_keys=True) + "\n"; args.output.write_text(raw); print(raw, end="")


if __name__ == "__main__": main()
