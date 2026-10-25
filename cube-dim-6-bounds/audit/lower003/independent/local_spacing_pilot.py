#!/usr/bin/env python3
"""Independent local audit of the 96-box R6 lower construction.

This uses mpmath only for a finite *numerical* pilot.  It independently
re-enumerates all 4^6 common-radius witnesses for twelve registered rational
spacing directions.  A forward-mode jet differentiates the union expression
inside the baseline's fixed combinatorial cell.  The output deliberately
distinguishes this from a rigorous certificate.
"""
from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from functools import lru_cache
from itertools import product
import hashlib
import json
from pathlib import Path
import platform
import sys

import mpmath as mp

mp.mp.dps = 70
N = 6
R6_CERTIFIED_UPPER = mp.mpf(206513594276569) / 10**14
HERE = Path(__file__).resolve().parent


def rkey(x):
    """Stable diagnostic key; this is not used for mathematical rounding."""
    return mp.nstr(x, 68)


def root_mass(mass: int) -> mp.mpf:
    return mp.root(mass, N) / 2


def make_records(spacing):
    """Return every nonempty box and its endpoint-choice signature."""
    records = []
    for count in product(range(1, 5), repeat=N):
        mass = 1
        for ell in count:
            mass *= ell
        radius = root_mass(mass)
        box, signature = [], []
        for a, ell in zip(spacing, count):
            if ell % 2:
                left, x, y = mp.mpf(0), a / 2, radius - (ell - 1) * a / 2
                right, tag = (x, "O-half") if x <= y else (y, "O-radius")
            else:
                x, y, right = mp.mpf(0), ell * a / 2 - radius, a / 2
                left, tag = (x, "E-zero") if x >= y else (y, "E-radius")
            if right <= left:
                break
            box.extend((left, right))
            signature.append(tag)
        else:
            records.append((count, tuple(box), tuple(signature)))
    return records


def contained(inner, outer):
    return all(outer[2*j] <= inner[2*j] and inner[2*j+1] <= outer[2*j+1]
               for j in range(N))


def maximal_records(records):
    # Numeric values arise from the same mpmath expressions; key preserves
    # duplicate removal independently of the source verifier's field type.
    unique = {}
    for count, box, signature in records:
        unique[tuple(rkey(x) for x in box)] = (count, box, signature)
    values = list(unique.values())
    kept = []
    for count, box, signature in values:
        if not any(other_box != box and contained(box, other_box)
                   for _, other_box, _ in values):
            kept.append((count, box, signature))
    return kept


def union_numeric(boxes, trace=None):
    @lru_cache(None)
    def rec(B):
        if not B:
            return mp.mpf(0)
        if len(B[0]) == 2:
            ordered = sorted(B)
            a, b = ordered[0]
            total = mp.mpf(0)
            merges = []
            for c, d in ordered[1:]:
                if b < c:
                    total += b-a
                    merges.append("gap")
                    a, b = c, d
                elif d > b:
                    b = d
                    merges.append("extend")
                else:
                    merges.append("covered")
            if trace is not None:
                trace.append(("leaf", len(B), tuple(merges)))
            return total + b-a
        cuts = sorted({x for box in B for x in box[:2]})
        total = mp.mpf(0)
        child_sizes = []
        for a, b in zip(cuts, cuts[1:]):
            child = tuple(sorted({box[2:] for box in B if box[0] <= a and b <= box[1]}))
            child_sizes.append(len(child))
            total += (b-a) * rec(child)
        if trace is not None:
            # This contains the complete finite branch pattern (not endpoint
            # values) selected by the recursion.  Equal hashes are useful
            # diagnostics, but not a standalone proof of a common open cell.
            trace.append(("node", len(B[0]) // 2, len(B), len(cuts), tuple(child_sizes)))
        return total
    return rec(tuple(sorted(boxes)))


@dataclass(frozen=True)
class Jet:
    value: mp.mpf
    derivative: tuple

    @classmethod
    def const(cls, value):
        return cls(mp.mpf(value), (mp.mpf(0),) * N)

    def __add__(self, other):
        other = asjet(other)
        return Jet(self.value + other.value,
                   tuple(x+y for x, y in zip(self.derivative, other.derivative)))
    __radd__ = __add__
    def __neg__(self):
        return Jet(-self.value, tuple(-x for x in self.derivative))
    def __sub__(self, other):
        return self + -asjet(other)
    def __rsub__(self, other):
        return asjet(other) + -self
    def __mul__(self, other):
        other = asjet(other)
        return Jet(self.value * other.value,
                   tuple(self.value*y + other.value*x
                         for x, y in zip(self.derivative, other.derivative)))
    __rmul__ = __mul__
    def __truediv__(self, scalar):
        scalar = mp.mpf(scalar)
        return Jet(self.value / scalar, tuple(x/scalar for x in self.derivative))
    def __lt__(self, other):
        return self.value < asjet(other).value
    def __le__(self, other):
        return self.value <= asjet(other).value


def asjet(x):
    return x if isinstance(x, Jet) else Jet.const(x)


def jet_box_for_count(count, baseline_signatures, spacing):
    mass = 1
    for ell in count:
        mass *= ell
    radius = Jet.const(root_mass(mass))
    box = []
    for a, ell, tag in zip(spacing, count, baseline_signatures):
        if tag == "O-half":
            left, right = Jet.const(0), a/2
        elif tag == "O-radius":
            left, right = Jet.const(0), radius - (ell-1)*a/2
        elif tag == "E-zero":
            left, right = Jet.const(0), a/2
        elif tag == "E-radius":
            left, right = ell*a/2-radius, a/2
        else:
            raise ValueError(tag)
        box.extend((left, right))
    return tuple(box)


def union_jet(boxes):
    @lru_cache(None)
    def rec(B):
        if not B:
            return Jet.const(0)
        if len(B[0]) == 2:
            ordered = sorted(B, key=lambda b: (b[0].value, b[1].value))
            a, b = ordered[0]
            total = Jet.const(0)
            for c, d in ordered[1:]:
                if b < c:
                    total, a, b = total + b-a, c, d
                elif b < d:
                    b = d
            return total + b-a
        cuts = sorted({x for box in B for x in box[:2]}, key=lambda x: x.value)
        total = Jet.const(0)
        for a, b in zip(cuts, cuts[1:]):
            child = tuple(sorted({box[2:] for box in B if box[0] <= a and b <= box[1]},
                                 key=lambda z: tuple(x.value for x in z)))
            total += (b-a) * rec(child)
        return total
    return rec(tuple(sorted(boxes, key=lambda z: tuple(x.value for x in z))))


def signature(records, kept):
    return {
        "nonempty_counts": [list(c) for c, _, _ in records],
        "retained_counts": [list(c) for c, _, _ in kept],
        "all_endpoint_branches": [list(s) for _, _, s in records],
    }


def evaluate(spacing):
    records = make_records(spacing)
    kept = maximal_records(records)
    trace = []
    ratio = 2**N * union_numeric(tuple(b for _, b, _ in kept), trace)
    trace_hash = hashlib.sha256(repr(tuple(trace)).encode()).hexdigest()
    return records, kept, ratio, trace_hash


def main():
    u, v = mp.root(2, N), mp.root(3, N)
    base = (1 + u**4/2, (2 + u**5*v)/4, *((1+u)/2 for _ in range(4)))
    records, kept, base_ratio, base_trace_hash = evaluate(base)
    base_signature = signature(records, kept)
    # Start the symbolic differential calculation only after re-enumerating
    # the baseline; its retained count is an acceptance condition.
    if (len(records), len(kept)) != (127, 96):
        raise RuntimeError("baseline count differs from registered 127 -> 96")
    jets = []
    for i, a in enumerate(base):
        d = [mp.mpf(0)] * N
        d[i] = 1
        jets.append(Jet(a, tuple(d)))
    jboxes = tuple(jet_box_for_count(c, s, jets) for c, _, s in kept)
    jratio = union_jet(jboxes) * (2**N)
    if abs(jratio.value - base_ratio) > mp.mpf("1e-55"):
        raise RuntimeError("jet and numeric union disagree at baseline")

    h = Fraction(1, 10000)
    perturbations = []
    for i in range(N):
        for sign in (-1, 1):
            multiplier = Fraction(1, 1) + sign*h
            trial = list(base)
            trial[i] *= mp.mpf(multiplier.numerator) / multiplier.denominator
            rr, kk, value, trace_hash = evaluate(tuple(trial))
            trial_signature = signature(rr, kk)
            same_box_signature = trial_signature == base_signature
            perturbations.append({
                "coordinate_1_based": i+1,
                "multiplier": f"{multiplier.numerator}/{multiplier.denominator}",
                "candidate_boxes": len(rr), "retained_boxes": len(kk),
                "same_box_signature": same_box_signature,
                "union_recursion_trace_sha256": trace_hash,
                "same_union_recursion_trace": trace_hash == base_trace_hash,
                "ratio": mp.nstr(value, 45),
                "delta_from_baseline": mp.nstr(value-base_ratio, 35),
                "delta_from_current_R6_certificate_upper": mp.nstr(value-R6_CERTIFIED_UPPER, 35),
                "first_order_prediction": mp.nstr(jratio.derivative[i]*(trial[i]-base[i]), 35),
            })

    report = {
        "status": "numerical_local_pilot_only; no new bound certified",
        "dimension": 6,
        "working_precision_decimal_digits": mp.mp.dps,
        "current_complete_R6_certificate_interval": [
            "206513594276568/10^14", "206513594276569/10^14"],
        "baseline_spacing": [mp.nstr(a, 50) for a in base],
        "baseline_candidate_boxes": len(records), "baseline_retained_boxes": len(kept),
        "baseline_union_recursion_trace_sha256": base_trace_hash,
        "baseline_ratio": mp.nstr(base_ratio, 50),
        "jet_ratio": mp.nstr(jratio.value, 50),
        "derivative_d_ratio_d_spacing": [mp.nstr(x, 45) for x in jratio.derivative],
        "elasticity": [mp.nstr(a*x/jratio.value, 40)
                       for a, x in zip(base, jratio.derivative)],
        "directions_registration": "twelve: a_i -> (1 +/- 1/10000) a_i",
        "perturbations": perturbations,
        "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "environment": {"python": sys.version, "mpmath": mp.__version__, "platform": platform.platform()},
        "scope_warning": "Fixed-cell derivative and numerical perturbations are not exact certificates; active-cell equality is tested only at these twelve finite samples.",
    }
    (HERE / "local_spacing_pilot_results.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
