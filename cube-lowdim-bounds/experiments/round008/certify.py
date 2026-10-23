#!/usr/bin/env python3
"""Rational-interval audit of the conservative round-008 two-bump kernel.

This program deliberately certifies only the one-dimensional lower bound
supplied by ``angular_reduction.md``.  It is not a sampled
grid: every accepted rational interval is bounded using convergent series and
outward rational enclosures for logarithms.
"""
from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction as Q
from functools import lru_cache
from hashlib import sha256
import json
from decimal import Decimal, localcontext
from pathlib import Path
from time import monotonic

ROOT = Path(__file__).resolve().parent
R = Q(4_746_319, 2_500_000)
E2 = Q(13, 20)
E3 = Q(13, 50)
PREVIOUS_UPPER = Q(27_226_955_153_507, 7_019_746_875_000)
MAX_BOXES = 100_000
SECONDS = 285.0


@dataclass(frozen=True)
class I:
    lo: Q
    hi: Q
    def __post_init__(self):
        if self.lo > self.hi:
            raise ValueError((self.lo, self.hi))
    @classmethod
    def point(cls, x: Q | int) -> "I":
        return cls(Q(x), Q(x))
    def __add__(self, other):
        other = interval(other)
        return I(self.lo + other.lo, self.hi + other.hi)
    __radd__ = __add__
    def __neg__(self): return I(-self.hi, -self.lo)
    def __sub__(self, other): return self + -interval(other)
    def __rsub__(self, other): return interval(other) + -self
    def __mul__(self, other):
        other = interval(other)
        vs = [self.lo * other.lo, self.lo * other.hi,
              self.hi * other.lo, self.hi * other.hi]
        return I(min(vs), max(vs))
    __rmul__ = __mul__
    def __truediv__(self, other):
        other = interval(other)
        if other.lo <= 0 <= other.hi:
            raise ZeroDivisionError(other)
        return self * I(Q(1) / other.hi, Q(1) / other.lo)
    def __rtruediv__(self, other): return interval(other) / self
    def square(self): return self * self
    def pair(self): return [str(self.lo), str(self.hi)]


def interval(x) -> I:
    return x if isinstance(x, I) else I.point(x)


@lru_cache(maxsize=None)
def _log_unit(x: Q, n: int = 72) -> I:
    """Outward rational enclosure of log(x), via the atanh series."""
    # Called only for 1 <= x <= 2, where the atanh series is fast.
    z = (x - 1) / (x + 1)
    partial = sum((z ** (2 * j + 1) / (2 * j + 1) for j in range(n)), Q(0)) * 2
    error = 2 * abs(z) ** (2 * n + 1) / ((2 * n + 1) * (1 - z * z))
    # Compress outward, keeping certificates and runtime manageable.
    D = 10 ** 32
    lo0, hi0 = partial - error, partial + error
    lo = Q(lo0.numerator * D // lo0.denominator, D)
    hi = Q(-((-hi0.numerator * D) // hi0.denominator), D)
    return I(lo, hi)


@lru_cache(maxsize=None)
def log_point(x: Q) -> I:
    """Outward log enclosure after exact binary range reduction to [1,2]."""
    x = Q(x)
    if x <= 0:
        raise ValueError(x)
    y, k = x, 0
    while y < 1:
        y *= 2; k -= 1
    while y > 2:
        y /= 2; k += 1
    if y == 2:
        y = Q(1); k += 1
    return _log_unit(y) + k * _log_unit(Q(2))


def log_i(x: I) -> I:
    return I(log_point(x.lo).lo, log_point(x.hi).hi)


def xlogx(x: I) -> I:
    """Interval for x log x on [0,1], including its continuous value at 0."""
    if x.lo < 0 or x.hi > 1:
        raise ValueError(x)
    vals_lo: list[Q] = []
    vals_hi: list[Q] = []
    if x.lo:
        a = interval(x.lo) * log_point(x.lo)
        vals_lo.append(a.lo); vals_hi.append(a.hi)
    else:
        vals_lo.append(Q(0)); vals_hi.append(Q(0))
    if x.hi:
        b = interval(x.hi) * log_point(x.hi)
        vals_lo.append(b.lo); vals_hi.append(b.hi)
    else:
        vals_lo.append(Q(0)); vals_hi.append(Q(0))
    # x log x reaches -1/e in (1/3,3/8).  -3/8 is a strict rational lower bound.
    if x.lo <= Q(3, 8) and x.hi >= Q(1, 3):
        vals_lo.append(Q(-3, 8))
    return I(min(vals_lo), max(vals_hi))


def tail_series(r: I, terms: int = 28) -> I:
    """Enclose f(r)=sum r^n/((n+2)R^(n+1)) on a nonnegative interval."""
    if r.lo < 0 or r.hi > 1:
        raise ValueError(r)
    def partial(t: Q) -> Q:
        return sum((t ** n / ((n + 2) * R ** (n + 1)) for n in range(terms + 1)), Q(0))
    low = partial(r.lo)
    high = partial(r.hi)
    # Bound n>=terms+1 by a geometric series after replacing n+2 by terms+3.
    rem = r.hi ** (terms + 1) / ((terms + 3) * R ** (terms + 2) * (1 - r.hi / R))
    return I(low, high + rem)


def b2_diag(r: I) -> I:
    """A psi_2 on the diagonal u=v=r/2."""
    s = I.point(1) - r
    return 8 * (2 - r - xlogx(s) + (2 - r) * log_i(r / 2))


def b3_diag(r: I) -> I:
    """A psi_3 on the diagonal u=v=r/2."""
    s = I.point(1) - r
    return 12 * (3 * (I.point(1) - r / 2).square()
                 + (s.square() + 1) * log_i(r / 2) - s * xlogx(s))


def g_diag(r: I) -> I:
    tail = 2 * (2 * tail_series(r) + I.point(1) / R) / (R - 1)
    return tail - E2 * b2_diag(r) - E3 * b3_diag(r)


def c_bound(r: I) -> I:
    # C(r)=3-4log2+4logr+1/(2r²), increasing on [3/4,1].
    return 3 - 4 * log_point(2) + 4 * log_i(r) + I.point(1) / (2 * r.square())


def k_bound(r: I) -> I:
    return (E2 * 64 * (I.point(1) - r).square() / (2 - r).square() / (2 - r)
            + I.point(1) / (R ** 3 * (R - 1)) - 3 * E3 * c_bound(r))


def lower_on_box(r: I) -> I:
    g = g_diag(r)
    if r.hi <= Q(3, 4):
        return g
    if r.lo < Q(3, 4):
        raise ValueError("box straddles 3/4")
    k = k_bound(r)
    # For d² in [0,r²], k d² is bounded below by min(k.lo,0)*r.hi².
    correction = min(k.lo, Q(0)) * r.hi * r.hi
    return g + correction


def small_r_certificate() -> dict:
    """Rational enclosure of the sign constant used for 0<r<=1/10.

    The angular-reduction proof supplies Apsi2<=0 and Apsi3<=0 if this is <0;
    the positive tail then gives G>=0 in that omitted endpoint strip.
    """
    h3 = 6 * (3 + Q(81, 100) * log_point(Q(1, 9)) + log_point(Q(1, 11)))
    return {"h3_upper": str(h3.hi), "h3_upper_display": float(h3.hi),
            "strictly_negative": h3.hi < 0}


def b3_slice_regression() -> list[dict]:
    """Numerical regression against the stated two-slice definition.

    This is a guard for the diagonal closed form, not part of the interval proof.
    """
    out = []
    for rr in (Q(1, 5), Q(3, 4), Q(99, 100)):
        with localcontext() as ctx:
            ctx.prec = 80
            r = Decimal(rr.numerator) / Decimal(rr.denominator)
            u = r / 2; a = 1 - u
            one_slice = 6 * (3 * a * a + ((a - u) ** 2 + (a + u) ** 2) * u.ln()
                             - (a - u) ** 2 * (a - u).ln()
                             - (a + u) ** 2 * (a + u).ln())
            direct = 2 * one_slice
            closed = b3_diag(I.point(rr))
            lo = Decimal(closed.lo.numerator) / Decimal(closed.lo.denominator)
            hi = Decimal(closed.hi.numerator) / Decimal(closed.hi.denominator)
            midpoint = (lo + hi) / 2
            out.append({"r": str(rr), "two_slice_decimal": str(direct),
                        "closed_interval_decimal": [str(lo), str(hi)],
                        "contained": lo <= direct <= hi,
                        "midpoint_difference": str(direct - midpoint)})
    return out


def audit_inner() -> dict:
    """Adaptive exact interval bisection; output is reproducible proof data."""
    start = monotonic()
    todo = [I(Q(1, 10), Q(3, 4)), I(Q(3, 4), Q(1))]
    accepted = 0
    splits = 0
    worst: I | None = None
    best_accepted: Q | None = None
    leaves: list[I] = []
    unresolved: list[I] = []
    while todo:
        if accepted + splits >= MAX_BOXES or monotonic() - start > SECONDS:
            unresolved.extend(todo)
            break
        box = todo.pop()
        lower = lower_on_box(box)
        worst = lower if worst is None or lower.lo < worst.lo else worst
        if lower.lo >= 0:
            accepted += 1
            leaves.append(box)
            best_accepted = lower.lo if best_accepted is None else min(best_accepted, lower.lo)
            if accepted % 250 == 0:
                print(f"progress accepted={accepted} splits={splits} pending={len(todo)}", flush=True)
            continue
        mid = (box.lo + box.hi) / 2
        # Rational bisection eventually resolves continuous positive margins;
        # a 2^-80 width floor documents a genuine enclosure limitation.
        if box.hi - box.lo <= Q(1, 2 ** 80):
            unresolved.append(box)
            continue
        todo.extend([I(box.lo, mid), I(mid, box.hi)])
        splits += 1
        if splits % 250 == 0:
            print(f"progress accepted={accepted} splits={splits} pending={len(todo)}", flush=True)
    elapsed = monotonic() - start
    leaves.sort(key=lambda b: b.lo)
    complete_cover = (not unresolved and leaves and leaves[0].lo == Q(1, 10)
                      and leaves[-1].hi == 1
                      and all(leaves[j].hi == leaves[j + 1].lo for j in range(len(leaves) - 1)))
    leaf_text = "\n".join(f"{b.lo}/{b.hi}" for b in leaves)
    return {
        "status": "certified" if not unresolved else "unresolved",
        "accepted_boxes": accepted, "splits": splits, "unresolved_boxes": len(unresolved),
        "elapsed_seconds": elapsed,
        "worst_raw_lower": None if worst is None else str(worst.lo),
        "worst_raw_lower_display": None if worst is None else float(worst.lo),
        "minimum_accepted_lower": None if best_accepted is None else str(best_accepted),
        "minimum_accepted_lower_display": None if best_accepted is None else float(best_accepted),
        "accepted_leaf_sha256": sha256(leaf_text.encode()).hexdigest(),
        "accepted_leaf_cover_verified": complete_cover,
        "unresolved_intervals": [b.pair() for b in unresolved[:20]],
    }


def outer_checks() -> dict:
    b0 = 2 - 2 * log_point(2)
    outer = I.point(1) - (E2 + E3) * R * (R - 1) * b0
    kernel = R - (E2 + E3) * (R - 1) / 4
    cost = R * R / (R - 1) - E2 / 6 - E3 / 10
    return {"outer_slack": outer.pair(), "outer_slack_display": [float(outer.lo), float(outer.hi)],
            "kernel_dominance_margin": str(kernel), "kernel_dominance_margin_display": float(kernel),
            "cost": str(cost), "cost_display": float(cost),
            "previous_upper": str(PREVIOUS_UPPER), "previous_upper_display": float(PREVIOUS_UPPER),
            "improvement": str(PREVIOUS_UPPER - cost), "improvement_display": float(PREVIOUS_UPPER - cost),
            "strict": outer.lo > 0 and kernel > 0 and cost < PREVIOUS_UPPER}


def main() -> None:
    inner = audit_inner()
    data = {
        "method": "exact Fraction intervals; logarithm atanh enclosures; adaptive 1D bisection",
        "parameters": {"R": str(R), "e2": str(E2), "e3": str(E3), "R_display": float(R)},
        "angular_reduction_dependency": "angular_reduction.md",
        "small_r": small_r_certificate(),
        "b3_diagonal_slice_regression": b3_slice_regression(),
        "inner": inner,
        "outer_and_cost": outer_checks(),
        "limits": {"max_boxes": MAX_BOXES, "seconds": SECONDS, "tail_terms": 28, "log_terms": 72},
    }
    data["scalar_status"] = ("certified_under_listed_analytic_reduction"
                             if inner["status"] == "certified" and inner["accepted_leaf_cover_verified"]
                             and data["small_r"]["strictly_negative"]
                             and data["outer_and_cost"]["strict"] else "not_closed")
    (ROOT / "certificate.json").write_text(json.dumps(data, indent=2, ensure_ascii=False) + "\n")
    print(json.dumps(data, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
