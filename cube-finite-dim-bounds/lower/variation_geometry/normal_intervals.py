"""Exact normal velocity coverage, with finite endpoints represented by Fraction.

left: velocities beta of old intervals ending at the face, hence (-inf,beta].
right: velocities alpha of old intervals starting there, hence [alpha,inf).
slabs: (alpha,beta) of newborn intervals. Reversed pairs are empty.
crossing: an old interval crosses the face with positive baseline margin.

This interface takes endpoint velocities, not outward velocities. In particular
the right endpoint in the n4 model is R=-minus, whereas L=plus.
"""
from fractions import Fraction
from itertools import product


def _q(value):
    if isinstance(value, float):
        raise TypeError('Pass exact integers/Fraction/decimal strings, not floats')
    return Fraction(value)


def merged_intervals(intervals):
    """Union components; ties merge, reversed pairs and points have zero length."""
    pairs = sorted((_q(a), _q(b)) for a,b in intervals)
    out = []
    for a,b in pairs:
        if a >= b:
            continue
        if out and a <= out[-1][1]:
            out[-1] = (out[-1][0], max(out[-1][1], b))
        else:
            out.append((a,b))
    return out


def signed_length(left=(), right=(), slabs=(), *, crossing=False):
    """Integral_R (1_velocity_union - 1_baseline_sides), as Fraction.

    A nonempty left/right list means that baseline half-line is present.
    Empty lists do not mean an endpoint at zero. With crossing=True the union
    is locally full and the answer is zero. No configuration counting is used.
    """
    left = tuple(map(_q,left))
    right = tuple(map(_q,right))
    finite = merged_intervals(slabs)
    if crossing:
        return Fraction(0)
    L = max(left) if left else None
    R = min(right) if right else None
    if L is not None and R is not None:
        if L >= R:
            return Fraction(0)
        filled = sum((max(Fraction(0), min(b,R)-max(a,L))
                      for a,b in finite), Fraction(0))
        return -(R-L)+filled
    if L is not None:
        return L + sum((max(Fraction(0),b-max(a,L))
                        for a,b in finite), Fraction(0))
    if R is not None:
        return -R + sum((max(Fraction(0),min(b,R)-a)
                         for a,b in finite), Fraction(0))
    return sum((b-a for a,b in finite), Fraction(0))


def pure_satellite_term(center, mass, old_boxes):
    """Separate N0=0 4D singleton term m*theta; not an O(epsilon) slab.

    old_boxes is a finite list of (lower_4tuple,upper_4tuple), in physical or
    fixed normalized coordinates consistently. Returns the physical volume
    coefficient if mass is the physical singleton mass coefficient. Old boxes
    must have O(epsilon) endpoint motion; distinct singleton centers must have
    a fixed positive separation, or co-located masses must first be aggregated.
    Caller must supply all relevant periodic images. See the proof addendum.
    """
    c = tuple(map(_q,center))
    if len(c) != 4:
        raise ValueError('Expected four coordinates')
    m = _q(mass)
    if m < 0:
        raise ValueError('Mass must be nonnegative')
    boxes = [(tuple(map(_q,lo)),tuple(map(_q,hi))) for lo,hi in old_boxes]
    if any(len(lo) != 4 or len(hi) != 4 for lo,hi in boxes):
        raise ValueError('Expected four-coordinate box endpoints')
    active = [(lo,hi) for lo,hi in boxes
              if all(lo[j] < hi[j] and lo[j] <= c[j] <= hi[j] for j in range(4))]
    uncovered = 0
    for signs in product((-1,1),repeat=4):
        covered = any(all((c[j] > lo[j] or signs[j] > 0)
                          and (c[j] < hi[j] or signs[j] < 0)
                          for j in range(4)) for lo,hi in active)
        uncovered += not covered
    return m * Fraction(uncovered,16)
