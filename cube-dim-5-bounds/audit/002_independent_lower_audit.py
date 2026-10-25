#!/usr/bin/env python3
"""Independent exact audit for the selected n=5 rectangular witness.

This implementation does not import the producing script.  It rebuilds the
degree-25 field, receiver boxes, dominated-box removal, and recursive sweep
union.  It compares its lower interval against the archived exact R4 upper
endpoint, never against a displayed decimal.
"""
from fractions import Fraction
from functools import lru_cache
from itertools import product
from math import prod
from pathlib import Path
import hashlib
import json


DEN = 10**70


def floor_kth_root(value, degree):
    low, high = 0, 1 << ((value.bit_length() + degree - 1) // degree)
    while high - low > 1:
        mid = (low + high) // 2
        if mid**degree <= value:
            low = mid
        else:
            high = mid
    return low


u_floor = floor_kth_root(2 * DEN**5, 5)
v_floor = floor_kth_root(3 * DEN**5, 5)
u_interval = (Fraction(u_floor, DEN), Fraction(u_floor + 1, DEN))
v_interval = (Fraction(v_floor, DEN), Fraction(v_floor + 1, DEN))


class Alg:
    """Q[u,v]/(u^5-2,v^5-3), stored in the u^i v^j basis."""

    __slots__ = ("coeff",)

    def __init__(self, value=0, coeff=None):
        self.coeff = tuple(coeff) if coeff is not None else (Fraction(value),) + (Fraction(0),) * 24

    def __hash__(self):
        return hash(self.coeff)

    def __eq__(self, other):
        return self.coeff == as_alg(other).coeff

    def __add__(self, other):
        other = as_alg(other)
        return Alg(coeff=tuple(x + y for x, y in zip(self.coeff, other.coeff)))

    __radd__ = __add__

    def __neg__(self):
        return Alg(coeff=tuple(-x for x in self.coeff))

    def __sub__(self, other):
        return self + -as_alg(other)

    def __rsub__(self, other):
        return as_alg(other) + -self

    def __mul__(self, other):
        other = as_alg(other)
        answer = [Fraction(0)] * 25
        for index, left in enumerate(self.coeff):
            if not left:
                continue
            i, j = divmod(index, 5)
            for other_index, right in enumerate(other.coeff):
                if not right:
                    continue
                k, ell = divmod(other_index, 5)
                ii, jj = i + k, j + ell
                answer[(ii % 5) * 5 + (jj % 5)] += left * right * 2 ** (ii // 5) * 3 ** (jj // 5)
        return Alg(coeff=answer)

    __rmul__ = __mul__

    def __truediv__(self, rational):
        return Alg(coeff=tuple(c / Fraction(rational) for c in self.coeff))

    def __pow__(self, exponent):
        result = Alg(1)
        for _ in range(exponent):
            result *= self
        return result

    def interval(self):
        low = high = Fraction(0)
        for slot, coefficient in enumerate(self.coeff):
            i, j = divmod(slot, 5)
            a, b = u_interval[0] ** i * v_interval[0] ** j, u_interval[1] ** i * v_interval[1] ** j
            if coefficient >= 0:
                low += coefficient * a
                high += coefficient * b
            else:
                low += coefficient * b
                high += coefficient * a
        return low, high

    def __lt__(self, other):
        other = as_alg(other)
        if self == other:
            return False
        low, high = (other - self).interval()
        if low > 0:
            return True
        if high < 0:
            return False
        raise ArithmeticError("unseparated algebraic comparison")

    def __le__(self, other):
        return self == other or self < other


def as_alg(value):
    return value if isinstance(value, Alg) else Alg(value)


u = Alg(coeff=[Fraction(0)] * 5 + [Fraction(1)] + [Fraction(0)] * 19)
v = Alg(coeff=[Fraction(0), Fraction(1)] + [Fraction(0)] * 23)
assert u**5 == Alg(2) and v**5 == Alg(3)


def fifth_root_of_2_3_smooth(number):
    exponent_2 = exponent_3 = 0
    while number % 2 == 0:
        number //= 2
        exponent_2 += 1
    while number % 3 == 0:
        number //= 3
        exponent_3 += 1
    assert number == 1
    return u**exponent_2 * v**exponent_3


def selected_boxes():
    spacing = [1 + u**3 / 2, (2 + u**4 * v) / 4 + Fraction(1, 64)] + [(1 + u) / 2] * 3
    raw = []
    for counts in product(range(1, 5), repeat=5):
        radius = fifth_root_of_2_3_smooth(prod(counts)) / 2
        limits = []
        for length, step in zip(counts, spacing):
            if length % 2:
                left, right = Alg(0), min(step / 2, radius - (length - 1) * step / 2)
            else:
                left, right = max(Alg(0), length * step / 2 - radius), step / 2
            if right <= left:
                break
            limits.extend((left, right))
        else:
            raw.append((counts, tuple(limits)))
    unique = {box: counts for counts, box in raw}
    retained = []
    for box, counts in unique.items():
        contained = any(
            other != box and all(other[2 * j] <= box[2 * j] and box[2 * j + 1] <= other[2 * j + 1] for j in range(5))
            for other in unique
        )
        if not contained:
            retained.append((counts, box))
    return retained


@lru_cache(maxsize=None)
def union_measure(boxes):
    """Exact recursive coordinate sweep, using a fresh projection at each slab."""
    if not boxes:
        return Alg(0)
    dimension = len(boxes[0]) // 2
    if dimension == 1:
        total = Alg(0)
        current_left = current_right = None
        for left, right in sorted(boxes):
            if current_left is None:
                current_left, current_right = left, right
            elif current_right < left:
                total += current_right - current_left
                current_left, current_right = left, right
            elif current_right < right:
                current_right = right
        return total + current_right - current_left
    endpoints = sorted({endpoint for box in boxes for endpoint in box[:2]})
    total = Alg(0)
    for left, right in zip(endpoints, endpoints[1:]):
        projected = tuple(sorted({box[2:] for box in boxes if box[0] <= left and right <= box[1]}))
        total += (right - left) * union_measure(projected)
    return total


def coefficients(element):
    return {f"{i},{j}": str(value) for i, value in enumerate(element.coeff) for i, j in [divmod(i, 5)] if value}


def main():
    retained = selected_boxes()
    ratio = 32 * union_measure(tuple(sorted(box for _, box in retained)))
    lower, upper = ratio.interval()
    package_root = Path(__file__).resolve().parents[2] / "cube-finite-dim-bounds" / "baseline" / "original_certificates.json"
    r4 = json.loads(package_root.read_text(encoding="utf-8"))[1]
    r4_upper = Fraction(r4["lower_bound_interval"][1])
    assert lower > r4_upper
    output = {
        "status": "independent exact reconstruction passed",
        "model": "rectangular n=5 one-orbit lattice; delta_a1=0, delta_a2=1/64",
        "root_interval_denominator": str(DEN),
        "retained_witness_count": len(retained),
        "retained_witness_counts": [list(counts) for counts, _ in retained],
        "coefficient_basis": "u^i v^j, u^5=2, v^5=3",
        "ratio_coefficients": coefficients(ratio),
        "ratio_interval": [str(lower), str(upper)],
        "r4_exact_upper": str(r4_upper),
        "strictly_above_r4": True,
        "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "r4_source_sha256": hashlib.sha256(package_root.read_bytes()).hexdigest(),
    }
    Path(__file__).with_name("002_independent_lower_results.json").write_text(
        json.dumps(output, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(output, indent=2))


if __name__ == "__main__":
    main()
