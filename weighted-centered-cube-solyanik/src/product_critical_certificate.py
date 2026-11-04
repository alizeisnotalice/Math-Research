"""Exact finite certificates for the periodic-block algebra in round 10.

This program checks rational instances of the two identities used in the
proof.  It is a reproducible arithmetic check, not a substitute for the
measure-theoretic proof in paper/general_critical.tex (historical special proof: proofs/README.md).
"""

from fractions import Fraction


def floor_fraction(x: Fraction) -> int:
    return x.numerator // x.denominator


def stripe_prefix(x: Fraction, cuts: tuple[Fraction, ...], i: int) -> Fraction:
    """Mass of stripe i in [0,x] for the period-one stripe partition."""
    q = floor_fraction(x)
    r = x - q
    left = sum(cuts[:i])
    right = left + cuts[i]
    partial = max(Fraction(0), min(r, right) - left)
    return q * cuts[i] + partial


def stripe_mass(start: Fraction, length: Fraction,
                cuts: tuple[Fraction, ...], i: int) -> Fraction:
    return stripe_prefix(start + length, cuts, i) - stripe_prefix(start, cuts, i)


def check_periodic_identity(cuts: tuple[Fraction, ...]) -> None:
    starts = (Fraction(-17, 6), Fraction(-1, 7), Fraction(0),
              Fraction(5, 13), Fraction(31, 10))
    for start in starts:
        for integer_length in range(1, 10):
            for i, d_i in enumerate(cuts):
                observed = stripe_mass(start, Fraction(integer_length), cuts, i)
                assert observed == integer_length * d_i, (start, integer_length, i,
                                                           observed, integer_length * d_i)


def check_product_gap() -> None:
    # Equal exponents: the common F adds u to exactly one complement average.
    for m in (2, 3, 4, 8):
        c = Fraction(m - 1, m)
        for u in (Fraction(1, 101), Fraction(1, 17), Fraction(3, 10)):
            lhs = (c + u) * c ** (m - 1)
            critical = c ** m
            assert lhs - critical == u * c ** (m - 1) > 0

    # Unequal exponents (a1,a2,a3)=(1,1,2): its water-fill maximum at s=1
    # is d=(1/2,1/2,0), hence P=1/4.  F is placed in stripe 1.
    d = (Fraction(1, 2), Fraction(1, 2), Fraction(0))
    u = Fraction(1, 97)
    lhs = (1 - d[0] + u) * (1 - d[1]) * (1 - d[2]) ** 2
    critical = (1 - d[0]) * (1 - d[1]) * (1 - d[2]) ** 2
    assert critical == Fraction(1, 4)
    assert lhs - critical == u * Fraction(1, 2) > 0


def main() -> None:
    for cuts in ((Fraction(1, 2), Fraction(1, 2)),
                 (Fraction(1, 3), Fraction(1, 3), Fraction(1, 3)),
                 (Fraction(1, 2), Fraction(1, 2), Fraction(0)),
                 (Fraction(1, 5), Fraction(3, 10), Fraction(1, 2))):
        assert sum(cuts) == 1
        check_periodic_identity(cuts)
    check_product_gap()
    print("PASS: rational periodic identities and strict product gaps verified")


if __name__ == "__main__":
    main()
