"""Exact finite algebra check for the round-11 Boolean-atom invariant.

This is evidence for the displayed rational algebra only; the general Radon
construction is proved in general_critical.tex.
"""

from fractions import Fraction


def colours(total, d, bias):
    return tuple(di * total + bi for di, bi in zip(d, bias))


def main():
    # Endpoint for exponents (1,1,2), with d=(1/2,1/2,0), F in colour 1.
    d = (Fraction(1, 2), Fraction(1, 2), Fraction(0))
    f = Fraction(7, 31)
    bias = (f / 2, -f / 2, Fraction(0))
    assert sum(bias) == 0

    total = f  # D_0=F
    assert colours(total, d, bias) == (f, Fraction(0), Fraction(0))
    for atom in (Fraction(5, 7), Fraction(11, 13), Fraction(17, 19)):
        total += atom  # exact proportional splitting of a new Boolean atom
        now = colours(total, d, bias)
        assert sum(now) == total
        assert tuple(ci - di * total for ci, di in zip(now, d)) == bias

    # A reference cube contains all old D and one new atom: same bias.
    q_mass = total + Fraction(23, 29)
    q_colours = colours(q_mass, d, bias)
    e_mass = tuple(q_mass - ci + (f if i == 0 else 0)
                   for i, ci in enumerate(q_colours))
    slack = tuple(ei - (1 - di) * q_mass for ei, di in zip(e_mass, d))
    assert slack == tuple(di * f for di in d)

    # If mu(Q triangle Q')<eps, the proof uses d_i f - 2 eps.
    eps = f / 8
    assert all(di * f - 2 * eps > 0 for di in d if di > 0)
    endpoint = Fraction(1, 4)
    assert (e_mass[0] / q_mass) * (e_mass[1] / q_mass) > endpoint
    print("PASS: exact bias invariant, reference slack, and error margin")


if __name__ == "__main__":
    main()
