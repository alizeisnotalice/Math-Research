"""Rational certificate for the weighted interval bound."""
from fractions import Fraction


def log_bounds(x):
    assert x >= 1
    z = (x - 1) / (x + 1)
    lower = sum(2 * z ** (2 * k + 1) / (2 * k + 1) for k in range(40))
    remainder = 2 * z**81 / (81 * (1 - z*z))
    return lower, lower + remainder


weight = Fraction(127, 100)
lower, upper = log_bounds(weight)
_, log_inverse_upper = log_bounds(1 / lower)
assert log_inverse_upper - 1 + (2 * weight + 1) * upper < Fraction(12774, 10000)
threshold = Fraction(1278, 1000)
assert log_bounds(1 / (threshold - 1))[0] > threshold
print("Certified c < 1.2774 < 1.278 < mu")
