"""Exact rational checks for the two source geometries and Gaussian tail budget."""
from fractions import Fraction as F

source = F(2, 5) * 2 * 2
witness_area = 2
hit = F(2, 5) * (2 + 1)
joint = 2 * hit
independent = 2 * F(2, 5) * 2 * 2
assert hit / witness_area == F(3, 5)
assert joint == F(12, 5) and independent == F(16, 5)
K, h, q = 2, 16, 8
relaxed = h / (K + h + F(1, 4)) * (q - 1) / (q + F(1, 4))
assert relaxed == F(1792, 2409) and F(1, K + 1) < F(1, 2) < relaxed
assert F(1, 3) - F(384, 34**2) == F(1, 867)
print(dict(source=str(source), union=3, density=str(hit / witness_area),
           joint=str(joint), independent=str(independent), relaxed=str(relaxed),
           gaussian_margin_factor=str(F(1, 867))))
