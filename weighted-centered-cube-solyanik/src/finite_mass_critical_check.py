"""Exact rational check for the finite-mass intrinsic critical construction.

Example: m=3, exponents (1,2,3), and total defects (2/3,1/3,0).
For the probability measure 1_[0,1] dx, A1=[0,2/3], A2=(2/3,1],
A3=empty and F=[0,1/N], E_i=A_i^c union F.  Any sufficiently large
centred interval contains [0,1], so its product is the exact global product.
"""
from fractions import Fraction

d = (Fraction(2, 3), Fraction(1, 3), Fraction(0, 1))
a = (1, 2, 3)
critical = Fraction(1, 1)
for di, ai in zip(d, a):
    critical *= (1 - di) ** ai

assert critical == Fraction(4, 27)
assert sum(d) == 1 and all(0 <= di < 1 for di in d)

for N in (10, 100, 1000):
    f = Fraction(1, N)
    # The F term is genuinely new only in E_1.
    product = (1 - d[0] + f) ** a[0]
    for di, ai in zip(d[1:], a[1:]):
        product *= (1 - di) ** ai
    assert product > critical
    print(
        f"N={N}: P_a(1)={critical}; global witness product={product}; "
        f"excess={product-critical}; halo quotient=1/f={1/f}"
    )
