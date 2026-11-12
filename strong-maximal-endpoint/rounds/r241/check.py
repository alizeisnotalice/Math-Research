from fractions import Fraction
from itertools import permutations
from math import factorial

for leaves in range(1, 7):
    # Label 0 is R, label 1 is the hub, and all other labels are leaves.
    kept = 0
    for order in permutations(range(leaves + 2)):
        earlier = order[:order.index(0)]
        kept += 1 not in earlier or len(earlier) == 1
    actual = Fraction(kept, factorial(leaves + 2))
    expected = Fraction(1, 2) + Fraction(1, (leaves + 1) * (leaves + 2))
    assert actual == expected
    print(leaves, actual)
