from fractions import Fraction
from itertools import product

# Capacities are measured in quarter units.
capacities = (2, 3, 4)
multiplicities = (1, 4, 2)
budget = 5
values = (sum(m*x for m, x in zip(multiplicities, allocation))
          for allocation in product(*(range(c+1) for c in capacities))
          if sum(allocation) <= budget)
assert Fraction(max(values), 4) == Fraction(3*4 + 2*2, 4)
