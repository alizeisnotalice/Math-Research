"""Small exact polynomial/rational-power checks; no production calculation."""
from fractions import Fraction as Q
from math import comb
from pathlib import Path
import json

integer_cases = []
for p in range(1, 9):
    # Coefficients of E ell^p = 1 + gamma sum_{j>=1} coeff_j u^j,
    # u=g/(2r), gamma=g/rho. Compare direct integration and closed formula.
    direct = {j: Q((-1)**j * comb(p, j), j + 1) for j in range(1, p + 1)}
    closed = {j: Q((-1)**j * comb(p + 1, j + 1), p + 1) for j in range(1, p + 1)}
    assert direct == closed
    integer_cases.append(p)

fractional_cases = []
q = Q(15, 16)
for p in (Q(3, 2), Q(5, 4)):
    a = q**p.denominator
    u = 1 - a
    gamma = Q(1, 5)
    # Direct change y=z^denominator in the tail integral.
    direct_tail = Q(p.denominator, p.numerator + p.denominator) * (1 - q**(p.numerator + p.denominator))
    closed_tail = (1 - q**(p.numerator + p.denominator)) / (p + 1)
    assert direct_tail == closed_tail
    moment_direct = 1 - gamma + gamma * direct_tail/u
    moment_closed = 1 - gamma + gamma * closed_tail/u
    assert moment_direct == moment_closed and 0 < moment_direct <= 1
    assert u <= Q(1, 4)
    fractional_cases.append({"p": str(p), "u": str(u), "moment": str(moment_direct)})

r, rho, g = Q(2), Q(7), Q(1, 2)
mean = 1 - g*g/(4*r*rho)
second = 1 - g*g/(2*r*rho) + g**3/(12*r*r*rho)
variance = second - mean*mean
assert variance == g**3/(12*r*r*rho) - g**4/(16*r*r*rho*rho)
assert variance >= g**3/(24*r*r*rho)
assert Q(3, 8)*Q(1, 24) == Q(1, 64)
assert 128**2 == 16384

result = {"status": "PASS", "scope": "small exact identity checks only",
          "integer_moments": integer_cases, "fractional_moments": fractional_cases,
          "variance_check": str(variance), "sympy_available": False}
Path(__file__).with_name("identity_checks.json").write_text(json.dumps(result, indent=2) + "\n")
print(json.dumps({"status": "PASS", "integer_cases": len(integer_cases), "fractional_cases": len(fractional_cases)}))
