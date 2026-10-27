"""Exact small diagnostics for the all-n outer-box proof; no optimization."""
from fractions import Fraction as Q
from math import comb
from pathlib import Path
import json

count = 0
guards = []
for n in (16, 64, 256):
    for q in range(1, n):
        assert Q(1, (2**q) * comb(n - 1, q)) <= Q(1, 2 * (n - 1))
        count += 1
    t = Q(n, 4)
    beta = 1 + 1/t
    r, s, R = Q(5, 4), Q(3, 2), Q(3)
    a = 10000 * t
    alpha = a/(a + 1)
    c = (a * beta - 1)/(a + 1)
    kappa = (1 + beta)/2
    for h in (Q(1), Q(6, 5)):
        assert 1 < r < s < t and h <= s
        assert r > beta * s - 1
        assert R > (1 + s)/beta and R < t
        assert R >= (alpha + s)/c
        assert kappa + 2/(a + 1) <= alpha * beta
        guards.append({"n": n, "t": str(t), "r": str(r), "s": str(s),
                       "h": str(h), "R": str(R), "kappa": str(kappa), "a": str(a)})
result = {"status": "PASS", "scope": "exact algebra diagnostics, not an all-n numerical proof",
          "integer_inequalities_checked": count, "finite_a_guard_cases": guards}
Path(__file__).with_name("box_algebra_check.json").write_text(json.dumps(result, indent=2) + "\n")
print(json.dumps({"status": "PASS", "integer_cases": count, "guard_cases": len(guards)}))
