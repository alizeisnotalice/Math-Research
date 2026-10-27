"""Small exact algebra diagnostics; not a substitute for the all-parameter proof."""
from fractions import Fraction as Q
from pathlib import Path
import json

rows = []
for n in (1, 2, 7):
    for t in (Q(2), Q(7)):
        left = 1 / t
        sigma = (left + 1) / 2
        beta = 1 + left
        core = 1 - (beta / 2) ** n
        D = lambda s: s**n - ((s + left) / 2) ** n
        Dprime = n * sigma ** (n - 1) - Q(n, 2) * ((sigma + left) / 2) ** (n - 1)
        ell = 1 + 1 / (2 * t) - 1 / (2 * t * sigma)
        p_bracket = 1 - Q(1, 2) * ell ** (n - 1)
        assert D(left) == 0 and D(Q(1)) == core
        assert 0 < D(sigma) < core < 1
        assert Dprime > 0 and p_bracket > 0
        rows.append({"n": n, "t": str(t), "sigma": str(sigma),
                     "D_prime": str(Dprime), "P_prime_bracket": str(p_bracket)})
result = {"status": "PASS", "scope": "six exact rational algebra checks only", "rows": rows}
Path(__file__).with_name("algebra_check.json").write_text(json.dumps(result, indent=2) + "\n")
print(json.dumps({"status": "PASS", "checked_cases": len(rows)}))
