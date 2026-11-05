"""Finite exact cross-checks; all-n formulas are proved in proof_record.md."""
from pathlib import Path
from math import comb
from fractions import Fraction as F
import json

rows = []
for n in range(2, 13):
    k = F(comb(n, 2) * n * n, 27)
    base = sum(
        (-1) ** (7 - (j - 1) * n) * comb(n + 6, 7 - (j - 1) * n)
        for j in range(1, 2 + 7 // n)
    )
    rows.append({'n': n, 'normalized_b7': str(F(base) - (F(n) + F(26, 9)) * k)})
assert [r['normalized_b7'] for r in rows[:6]] == [
    '-31280/243', '676/9', '-21424/81', '-84575/243', '-8836/9', '-169442/81',
]
assert sum(
    (-1) ** (5 - (j - 1) * 3) * comb(7, 5 - (j - 1) * 3)
    for j in range(1, 2 + 5 // 3)
) == 0
exponents = []
for q in range(2, 25):
    z = F(q, 2 * q - 1)
    err = F(q) / (1 - z)
    assert err - 6 == F((q - 2) * (2 * q - 3), q - 1)
    exponents.append({'q': q, 'young_exponent': str(err)})
result = {
    'tool': 'Python3 fractions.Fraction and math.comb, actual new finite algebra run',
    'b7_rows': rows, 'b35_zero': True, 'young_exponents': exponents,
    'scope': 'finite exact checks only; all-n formula and geometry proved analytically in proof_record.md',
}
Path(__file__).with_suffix('.json').write_text(
    json.dumps(result, ensure_ascii=False, indent=2) + '\n'
)
print(json.dumps({
    'b7_first_six': rows[:6], 'b35_zero': True,
    'young_q2': exponents[0], 'young_q3': exponents[1], 'checks_passed': True,
}, ensure_ascii=False, indent=2))
