"""Modified portable copy: exact conditional coefficient check; no proof of H."""
from fractions import Fraction
from math import comb
from pathlib import Path
import json
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--input", type=Path, default=Path(__file__).resolve().parent / "coefficients.json")
source = parser.parse_args().input
rows = json.loads(source.read_text())
assert len(rows) == 15
print("Conditional arithmetic only; H is not checked.")
print("n | L coefficient in v^(n+7) | b_n,7 / 2^(n+7)")
for row in rows:
    n = row["n"]
    baseline = sum((-1) ** (7 - (q - 1) * n) * comb(n + 6, 7 - (q - 1) * n)
                   for q in range(1, 2 + 7 // n))
    kn = Fraction(comb(n, 2) * n * n, 27)
    normalized = baseline - (n + Fraction(26, 9)) * kn
    coefficient = 2 ** (n + 7) * normalized
    assert baseline == int(row["ln_formula"]) == int(row["ln_series"])
    assert row["series_matches"] is True
    assert coefficient == Fraction(row["conditional_coefficient"])
    assert row["zero"] is (coefficient == 0)
    print(f"{n} | {baseline} | {normalized}")
print("All 15 exact comparisons against stored Mathematica JSON: True")
print("Zero dimensions in stored cases:", [r["n"] for r in rows if r["zero"]])
assert Fraction(rows[0]["conditional_coefficient"]) == Fraction(-16015360, 243)
print("n=2 regression: -16015360/243")
