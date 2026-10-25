#!/usr/bin/env python3
"""Small complexity audit for the 6D body-frame Fourier fan.

On the sorted positive chamber a0=1>=a1>=...>=a5>=0, this checks exactly
which of the 32 body forms can change sign at a full-dimensional point.  It
then uses finite LP feasibility only as a planning diagnostic for the number
of open sign cells; it is not an exact fan coverage certificate.
"""
from __future__ import annotations

import hashlib
import json
from fractions import Fraction
from itertools import product
from pathlib import Path

import numpy as np
from scipy.optimize import linprog

ROOT = Path(__file__).resolve().parent
FORMS = [(1,) + e for e in product((-1, 1), repeat=5)]
# The six vertices of the sorted chamber section a0=1.
VERTICES = [(1,) + (1,) * k + (0,) * (5 - k) for k in range(6)]


def form_value(form, point):
    return sum(Fraction(u) * v for u, v in zip(form, point, strict=True))


def crossing_forms():
    records = []
    for form in FORMS:
        values = [form_value(form, point) for point in VERTICES]
        records.append({"form": list(form), "minimum_vertex": str(min(values)), "maximum_vertex": str(max(values)), "strictly_crossing": min(values) < 0 < max(values)})
    return records


def max_margin(signs, forms):
    """Float LP only for a finite complexity estimate, never a proof claim."""
    # Variables are a1,...,a5,delta.  Chamber inequalities are a1<=1,
    # a2-a1<=0,...,a5-a4<=0,-a5<=0 and delta>=0.
    rows, rhs = [], []
    first = [0.0] * 6
    first[0] = 1.0
    rows.append(first); rhs.append(1.0)
    for j in range(1, 5):
        row = [0.0] * 6
        row[j] = 1.0; row[j - 1] = -1.0
        rows.append(row); rhs.append(0.0)
    row = [0.0] * 6; row[4] = -1.0
    rows.append(row); rhs.append(0.0)
    for sign, form in zip(signs, forms, strict=True):
        # f=1+sum form[j+1] a_j.  sign*f >= delta.
        row = [-sign * form[j + 1] for j in range(5)] + [1.0]
        rows.append(row); rhs.append(float(sign * form[0]))
    result = linprog([0, 0, 0, 0, 0, -1], A_ub=np.asarray(rows), b_ub=np.asarray(rhs), bounds=[(None, None)] * 5 + [(0, None)], method="highs")
    return float(result.x[-1]) if result.success else 0.0


def source_hash():
    return hashlib.sha256(Path(__file__).read_bytes()).hexdigest()


def main():
    records = crossing_forms()
    crossing = [tuple(record["form"]) for record in records if record["strictly_crossing"]]
    margins = []
    for signs in product((-1, 1), repeat=len(crossing)):
        margin = max_margin(signs, crossing)
        if margin > 1e-9:
            margins.append(margin)
    out = {
        "scope": "exact crossing-form count plus floating LP planning diagnostic, not a fan certificate",
        "chamber": "a0=1>=a1>=a2>=a3>=a4>=a5>=0",
        "form_records": records,
        "forms_total": len(FORMS),
        "strictly_crossing_forms_exact": len(crossing),
        "strictly_crossing_forms": [list(form) for form in crossing],
        "open_sign_cells_lp_diagnostic": len(margins),
        "smallest_positive_margin_lp": min(margins),
        "full_fan_estimate_if_each_cell_triangulates_once": len(margins) * (2**6) * math_factorial_6(),
        "next_algorithm": "replace LP diagnostic with exact rational polyhedral cell enumeration, then triangulate each cell and prove a.e. coverage before using Fourier sums",
        "source_sha256": source_hash(),
    }
    (ROOT / "body32_fan_complexity_results.json").write_text(json.dumps(out, indent=2) + "\n")
    print(json.dumps({k: out[k] for k in out if k not in {"form_records", "strictly_crossing_forms"}}, indent=2))


def math_factorial_6():
    return 720


if __name__ == "__main__":
    main()
