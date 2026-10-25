#!/usr/bin/env python3
"""Finite raw-endpoint screen around the R6-003 breakpoint construction."""
from __future__ import annotations

from fractions import Fraction
from functools import lru_cache
from itertools import product
import hashlib
import json
from pathlib import Path
import time

import mpmath as mp

mp.mp.dps = 70
N = 6
HERE = Path(__file__).resolve().parent
MAIN = Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube-dim-6/cube-dim-6-bounds')


def root_mass(m):
    return mp.root(m, N) / 2


def interval(a, ell, radius):
    if ell % 2:
        return mp.mpf(0), min(a/2, radius-(ell-1)*a/2)
    return max(mp.mpf(0), ell*a/2-radius), a/2


def records(spacing):
    out = []
    for count in product(range(1, 5), repeat=N):
        m = 1
        for e in count:
            m *= e
        radius = root_mass(m)
        box = []
        for a, ell in zip(spacing, count):
            lo, hi = interval(a, ell, radius)
            if hi <= lo:
                break
            box.extend((lo, hi))
        else:
            out.append((count, tuple(box)))
    return out


def contained(inner, outer):
    return all(outer[2*j] <= inner[2*j] and inner[2*j+1] <= outer[2*j+1]
               for j in range(N))


def kept_boxes(raw):
    unique = {}
    for count, box in raw:
        unique[tuple(mp.nstr(x, 68) for x in box)] = (count, box)
    values = list(unique.values())
    return [(c, b) for c, b in values
            if not any(other != b and contained(b, other) for _, other in values)]


def union(boxes):
    @lru_cache(None)
    def rec(B):
        if not B:
            return mp.mpf(0)
        if len(B[0]) == 2:
            a, b = sorted(B)[0]
            total = mp.mpf(0)
            for c, d in sorted(B)[1:]:
                if b < c:
                    total += b-a
                    a, b = c, d
                elif d > b:
                    b = d
            return total+b-a
        cuts = sorted({x for box in B for x in box[:2]})
        total = mp.mpf(0)
        for a, b in zip(cuts, cuts[1:]):
            child = tuple(sorted({box[2:] for box in B if box[0] <= a and b <= box[1]}))
            total += (b-a)*rec(child)
        return total
    return rec(tuple(sorted(boxes)))


def evaluate(spacing):
    raw = records(spacing)
    kept = kept_boxes(raw)
    return len(raw), len(kept), (2**N)*union(tuple(b for _, b in kept))


def base_spacing():
    u, v = mp.root(2, N), mp.root(3, N)
    return (1+u**4/2, (u**3+u**4*v**3)/4, *((1+u)/2 for _ in range(4)))


def fixed_other_nonempty(count, varying, spacing):
    m = 1
    for e in count:
        m *= e
    r = root_mass(m)
    for j, (a, ell) in enumerate(zip(spacing, count)):
        if j == varying:
            continue
        lo, hi = interval(a, ell, r)
        if hi <= lo:
            return False
    return True


def form(alpha, beta_kind, mass, ell, label):
    # beta_kind is 0, +r_m, or -r_m.  Keeping this tuple allows an exact
    # checker to reconstruct an algebraic event without trusting a decimal.
    return (Fraction(alpha), beta_kind, mass, ell, label)


def form_value(f):
    alpha, sign, m, _, _ = f
    beta = mp.mpf(0) if sign == 0 else sign*root_mass(m)
    return mp.mpf(alpha.numerator)/alpha.denominator, beta


def form_json(f):
    alpha, sign, m, ell, label = f
    return {"slope": f"{alpha.numerator}/{alpha.denominator}",
            "radius_sign": sign, "mass": m, "ell": ell, "role": label}


def endpoint_forms(varying, spacing):
    forms = set()
    for count in product(range(1, 5), repeat=N):
        if not fixed_other_nonempty(count, varying, spacing):
            continue
        m = 1
        for e in count:
            m *= e
        ell = count[varying]
        forms.add(form(0, 0, m, ell, "zero"))
        forms.add(form(Fraction(1, 2), 0, m, ell, "half"))
        if ell % 2:
            forms.add(form(-Fraction(ell-1, 2), +1, m, ell, "odd-radius"))
        else:
            forms.add(form(Fraction(ell, 2), -1, m, ell, "even-radius"))
    return sorted(forms)


def event_list(varying, spacing):
    forms = endpoint_forms(varying, spacing)
    events = []
    for i, f in enumerate(forms):
        af, bf = form_value(f)
        for g in forms[i+1:]:
            ag, bg = form_value(g)
            if af == ag:
                continue
            s = (bg-bf)/(af-ag)
            if s > 0:
                events.append((s, f, g))
    events.sort(key=lambda x: x[0])
    # Numerically merge only for selection; every retained description is
    # preserved so an exact audit can decide equality later.
    merged = []
    for s, f, g in events:
        if merged and abs(s-merged[-1][0]) < mp.mpf('1e-55')*max(1, abs(s)):
            merged[-1][1].append((f, g))
        else:
            merged.append([s, [(f, g)]])
    base = spacing[varying]
    # A symbolic base spacing can itself be an endpoint equality.  mpmath
    # evaluates identical algebraic expressions through different operation
    # trees, so exclude a symmetric high-precision neighbourhood rather than
    # accidentally replaying the already certified baseline as a "new" side.
    equal_tol = mp.mpf('1e-50') * max(1, abs(base))
    below = [x for x in merged if x[0] < base-equal_tol]
    above = [x for x in merged if x[0] > base+equal_tol]
    chosen = list(reversed(below[-8:])) + above[:8]
    return chosen, len(forms), len(merged)


def run_axis(varying, spacing, threshold):
    selected, form_count, event_count = event_list(varying, spacing)
    rows = []
    for s, pairs in selected:
        trial = list(spacing)
        trial[varying] = s
        raw, kept, value = evaluate(tuple(trial))
        rows.append({
            "event_value_decimal": mp.nstr(s, 50),
            "side": "below" if s < spacing[varying] else "above",
            "event_pairs_sample": [[form_json(a), form_json(b)] for a, b in pairs[:4]],
            "merged_pair_count": len(pairs),
            "raw_nonempty_boxes": raw,
            "maximal_boxes": kept,
            "ratio_numerical": mp.nstr(value, 50),
            "above_current_R6_003_upper_numerically": bool(value > threshold),
            "difference_from_current_upper": mp.nstr(value-threshold, 35),
        })
    return {"coordinate_1_based": varying+1, "base_spacing": mp.nstr(spacing[varying], 50),
            "endpoint_form_count": form_count, "positive_event_count_before_selection": event_count,
            "evaluated_events": rows}


def main():
    t0 = time.monotonic()
    ledger_path = MAIN/'results/bounds_ledger.json'
    ledger = json.loads(ledger_path.read_text())
    current = ledger['lower']
    assert current['id'] == 'R6-003-breakpoint'
    upper = Fraction(current['exact_interval'][1])
    threshold = mp.mpf(upper.numerator)/upper.denominator
    spacing = base_spacing()
    baseline = evaluate(spacing)
    out = {
        "status": "numerical_finite_endpoint_screen_only; no certification",
        "base_id": current['id'], "base_exact_upper": str(upper),
        "base_ratio_numerical": mp.nstr(baseline[2], 50),
        "base_raw_nonempty_boxes": baseline[0], "base_maximal_boxes": baseline[1],
        "axes": [run_axis(0, spacing, threshold), run_axis(2, spacing, threshold)],
        "selection_rule": "at most eight nearest lower and eight nearest upper distinct pre-branch endpoint events per axis",
        "scope_warning": "Numeric results at events do not certify a lower bound or a stable active cell; identical box counts are not interpreted as identical union combinatorics.",
        "ledger_sha256": hashlib.sha256(ledger_path.read_bytes()).hexdigest(),
        "elapsed_seconds": time.monotonic()-t0,
        "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    }
    (HERE/'endpoint_event_screen_results.json').write_text(json.dumps(out, indent=2)+'\n')
    print(json.dumps(out, indent=2))


if __name__ == '__main__':
    main()
