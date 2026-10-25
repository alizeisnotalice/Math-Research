#!/usr/bin/env python3
"""Independent finite range check for the reflection-compressed period-2 interface."""
from fractions import Fraction
from itertools import product
import hashlib
import json
from pathlib import Path
import mpmath as mp

mp.mp.dps = 70
N = 6
HERE = Path(__file__).resolve().parent


def main():
    u, v = mp.root(2, 6), mp.root(3, 6)
    a1 = 1+u**4/2
    out = {"status": "numerical_range_check_only", "a1": mp.nstr(a1, 50), "deltas": {}}
    for delta_q in (Fraction(0), Fraction(1,20), Fraction(1,10), Fraction(1,5)):
        delta = mp.mpf(delta_q.numerator)/delta_q.denominator
        starts, original_indices, reflected_indices = set(), set(), set()
        for ell in product(range(1,5), repeat=N):
            rest = mp.mpf(1)
            for e in ell[1:]:
                rest *= e
            k = ell[0]
            # Deliberately use a larger finite integer range than the source
            # implementation's r_max pre-range, then impose the literal
            # interval-intersection condition.
            for s in range(-6, 7):
                selected = k if k % 2 == 0 else k + (delta if s % 2 == 0 else -delta)
                r = mp.root(rest*selected, 6)/2
                lo, hi = (s+k-1)*a1-r, s*a1+r
                if lo < a1 and hi > 0:
                    starts.add(s)
                    original_indices.update(range(s, s+k))
                    # x -> 2a1-x sends integer source index q to 2-q.
                    reflected_indices.update(2-q for q in range(s, s+k))
        union_indices = original_indices | reflected_indices
        out["deltas"][str(delta_q)] = {
            "first_axis_starts": sorted(starts),
            "original_atom_indices": [min(original_indices), max(original_indices)],
            "after_first_reflection_atom_indices": [min(union_indices), max(union_indices)],
            "coarse_offset_cells_after_reflection": [-1, 1],
        }
        assert starts == {-1,0,1}
        assert (min(original_indices), max(original_indices)) == (-1,2)
        assert (min(union_indices), max(union_indices)) == (-1,3)
    out["source_sha256"] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    (HERE/'start_range_check_results.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out, indent=2))


if __name__ == '__main__':
    main()
