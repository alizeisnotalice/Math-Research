#!/usr/bin/env python3
"""Exact normalized-coordinate check of staggered block extrema/symmetries."""
from fractions import Fraction
from itertools import product
from pathlib import Path
import hashlib
import json


def sign(q):
    return 1 if q%2 == 0 else -1


def main():
    checked = 0
    for d in map(Fraction, ('1/32','1/16','1/8')):
        for qi,qj,li,lj in product(range(-3,4),range(-3,4),range(1,5),range(1,5)):
            atoms = [(ki,Fraction(kj)+d*sign(ki))
                     for ki in range(qi,qi+li) for kj in range(qj,qj+lj)]
            offsets = [d*sign(ki) for ki in range(qi,qi+li)]
            assert len(set(atoms)) == li*lj
            assert min(y for _,y in atoms) == qj+min(offsets)
            assert max(y for _,y in atoms) == qj+lj-1+max(offsets)
            reflected_i = {(-ki,y) for ki,y in atoms}
            q_ref = -qi-li+1
            expected_i = {(ki,Fraction(kj)+d*sign(ki))
                          for ki in range(q_ref,q_ref+li) for kj in range(qj,qj+lj)}
            assert reflected_i == expected_i
            reflected_j = {(ki,-y) for ki,y in atoms}
            expected_negative = {(ki,Fraction(kj)-d*sign(ki))
                                 for ki in range(qi,qi+li) for kj in range(-qj-lj+1,-qj+1)}
            assert reflected_j == expected_negative
            checked += 1
    here = Path(__file__).resolve()
    result = {'status':'all_exact_geometry_identity_checks_passed', 'block_cases':checked,
              'scope':'exact block extrema and source reflections; no volume/C6 certificate',
              'deltas':['1/32','1/16','1/8'], 'starts_each_axis':[-3,3],
              'lengths_each_axis':[1,4], 'code_sha256':hashlib.sha256(here.read_bytes()).hexdigest()}
    here.with_name('geometry_results.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__ == '__main__':
    main()
