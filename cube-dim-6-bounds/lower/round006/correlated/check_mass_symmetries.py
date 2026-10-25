#!/usr/bin/env python3
"""Small exact-arithmetic check of the derived correlated mass identities."""
from fractions import Fraction
from itertools import product
from pathlib import Path
import hashlib
import json


def ps(q, length):
    return 0 if length % 2 == 0 else (-1 if q % 2 else 1)


def main():
    checks = 0
    for d in map(Fraction, ('1/20', '1/10', '1/5')):
        assert sum(1 + d * (1 if (ki+kj)%2 == 0 else -1)
                   for ki,kj in product(range(2), repeat=2)) == 4
        for qi,qj,li,lj in product(range(-3,4),range(-3,4),range(1,5),range(1,5)):
            direct = sum(1 + d * (1 if (ki+kj)%2 == 0 else -1)
                         for ki in range(qi,qi+li) for kj in range(qj,qj+lj))
            formula = li*lj+d*ps(qi,li)*ps(qj,lj)
            assert direct == formula
            # Reflect all indices in axis i: q'=-q-li+1.
            assert ps(-qi-li+1,li) == ps(qi,li)
            # x_i -> a_i-x_i: q'=2-q-li and delta flips.
            assert li*lj-d*ps(2-qi-li,li)*ps(qj,lj) == formula
            checks += 1
    here = Path(__file__).resolve()
    result = {'status':'all_exact_identity_checks_passed', 'block_cases':checks,
              'deltas':['1/20','1/10','1/5'], 'starts_each_axis':[-3,3],
              'lengths_each_axis':[1,4], 'period_masses':[4,4,4],
              'scope':'exact identity checks; not a union volume or C6 certificate',
              'code_sha256':hashlib.sha256(here.read_bytes()).hexdigest()}
    here.with_name('mass_symmetries_results.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__ == '__main__':
    main()
