#!/usr/bin/env python3
"""Exact motif identity checks plus a numeric merged-endpoint sanity check."""
from fractions import Fraction as Q
from pathlib import Path
import hashlib
import json
import os
import time
import dimer_screen as screen

HERE=Path(__file__).resolve().parent


def main():
    started=time.monotonic()
    tests=0
    b=Q(7,5)
    for e in (Q(0),Q(1,8),Q(1,3),Q(2,3),Q(7,8),Q(1)):
        for d in (Q(0),Q(1,5)):
            def p(k):
                return b*(k+e*(k%2))
            def w(k):
                return 1+d*(-1 if k%2 else 1)
            for q in range(-5,6):
                for ell in range(1,9):
                    points=[p(k) for k in range(q,q+ell)]
                    mass=sum(w(k) for k in range(q,q+ell))
                    expected_mass=ell+(d*(-1 if q%2 else 1) if ell%2 else 0)
                    expected_span=b*((ell-1)+(e*(1 if q%2==0 else -1) if ell%2==0 else 0))
                    assert mass==expected_mass
                    assert min(points)==p(q) and max(points)==p(q+ell-1)
                    assert max(points)-min(points)==expected_span
                    # Reflection about half the first gap exchanges parities.
                    reflected=[b*(1+e)-x for x in points]
                    transformed=[p(k) for k in range(2-q-ell,2-q)]
                    assert sorted(reflected)==transformed
                    assert sum(1-d*(-1 if k%2 else 1) for k in range(2-q-ell,2-q))==mass
                    # Reflection about zero changes only the signed displacement.
                    assert sorted([-x for x in points])==[
                        b*(k-e*(k%2)) for k in range(1-q-ell,1-q)]
                    tests+=1
    # Exact normalization follows from s^6=1+e: product(b_j)=product(a_j).
    # At e=1, a block of 2m starting at an odd q is m double atoms.
    endpoint_pairs=0
    for q in range(-5,6,2):
        for m in range(1,5):
            points=[k+(k%2) for k in range(q,q+2*m)]
            assert sorted(set(points))==list(range(q+1,q+2*m,2))
            assert all(points.count(x)==2 for x in set(points))
            endpoint_pairs+=1
    # This endpoint is prescribed by the theory, not a fresh search point.
    endpoint,_,_=screen.evaluate(2,Q(1),Q(1,5),started)
    baseline=json.loads((HERE/'dimer_results.json').read_text())
    endpoint['difference_from_complete_upper_float']=(
        endpoint['ratio_float']-float(Q(baseline['current_complete_upper'])))
    assert abs(endpoint['difference_from_complete_upper_float'])<2e-9
    out={'status':'PASS','pid':os.getpid(),'exact_Fraction_motif_cases':tests,
         'exact_double_atom_endpoint_cases':endpoint_pairs,
         'numeric_merged_endpoint':endpoint,
         'scope':'exact mass/span/reflection identities; endpoint volume float only, no C6 bound',
         'elapsed_seconds':time.monotonic()-started,
         'input_screen_sha256':hashlib.sha256((HERE/'dimer_screen.py').read_bytes()).hexdigest(),
         'input_result_sha256':hashlib.sha256((HERE/'dimer_results.json').read_bytes()).hexdigest(),
         'code_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
         'actual_process_status':'finished, no child process or MATLAB'}
    (HERE/'geometry_results.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out,indent=2))


if __name__=='__main__':
    main()
