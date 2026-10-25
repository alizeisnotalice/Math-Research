#!/usr/bin/env python3
"""E's sole finite-step ideal union: rebuild 4096 boxes, no D dual-union call.

Reuses the previously audited 36-coefficient field arithmetic only.
The parameter is fixed by primary instructions, not read from D's answer.
"""
from argparse import ArgumentParser
from fractions import Fraction as Q
from functools import lru_cache
from itertools import product
from math import prod
from pathlib import Path
import hashlib
import importlib.util
import json
import os
import sys
import time

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
PROJECT = HERE.parents[2]
FACTOR = Q(1000001, 1000000)
TARGET = Q('2.06692413818')
ROOT_DENOMINATOR = 10**80
NODE_CAP, WALL_CAP = 10**6, 1200


def main():
    parser = ArgumentParser(description=__doc__)
    parser.add_argument('--arithmetic', type=Path,
                        default=PROJECT/'lower/round001_rectangular_certificate.py')
    parser.add_argument('--main-project', type=Path, default=PROJECT)
    parser.add_argument('--candidate-provenance', type=Path)
    args = parser.parse_args()
    started = time.monotonic()
    nodes = 0
    def budget():
        if time.monotonic()-started > WALL_CAP:
            raise TimeoutError('1200-second sole candidate cap')
        if nodes > NODE_CAP:
            raise RuntimeError('million union-node cap')
    paths = [args.arithmetic, args.main_project/'rounds/007/frozen_lower_before_candidate.json',
             args.main_project/'coordination/phase_20261025/frozen_thresholds.json']
    if args.candidate_provenance:
        paths.append(args.candidate_provenance)
    hashes = {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    ledger = json.loads(paths[1].read_text())
    frozen = json.loads(paths[2].read_text())
    current_upper = Q(ledger['lower']['exact_interval'][1])
    frozen_upper = Q(frozen['lower']['exact_interval'][1])
    # Provenance file is hashed only: no answer coefficients or derivative imported.
    spec = importlib.util.spec_from_file_location('audited_candidate_field', args.arithmetic)
    A = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(A)
    A.D = ROOT_DENOMINATOR
    A.RU = A.floor_root(2*A.D**6, 6)
    A.RV = A.floor_root(3*A.D**6, 6)
    A.UL, A.UH = Q(A.RU, A.D), Q(A.RU+1, A.D)
    A.VL, A.VH = Q(A.RV, A.D), Q(A.RV+1, A.D)
    assert A.UL**6 <= 2 < A.UH**6
    assert A.VL**6 <= 3 < A.VH**6
    A.bounds.cache_clear()
    F, u, v = A.F, A.u, A.v
    zero = F(0)
    assert u**6 == F(2) and v**6 == F(3)
    spacing = (1+u**4/2, FACTOR*(u**3+u**4*v**3)/4,
               *((1+u)/2 for _ in range(4)))
    assert all(zero < a for a in spacing)
    out = {'status': 'running', 'pid': os.getpid(),
           'scope': 'one fixed Cartesian candidate ideal union; finite L1 pending primary',
           'candidate_count': 1, 'candidate': {'a2_factor': str(FACTOR),
               'a2': '(1000001/1000000)*(u^3+u^4*v^3)/4',
               'other_spacings': '1+u^4/2; (1+u)/2 x4', 'unit_atom_mass': True},
           'input_sha256': hashes, 'root_denominator': str(A.D),
           'root_lower_numerators': [str(A.RU), str(A.RV)],
           'relations': ['u^6=2, u>0', 'v^6=3, v>0'],
           'limits': {'candidates': 1, 'union_nodes': NODE_CAP, 'wall_seconds': WALL_CAP},
           'current_complete_old_upper': str(current_upper),
           'frozen_complete_upper': str(frozen_upper), 'proposed_target': str(TARGET)}
    result_path = HERE/'candidate_results.json'
    def save():
        out['elapsed_seconds'] = time.monotonic()-started
        out['union_nodes'] = nodes
        out['code_sha256'] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
        result_path.write_text(json.dumps(out, indent=2)+'\n')
    save()
    try:
        @lru_cache(None)
        def radius(m):
            r = A.root(m)/2
            assert (2*r)**6 == F(m)
            return r
        records = []
        examined = 0
        endpoint_checks = 0
        for lengths in product(range(1, 5), repeat=6):
            budget()
            examined += 1
            mass = prod(lengths)
            r = radius(mass)
            starts = tuple(-(ell//2) if ell % 2 else 1-ell//2 for ell in lengths)
            box = []
            for a, q, ell in zip(spacing, starts, lengths):
                low = max(zero, (q+ell-1)*a-r)
                high = min(a/2, q*a+r)
                if high <= low:
                    break
                box.extend((low, high))
            else:
                for h, (a, q, ell) in enumerate(zip(spacing, starts, lengths)):
                    for center in box[2*h:2*h+2]:
                        for source_index in (q, q+ell-1):
                            assert -r <= center-source_index*a <= r
                            endpoint_checks += 1
                    assert zero <= box[2*h] < box[2*h+1] <= a/2
                records.append((lengths, tuple(box)))
        assert examined == 4096
        unique = {box: lengths for lengths, box in records}
        kept = []
        for box, lengths in unique.items():
            budget()
            if not any(other != box and all(other[2*h] <= box[2*h]
                        and box[2*h+1] <= other[2*h+1] for h in range(6))
                       for other in unique):
                kept.append((lengths, box))
        out.update({'length_vectors_examined': examined,
                    'nonempty_raw_boxes': len(records), 'unique_boxes': len(unique),
                    'maximal_boxes': len(kept), 'source_endpoint_checks': endpoint_checks})
        save()
        # Static field-coordinate slab recursion, independent of D's dual ranks.
        @lru_cache(None)
        def union(B):
            nonlocal nodes
            nodes += 1
            budget()
            if not B:
                return zero
            if len(B[0]) == 2:
                ordered = sorted(B)
                a, b = ordered[0]
                total = zero
                for c, d in ordered[1:]:
                    if b < c:
                        total, a, b = total+b-a, c, d
                    else:
                        b = max(b, d)
                return total+b-a
            total = zero
            cuts = sorted({x for box in B for x in box[:2]})
            for a, b in zip(cuts, cuts[1:]):
                child = tuple(sorted({box[2:] for box in B if box[0] <= a and b <= box[1]}))
                total = total+(b-a)*union(child)
            return total
        value = 64*union(tuple(sorted(box for _, box in kept)))
        lo, hi = A.bounds(value.a)
        assert lo <= hi
        scale = 10**14
        out.update({'status': 'ideal_union_exact_pass_pending_primary_L1_and_independent_audit',
            'normalization': 'ordinary lattice: 64 reflected half-cells, period mass 1',
            'coefficients': {f'{i},{j}': str(value.a[6*i+j])
                             for i in range(6) for j in range(6) if value.a[6*i+j]},
            'exact_interval': [str(lo), str(hi)],
            'decimal_enclosure': [str(Q(lo.numerator*scale//lo.denominator, scale)),
                                  str(Q(-((-hi.numerator*scale)//hi.denominator), scale))],
            'strictly_above_current_complete_old_upper': lo > current_upper,
            'strictly_above_frozen_complete_upper': lo > frozen_upper,
            'strictly_above_proposed_target': lo > TARGET,
            'lower_minus_current_complete_old_upper': str(lo-current_upper),
            'lower_minus_proposed_target': str(lo-TARGET),
            'witness_lengths': [list(lengths) for lengths, _ in kept],
            'new_C6_bound_status': 'pending primary finite-L1 and F independent review',
            'no_D_dualunion_import': True, 'no_MATLAB': True,
            'actual_process_status': 'finished'})
    except (TimeoutError, RuntimeError, ArithmeticError, AssertionError) as exc:
        out.update(status='stopped', reason=repr(exc), new_C6_bound_status='not promoted',
                   actual_process_status='finished')
    save()
    print(json.dumps({k: v for k, v in out.items()
                      if k not in ('coefficients', 'witness_lengths', 'input_sha256',
                                   'current_complete_old_upper', 'frozen_complete_upper',
                                   'lower_minus_current_complete_old_upper', 'exact_interval')}, indent=2))


if __name__ == '__main__':
    main()
