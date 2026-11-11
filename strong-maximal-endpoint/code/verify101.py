#!/usr/bin/env python3
"""Validate the archived finite-grid private-cell search evidence exactly."""
from __future__ import annotations

import hashlib
import json
import math
import platform
import sys
from itertools import product
from pathlib import Path
from fractions import Fraction

PACKAGE = Path(__file__).resolve().parents[1]
ROUND = PACKAGE / 'rounds' / 'r101'
DATA = ROUND / 'data'


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def catalog(shape):
    points = list(product(*(range(m) for m in shape)))
    axes = [
        {tuple(sorted((o + i) % m for i in range(s)))
         for o in range(m) for s in range(1, m + 1)}
        for m in shape
    ]
    boxes = set()
    for arcs in product(*axes):
        mask = 0
        for idx, point in enumerate(points):
            if all(x in arc for x, arc in zip(point, arcs)):
                mask |= 1 << idx
        boxes.add(mask)
    return len(points), boxes


def t_value(private_mask: int, boxes: set[int]) -> int:
    """T(P)=sum_i max |R| over catalog boxes with R∩P={i}."""
    score = 0
    todo = private_mask
    while todo:
        bit = todo & -todo
        todo ^= bit
        options = [b.bit_count() for b in boxes if b & private_mask == bit]
        if not options:
            raise AssertionError(f'private cell {bit} has no singleton-intersection box')
        score += max(options)
    return score


def validate_case(path: Path):
    row = json.loads(path.read_text())
    shape = tuple(map(int, row['shape']))
    n, boxes = catalog(shape)
    full = (1 << n) - 1
    total_masks = 1 << n
    if n > 32:
        raise AssertionError(f'case exceeds compiled solver mask capacity: {shape}')
    if row['cell_count'] != n or row['number_of_masks'] != total_masks:
        raise AssertionError(f'grid size/mask count mismatch: {shape}')
    if row['rectangles'] != len(boxes):
        raise AssertionError(f'catalog rectangle count mismatch: {shape}')
    if row['score_denominator'] != n or row['workers'] != 12:
        raise AssertionError(f'normalization or worker count mismatch: {shape}')

    chunks = row['chunks']
    if len(chunks) != 12:
        raise AssertionError(f'expected twelve chunks: {shape}')
    expected_ranges = [(total_masks * i // 12, total_masks * (i + 1) // 12)
                       for i in range(12)]
    previous = 0
    chunk_best = []
    chunk_arg_count = 0
    for idx, (chunk, expected) in enumerate(zip(chunks, expected_ranges)):
        if len(chunk) != 6:
            raise AssertionError(f'malformed chunk record {idx}: {shape}')
        lo, hi, best, arg, ties, checksum = map(int, chunk)
        if (lo, hi) != expected or lo != previous or hi <= lo:
            raise AssertionError(f'chunk partition gap/overlap or wrong endpoint: {shape}, chunk {idx}')
        if not (lo <= arg < hi and ties >= 1 and ties <= hi - lo):
            raise AssertionError(f'chunk arg/tie count outside range: {shape}, chunk {idx}')
        if checksum < 0 or checksum > best * (hi - lo):
            raise AssertionError(f'chunk score checksum outside its maximum bound: {shape}, chunk {idx}')
        if t_value(arg, boxes) != best:
            raise AssertionError(f'chunk maximum witness does not attain its reported score: {shape}, chunk {idx}')
        previous = hi
        chunk_best.append(best)
        chunk_arg_count += ties
    if previous != total_masks:
        raise AssertionError(f'chunk coverage ends before full mask range: {shape}')

    private = int(row['private_mask'])
    if not (0 < private <= full):
        raise AssertionError(f'invalid winning private-cell mask: {shape}')
    score = int(row['score_numerator'])
    if max(chunk_best) != score or t_value(private, boxes) != score:
        raise AssertionError(f'global maximum score/witness mismatch: {shape}')
    if private not in [int(c[3]) for c in chunks if int(c[2]) == score]:
        raise AssertionError(f'winning mask not reported by a maximum chunk: {shape}')

    max_boxes = list(map(int, row['max_boxes']))
    private_bits = [1 << i for i in range(n) if private >> i & 1]
    if len(max_boxes) != len(private_bits) or len(set(max_boxes)) != len(max_boxes):
        raise AssertionError(f'max-box list does not have one distinct box per private cell: {shape}')
    assigned = []
    for b in max_boxes:
        if b not in boxes:
            raise AssertionError(f'max box absent from catalog: {shape}')
        hit = b & private
        if hit.bit_count() != 1:
            raise AssertionError(f'max box does not have exactly one private cell: {shape}')
        if b.bit_count() != max(x.bit_count() for x in boxes if x & private == hit):
            raise AssertionError(f'max box is not volume-maximal for its private cell: {shape}')
        assigned.append(hit)
    if len(set(assigned)) != len(private_bits) or set(assigned) != set(private_bits):
        raise AssertionError(f'max boxes fail to biject with the private cells: {shape}')
    max_union = 0
    for b in max_boxes:
        max_union |= b
    if max_union != full:
        raise AssertionError(f'max-box union is not the full torus: {shape}')
    if sum(b.bit_count() for b in max_boxes) != score:
        raise AssertionError(f'max-box volumes do not sum to T(P): {shape}')

    max_box_cells = max(b.bit_count() for b in max_boxes)
    eta = Fraction(score, n)
    return {
        'shape': list(shape),
        'cell_count': n,
        'catalog_rectangle_count': len(boxes),
        'mask_range': f'[0,{total_masks})',
        'chunk_count': len(chunks),
        'chunks_partition_full_mask_range': True,
        'chunk_maxima_and_arg_witnesses_replayed': True,
        'chunk_score_sum': sum(chunk_best),
        'winning_private_mask': private,
        'winning_private_cell_count': len(private_bits),
        'T_private_mask_numerator': score,
        'eta_star_fraction': f'{eta.numerator}/{eta.denominator}',
        'eta_star_decimal': float(eta),
        'max_box_count': len(max_boxes),
        'max_box_cell_volumes': [b.bit_count() for b in max_boxes],
        'max_box_volume_fraction': f'{max_box_cells}/{n}',
        'max_box_union_cell_count': max_union.bit_count(),
        'max_box_union_fraction': f'{max_union.bit_count()}/{n}',
        'max_boxes_union_is_full_torus': True,
        'maximum_boxes_have_distinct_private_cells': True,
        'reported_seconds': float(row['seconds']),
    }


def validate_resource(path: Path):
    data = json.loads(path.read_text())
    if data.get('exit_code') != 0 or data.get('stop_reason') is not None:
        raise AssertionError(f'resource guard did not pass: {path.name}')
    if data.get('reserve_GiB', 0) < 1 or data.get('owned_RSS_budget_GiB', 0) > 6.5:
        raise AssertionError(f'resource policy mismatch: {path.name}')
    if data.get('peak_owned_RSS_MiB', math.inf) > 6.5 * 1024:
        raise AssertionError(f'peak RSS budget exceeded: {path.name}')
    if data.get('min_system_available_MiB', -math.inf) < 1024:
        raise AssertionError(f'memory reserve estimate not met: {path.name}')
    return data


def main():
    files = sorted(DATA.glob('grid_*.json'))
    if len(files) != 14:
        raise AssertionError(f'expected 14 grid result files, found {len(files)}')
    cases = [validate_case(path) for path in files]
    resource_names = ('resources101.json', 'extend101-resources.json')
    resources = [validate_resource(ROUND / name) for name in resource_names]
    max_case = max(cases, key=lambda c: Fraction(c['eta_star_fraction']))
    binary5 = next(c for c in cases if c['shape'] == [2, 2, 2, 2, 2])
    if binary5['eta_star_fraction'] != '45/16' or binary5['T_private_mask_numerator'] != 90:
        raise AssertionError('5D binary-grid T(P*) witness mismatch')
    eta_all = Fraction(max_case['eta_star_fraction'])
    if eta_all != Fraction(13, 4) or max_case['shape'] != [4, 8]:
        raise AssertionError('overall finite-grid eta* mismatch')
    validation = {
        'pass_all': True,
        'goal_status': 'active',
        'general_continuous_upper_claimed': False,
        'finite_catalog_only': True,
        'case_count': len(cases),
        'all_12_chunks_per_case_partition_full_range': True,
        'all_chunk_argmax_witnesses_replayed': True,
        'all_max_box_volumes_and_unions_replayed': True,
        'max_eta_star_fraction': f'{eta_all.numerator}/{eta_all.denominator}',
        'max_eta_star_shape': max_case['shape'],
        'five_dimensional_binary_eta_star_fraction': binary5['eta_star_fraction'],
        'positivity_condition': 'theta > (1-p)*eta_star',
        'uniform_all_grid_full_coverage_strict_p_interval': '(9/13,1)',
        'cases': cases,
        'resources': [{
            'file': name,
            'exit_code': r['exit_code'],
            'stop_reason': r['stop_reason'],
            'peak_owned_RSS_MiB': r['peak_owned_RSS_MiB'],
            'min_system_available_MiB': r['min_system_available_MiB'],
            'reserve_GiB': r['reserve_GiB'],
            'seconds': r['seconds'],
            'scope': r['scope'],
        } for name, r in zip(resource_names, resources)],
        'python_version': platform.python_version(),
        'source_sha256': {
            'search101.py': sha256(PACKAGE / 'code' / 'search101.py'),
            'private101.cpp': sha256(PACKAGE / 'code' / 'private101.cpp'),
            'verify101.py': sha256(Path(__file__)),
        },
        'input_sha256': {
            f'data/{p.name}': sha256(p) for p in files
        } | {
            name: sha256(ROUND / name) for name in resource_names
        },
        'document_sha256': ({
            'report.md': sha256(ROUND / 'report.md')
        } if (ROUND / 'report.md').exists() else {}),
        'prior_dependency_sha256': {},
        'scope': 'Exact validation of finite periodic-grid exhaustive-search evidence, private-cell maximum-box witnesses, and resource logs. The finite catalogs are not a proof of any continuous general upper bound.',
    }
    out = ROUND / 'validation101.json'
    out.write_text(json.dumps(validation, indent=2, sort_keys=True) + '\n')
    print(json.dumps({
        'pass_all': True,
        'case_count': len(cases),
        'max_eta_star': validation['max_eta_star_fraction'],
        'max_eta_star_shape': validation['max_eta_star_shape'],
        'five_dimensional_binary_eta_star': validation['five_dimensional_binary_eta_star_fraction'],
        'validation_path': str(out),
    }, indent=2))


if __name__ == '__main__':
    main()
