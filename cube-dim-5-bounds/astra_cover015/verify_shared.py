#!/usr/bin/env python3
"""Second preregistered batch: critical shared cliques, entirely exact."""
from fractions import Fraction as F
from itertools import product, combinations
from pathlib import Path
from math import prod
import argparse
import hashlib
import json
import os
import time

BASE = Path(__file__).resolve().parent
START = time.monotonic()
parser = argparse.ArgumentParser()
parser.add_argument('--output-dir', type=Path, default=BASE)
args = parser.parse_args()
args.output_dir.mkdir(parents=True, exist_ok=True)


def enc(x):
    if isinstance(x, F):
        return f'{x.numerator}/{x.denominator}'
    if isinstance(x, dict):
        return {k: enc(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [enc(v) for v in x]
    return x


def save(name, value):
    (args.output_dir/name).write_text(json.dumps(enc(value), indent=2, sort_keys=True)+'\n')


def checkpoint(stage):
    elapsed = time.monotonic()-START
    save('shared_checkpoint.json', {'pid': os.getpid(), 'stage': stage,
                                  'elapsed_seconds': elapsed, 'deadline_seconds': 120,
                                  'recoverable_command': 'python -B verify_shared.py'})
    if elapsed > 120:
        raise TimeoutError('Recoverable checkpoint retained')


checkpoint('registered_start')
D, A = 65536, 65535
den = 2*D
a, eta = F(A, D), F(1, 2**24)
center_coeffs = [c for c in product((-3, -1, 1, 3), repeat=5)
                 if any(abs(x) == 3 for x in c)]
centers = [(0,)*5] + [tuple(A*x for x in c) for c in center_coeffs]
radii = [2*D] + [D]*len(center_coeffs)
assert len(centers) == 993
pair_checks = 0
minimum_margin = None
for i, k in combinations(range(len(centers)), 2):
    dist = max(abs(x-y) for x, y in zip(centers[i], centers[k]))
    margin = dist-max(radii[i], radii[k])
    assert margin > 0
    minimum_margin = margin if minimum_margin is None else min(minimum_margin, margin)
    pair_checks += 1
    if pair_checks % 100000 == 0:
        checkpoint(f'center_pairs_{pair_checks}')
checkpoint('center_pairs_complete')

signs = tuple(product((-1, 1), repeat=5))
source_centers = [tuple(2*A*s for s in e) for e in signs]
source_halfwidth = D-A
source_masks = []
for p in source_centers:
    mask = 0
    for i, (c, r) in enumerate(zip(centers, radii)):
        inside = all(abs(x-y)+source_halfwidth <= r for x, y in zip(c, p))
        separated = any(abs(x-y)-source_halfwidth > r for x, y in zip(c, p))
        assert inside or separated
        if inside:
            mask |= 1 << i
    assert mask.bit_count() == 32 and mask & 1
    source_masks.append(mask)
assert len(set(source_masks)) == 32
for i in range(993):
    number_of_sources = sum(bool(mask & (1 << i)) for mask in source_masks)
    assert number_of_sources == (32 if i == 0 else 1)
checkpoint('source_incidence_complete')

# All coordinates have the same ten endpoints. Open-cell membership is exact.
endpoints = sorted(set(x+s*r for c, r in zip(centers, radii)
                       for x in [c[0]] for s in (-1, 1)))
intervals = [(u, v) for u, v in zip(endpoints, endpoints[1:]) if u < v]
coordinate_masks = []
for j in range(5):
    masks = []
    for lo, hi in intervals:
        mask = 0
        for i, (c, r) in enumerate(zip(centers, radii)):
            if c[j]-r <= lo and hi <= c[j]+r:
                mask |= 1 << i
        masks.append(mask)
    coordinate_masks.append(masks)
union_numerator = 0
maximum_multiplicity = 0
all_32_masks = set()
cells = 0
for ix in product(range(len(intervals)), repeat=5):
    mask = coordinate_masks[0][ix[0]]
    for j in range(1, 5):
        mask &= coordinate_masks[j][ix[j]]
    count = mask.bit_count()
    maximum_multiplicity = max(maximum_multiplicity, count)
    if count:
        union_numerator += prod(intervals[t][1]-intervals[t][0] for t in ix)
    if count == 32:
        all_32_masks.add(mask)
    cells += 1
assert maximum_multiplicity == 32
assert all_32_masks == set(source_masks)
union_volume = F(union_numerator, den**5)
assert union_volume == (1+3*a)**5
assert all(mask & 1 for mask in all_32_masks)
checkpoint('cell_arrangement_complete')

mass_per_core = 1+eta
m = 32*mass_per_core
sum_volumes = F(32+992)
D_avg = 32*m-sum_volumes
D_mult = F(0)
assert D_avg/m == 32*eta/(1+eta)
assert union_volume/m > F(3199, 100)
assert D_avg/m < F(1, 500000)

# Every packet dilation includes 2Q_central, hence contains every source core.
assert all(all(abs(x)+source_halfwidth < 4*D for x in p) for p in source_centers)
packet_costs = []
for selected in source_masks:
    remaining_multiplicities = [(mask & ~selected).bit_count() for mask in source_masks]
    assert sorted(remaining_multiplicities) == [0]+[31]*31
    charge = F(127, 4)*m + sum(remaining_multiplicities)*mass_per_core
    packet_costs.append(charge)
assert len(set(packet_costs)) == 1
assert packet_costs[0]/m == F(1977, 32)
no_packet_cost = 32*m
optimal_packet_cost = min([no_packet_cost]+packet_costs)
assert optimal_packet_cost == 32*m

# Uniform genuine receiver loss, using every positive radius analytically.
# Any good cube has r < (1+eta)^(1/5) < 101/100. An outer corner box
# can therefore see at most its matching source. That source's mass forces
# r < (1+eta)^(1/5)/2 < (1+eta)/2, but its distance is >= (2-a)/2.
assert F(101, 100)**5 > 1+eta
assert F(3, 2)+a-(1-a)/2 > F(101, 100)
assert eta < 1-a
assert (2-a)/2 > (1+eta)/2
excluded_receiver_volume = (3*a-2)**5
assert excluded_receiver_volume/m > F(1, 33)

results = {'status': 'all_shared_clique_checks_passed', 'a': a, 'eta': eta,
           'cube_count': 993, 'pair_checks': pair_checks,
           'minimum_center_separation_margin': F(minimum_margin, den),
           'arrangement_cells': cells, 'maximum_multiplicity': maximum_multiplicity,
           'number_of_32_cliques': len(all_32_masks),
           'every_32_clique_contains_central_cube': True,
           'all_sources_absolute_continuous': True,
           'mass': m, 'witness_union_volume': union_volume,
           'witness_union_over_mass': union_volume/m,
           'sum_cube_volumes': sum_volumes,
           'D_avg_over_mass': D_avg/m, 'D_mult': D_mult,
           'all_packet_dilations_contain_all_sources': True,
           'one_packet_charge_over_mass': packet_costs[0]/m,
           'best_disjoint_index_complete_packet_charge_over_mass': optimal_packet_cost/m,
           'source_capacity_critical_at_eta_zero': True,
           'excluded_outer_corner_receiver_volume': excluded_receiver_volume,
           'D_rec_over_mass_strict_lower_bound': F(1, 33),
           'receiver_defect': 'at least the excluded corner volume; in particular greater than mass/33',
           'proves_global_strict_upper_bound': False}
save('shared_results.json', results)
checkpoint('complete')
receipt = {'pid': os.getpid(), 'elapsed_seconds': time.monotonic()-START,
           'exit_status': 0,
           'results_sha256': hashlib.sha256((args.output_dir/'shared_results.json').read_bytes()).hexdigest(),
           'script_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
           'registration_sha256': hashlib.sha256((BASE/'registration_shared.json').read_bytes()).hexdigest()}
save('shared_run_receipt.json', receipt)
print(json.dumps(receipt, sort_keys=True))
print('Exact shared-clique obstruction verified; C5 improvement not claimed.')
