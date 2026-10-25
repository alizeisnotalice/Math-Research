#!/usr/bin/env python3
"""Three preregistered exact checks. No radius or receiver-space sampling."""
from fractions import Fraction as F
from itertools import product, combinations
from pathlib import Path
import argparse
import hashlib
import json
import os
import time

START = time.monotonic()
BASE = Path(__file__).resolve().parent
parser = argparse.ArgumentParser()
parser.add_argument('--output-dir', type=Path, default=BASE)
args = parser.parse_args()
args.output_dir.mkdir(parents=True, exist_ok=True)


def encode(x):
    if isinstance(x, F):
        return str(x.numerator) + '/' + str(x.denominator)
    if isinstance(x, dict):
        return {k: encode(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [encode(v) for v in x]
    return x


def save(name, value):
    p = args.output_dir / name
    p.write_text(json.dumps(encode(value), indent=2, sort_keys=True) + '\n')


def checkpoint(stage):
    elapsed = time.monotonic() - START
    save('checkpoint.json', {'pid': os.getpid(), 'stage': stage,
                            'elapsed_seconds': elapsed,
                            'recoverable_command': 'python -B verify.py',
                            'time_limit_seconds': 120})
    if elapsed > 120:
        raise TimeoutError('Checkpoint deadline exceeded; checkpoint retained.')


def box(center, radius):
    return tuple((x-radius, x+radius) for x in center)


def contains(outer, inner):
    return all(a <= c and d <= b for (a, b), (c, d) in zip(outer, inner))


def disjoint(a, b):
    return any(v < s or t < u for (u, v), (s, t) in zip(a, b))


checkpoint('registered_start')
signs = tuple(product((-1, 1), repeat=5))
a = F(7, 8)
radii = [F(1) if e == (1,)*5 else F(2, 3) for e in signs]
centers = [tuple(a*r*s for s in e) for e, r in zip(signs, radii)]
qs = [box(c, r) for c, r in zip(centers, radii)]
separation_ratios = []
for i, k in combinations(range(32), 2):
    distance = max(abs(x-y) for x, y in zip(centers[i], centers[k]))
    assert distance > max(radii[i], radii[k])
    separation_ratios.append(distance / max(radii[i], radii[k]))
core = tuple((max(q[j][0] for q in qs), min(q[j][1] for q in qs))
             for j in range(5))
assert core == ((-F(1, 12), F(1, 12)),)*5
adjacencies = 0
for i, k in combinations(range(32), 2):
    differing = [j for j in range(5) if signs[i][j] != signs[k][j]]
    if len(differing) == 1:
        j = differing[0]
        gap = abs(centers[i][j]-centers[k][j])
        assert gap > max(radii[i], radii[k])
        assert radii[i]+radii[k]-gap < min(radii[i], radii[k])
        adjacencies += 1
assert adjacencies == 80
result_a = {'center_pairs_checked': 496,
            'orthant_edges_checked': adjacencies,
            'minimum_separation_ratio': min(separation_ratios),
            'core': core, 'maximum_radius': max(radii),
            'all_core_widths_less_than_max_radius': True}
checkpoint('case_A_complete')

# One new long packet family, checked through exact enclosing boxes.
J, q = 64, F(1, 8)
a, b = F(8191, 8192), F(11, 4)
eta, eps = F(1, 2**24), F(1, 2**12)
rs = [q**j for j in range(J+1)]
ps = [(b*r, F(0), F(0), F(0), F(0)) for r in rs]
cores = [box(p, (1-a)*r) for p, r in zip(ps, rs)]
unions = [box(p, (1+a)*r) for p, r in zip(ps, rs)]
expanded = [box(p, (2+a)*r) for p, r in zip(ps, rs)]
v = [(2*r)**5 for r in rs]
S = sum(v, F(0))
mu = [(1+eta)*x for x in v]
mu[-1] += eps*S
m = sum(mu, F(0))
assert m == (1+eta+eps)*S
assert 2*a > 1
assert (b-(1+a)) > q*(b+(1+a))
assert b < 2+a
assert b-(1-a) > q*(b+(2+a))
for i, k in combinations(range(J+1), 2):
    assert disjoint(unions[i], unions[k])
    assert contains(expanded[i], expanded[k])

incidences = []
for i in range(J+1):
    row = []
    for k in range(J+1):
        inside = contains(expanded[i], cores[k])
        assert inside == (k >= i)
        assert inside or disjoint(expanded[i], cores[k])
        row.append(int(inside))
    incidences.append(row)
assert all(mu_i > vi for mu_i, vi in zip(mu, v))
charge = sum((j+1)*mu_j for j, mu_j in enumerate(mu))
charge_direct = sum(sum(mu[k]*incidences[i][k] for k in range(J+1))
                    for i in range(J+1))
assert charge == charge_direct
D_avg = 32*(m-S)
D_mult = F(0)
union_volume = (1+a)**5*S
sum_cube_volume = 32*S
assert charge/m > F(128, 127)
assert D_avg/m < F(1, 100)
assert union_volume/m > F(3198, 100)
assert sum_cube_volume-union_volume >= 0
result_b = {'packets': J+1, 'cubes': 32*(J+1), 'q': q, 'a': a, 'b': b,
            'eta': eta, 'extra_mass_fraction_of_critical_sum': eps,
            'exact_packet_pair_checks': (J+1)*J//2,
            'expanded_core_incidence_checks': (J+1)**2,
            'mass': m, 'critical_mass_sum': S,
            'D_mult': D_mult, 'D_avg_over_mass': D_avg/m,
            'witness_union_over_mass': union_volume/m,
            'raw_expanded_charge_over_mass': charge/m,
            'receiver_defect': 'not small and not replaced by witness-union gap',
            'all_sources_absolute_continuous': True}
checkpoint('case_B_complete')

# Exact LP certificate: all 32 cube constraints in a packet are z_j <= v_j.
# Primal z_j=v_j. Dual selects one such cube constraint per packet, weight 1.
primal = list(v)
dual_cube_weights = [F(1)]*(J+1)
assert all(F(0) <= z <= mu_j and z <= vi
           for z, mu_j, vi in zip(primal, mu, v))
assert all(w >= 1 for w in dual_cube_weights)
primal_value = sum(primal, F(0))
dual_value = sum(w*vi for w, vi in zip(dual_cube_weights, v))
assert primal_value == dual_value == S
removed = m-S
retained_charge = sum((j+1)*z for j, z in enumerate(primal))
assert removed == D_avg/32
assert retained_charge/S < 1/(1-q**5)
assert F(127, 4)/(1-q**5) < 32
assert F(6, 5)**5 > 2
assert F(17, 10)**5 + 16 < F(127, 4)
result_c = {'lp_variables': J+1, 'original_cube_constraints': 32*(J+1),
            'primal_dual_equal': True, 'optimal_retained_mass': S,
            'removed_mass': removed, 'removed_mass_equals_D_avg_over_32': True,
            'retained_expanded_charge_over_retained_mass': retained_charge/S,
            'strict_charge_bound': 1/(1-q**5),
            'local_constant_times_charge_bound': F(127, 4)/(1-q**5),
            'is_global_C5_certificate': False}
checkpoint('case_C_complete')

results = {'status': 'all_registered_checks_passed',
           'method': 'Python standard-library Fraction; exact primal/dual LP certificate',
           'case_A': result_a, 'case_B': result_b, 'case_C': result_c,
           'local_constants': {'127_over_4': F(127, 4),
                               'crude_small_external_mass_bound': F(17, 10)**5+16},
           'global_upper_bound_proved': 32,
           'new_strict_global_bound_proved': False}
save('results.json', results)
checkpoint('complete')
receipt = {'pid': os.getpid(), 'elapsed_seconds': time.monotonic()-START,
           'exit_status': 0,
           'results_sha256': hashlib.sha256((args.output_dir/'results.json').read_bytes()).hexdigest(),
           'script_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
           'registration_sha256': hashlib.sha256((BASE/'registration.json').read_bytes()).hexdigest()}
save('run_receipt.json', receipt)
print(json.dumps(receipt, sort_keys=True))
print('Exact checks passed; no strict global C5 bound claimed.')
