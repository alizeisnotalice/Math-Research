#!/usr/bin/env python3
"""E's bounded dimer audit: 96 atom-interval cases, zero union evaluations.

No D module imports, MATLAB, old volume baseline or endpoint-volume rerun.
Only E's Fraction root/spacing and finite-loss helpers are reused.
"""
from fractions import Fraction as Q
from itertools import product
from math import prod
from pathlib import Path
import hashlib
import json
import sys
import time

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
from correlated_fraction_checker import root6, default_spacing, scale, ceilq, finite_transfer

MAIN = HERE.parents[2]
SOURCE = MAIN/'lower/round006/astra_construction'
INTERVAL_CAP, WALL_CAP = 96, 600


def active_position(k, epsilon):
    return Q(k)+epsilon*(k % 2)


def active_weight(k, delta):
    return 1+delta*(1-2*(k % 2))


def sufficient_collar(active_axis, starts, lengths):
    """For an arbitrary finite listed block family, no length assumption."""
    p = tuple(2 if h == active_axis else 1 for h in range(6))
    return tuple(ceilq(Q(max(abs(q), abs(q+ell-1)), ph))
                 for q, ell, ph in zip(starts, lengths, p))


def main():
    started = time.monotonic()
    inputs = [SOURCE/'report.md', SOURCE/'dimer_screen.py', SOURCE/'dimer_results.json',
              SOURCE/'registration.md', MAIN/'coordination/phase_20261025/frozen_thresholds.json',
              HERE.parent/'correlated_fraction_checker.py']
    hashes = {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs}
    screen_results = json.loads((SOURCE/'dimer_results.json').read_text())
    frozen = json.loads((MAIN/'coordination/phase_20261025/frozen_thresholds.json').read_text())
    assert screen_results['status'] == 'completed_no_promotion'
    assert len(screen_results['results']) == 24
    assert not any(x['promotion_trigger'] for x in screen_results['results'])
    assert Q(screen_results['frozen_complete_upper']) == Q(frozen['lower']['exact_interval'][1])

    # Global analytic constants, each checked in exact rational arithmetic.
    base = default_spacing(10**30)
    assert all(lo > 1 for lo, hi in base)
    # b_active=a_active/(1+epsilon)^(5/6)>a_active/2>1/2,
    # since (1+epsilon)^5<=2^5<2^6 for 0<=epsilon<=1.
    assert Q(2)**5 < Q(2)**6
    mass_cap = 8*4**5
    radius_cap = Q(9, 4)
    enumeration_radius_cap = Q(5, 2)
    assert (2*radius_cap)**6 > mass_cap
    assert (2*enumeration_radius_cap)**6 > 9*4**5

    cases, nonempty_outer, collar_endpoint_checks = 0, 0, 0
    normalized_scaling_checks = 0
    active_lengths_seen, inactive_lengths_seen = set(), set()
    arbitrary_collar_max = [0]*6
    sample_ranges = []
    for axis, epsilon, delta, ell in product(
            (0, 1, 2), (Q(1, 8), Q(7, 8)), (Q(0), Q(1, 5)), range(1, 9)):
        assert time.monotonic()-started < WALL_CAP
        assert cases < INTERVAL_CAP
        active_lengths_seen.add(ell)
        q = (-8, -5, -2, 0, 2, 6)[cases % 6]
        points = [active_position(k, epsilon) for k in range(q, q+ell)]
        weights = [active_weight(k, delta) for k in range(q, q+ell)]
        mass_active = Q(ell)+(delta*(1-2*(q % 2)) if ell % 2 else 0)
        span = Q(ell-1)+(epsilon*(1-2*(q % 2)) if ell % 2 == 0 else 0)
        assert min(points) == active_position(q, epsilon)
        assert max(points) == active_position(q+ell-1, epsilon)
        assert max(points)-min(points) == span
        assert sum(weights, Q(0)) == mass_active
        assert active_weight(0, delta)+active_weight(1, delta) == 2
        assert all(active_position(k+2, epsilon) == active_position(k, epsilon)+2
                   and active_weight(k+2, delta) == active_weight(k, delta)
                   for k in (q, q+ell-1))

        # Tail lengths exercise every permitted value without making a 6-D search.
        tails = tuple(1+(cases+h) % 4 for h in range(5))
        inactive_lengths_seen.update(tails)
        tail_mass = prod(tails)
        mass = mass_active*tail_mass
        indexed_atom_count = ell*tail_mass
        assert indexed_atom_count <= 8*4**5
        assert 0 < mass <= mass_cap
        r0, r1 = scale(root6(mass, 10**30), Q(1, 2))
        assert (2*r0)**6 <= mass < (2*r1)**6
        assert r1 < radius_cap

        # b scaling is enclosed independently; density identity is proved in report.
        s0, s1 = root6(1+epsilon, 10**30)
        assert s0**6/(1+epsilon) <= 1 < s1**6/(1+epsilon)
        ai = base[axis]
        bi = (ai[0]*s0/(1+epsilon), ai[1]*s1/(1+epsilon))
        assert bi[0] > Q(1, 2)
        normalized_scaling_checks += 1

        first = scale(bi, points[0])
        last = scale(bi, points[-1])
        outer_lo, outer_hi = max(Q(0), last[0]-r1), min(2*bi[1], first[1]+r1)
        if outer_lo < outer_hi:
            # This tests only interval feasibility; no length or union is summed.
            assert q >= -5 and q+ell-1 <= 6
            nonempty_outer += 1
        elif len(sample_ranges) < 3:
            sample_ranges.append({'start': q, 'length': ell, 'outer_empty': True})

        # Arbitrary finite-block collar, including sample blocks meeting no cell.
        starts, lengths = [0]*6, [0]*6
        starts[axis], lengths[axis] = q, ell
        htail = 0
        for h in range(6):
            if h == axis:
                continue
            length = tails[htail]
            starts[h] = -(length//2) if length % 2 else 1-length//2
            lengths[h] = length
            assert -2 <= starts[h] <= starts[h]+length-1 <= 2
            htail += 1
        collar = sufficient_collar(axis, starts, lengths)
        arbitrary_collar_max = [max(c, old) for c, old in zip(collar, arbitrary_collar_max)]
        periods = tuple(2 if h == axis else 1 for h in range(6))
        for h, (start, length, p, c) in enumerate(zip(starts, lengths, periods, collar)):
            signs = (1,) if h == axis else (-1, 1)
            for sign in signs:
                for endpoint in (start, start+length-1):
                    source_index = sign*endpoint
                    for receiver in (-3, 0, 3):
                        global_index = p*receiver+source_index
                        assert p*(-3-c) <= global_index <= p*(3+c+1)-1
                        collar_endpoint_checks += 1
        cases += 1
    assert cases == INTERVAL_CAP
    assert active_lengths_seen == set(range(1, 9))
    assert inactive_lengths_seen == set(range(1, 5))

    # Exact coarse mass with asymmetric periods; no physical-source enumeration.
    coarse_mass = {}
    k = 3
    for axis in (0, 1, 2):
        c = tuple(3 if h == axis else 2 for h in range(6))
        active_lo, active_hi = -2*(k+3), 2*(k+3+1)-1
        active_mass = sum((active_weight(n, Q(1, 5))
                           for n in range(active_lo, active_hi+1)), Q(0))
        exact = active_mass*prod(2*k+2*c[h]+1 for h in range(6) if h != axis)
        assert exact == 2*prod(2*k+2*x+1 for x in c)
        coarse_mass[str(axis)] = str(exact)
    # Conditional finite-loss lemma arithmetic, not a geometric C6 candidate.
    finite_example = finite_transfer(Q(1), Q(1), Q(1, 2),
                                    (3, 2, 2, 2, 2, 2), Q(1, 100), period_mass=2)
    assert finite_example['strictly_above_target']
    out = {'status': 'pass', 'new_C6_bound': False,
           'scope': 'dimerization model/normalization/range/collar/L1 interface only; zero union calls',
           'input_sha256': hashes,
           'source_numeric_batch_readonly': {'status': screen_results['status'], 'candidates': 24,
               'no_promoted_candidate': True, 'no_exact_candidate_volume_certification': True,
               'max_raw_boxes_readonly': max(x['raw_boxes'] for x in screen_results['results']),
               'frozen_complete_upper_exact_match': True},
           'limits': {'atom_interval_cases': INTERVAL_CAP, 'wall_seconds': WALL_CAP},
           'actual_atom_interval_cases': cases, 'atoms_per_interval_max': 8,
           'active_length_values': sorted(active_lengths_seen),
           'inactive_length_values': sorted(inactive_lengths_seen),
           'scaling_root_checks': normalized_scaling_checks,
           'outer_interval_possibly_meeting_cell': nonempty_outer,
           'collar_endpoint_comparisons': collar_endpoint_checks,
           'universal_bounds': {'mass_le': str(mass_cap), 'indexed_atoms_per_block_le': 8*4**5,
               'radius_lt': str(radius_cap),
               'active_spacing_gt': '1/2', 'retained_active_indices': [-5, 6],
               'retained_active_coarse_collar': 3, 'inactive_coarse_collar': 2,
               'enumerated_superset_active_indices': [-7, 8],
               'enumerated_superset_active_coarse_collar': 4},
           'exact_coarse_mass_checks_K3': coarse_mass,
           'finite_loss_conditional_example': {'assumed_raw_ratio': '1', 'comparison_target': '1/2',
               'is_C6_candidate': False, 'result': finite_example},
           'elapsed_seconds': time.monotonic()-started,
           'code_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
           'no_MATLAB': True, 'union_evaluations': 0}
    (HERE/'interface_audit_results.json').write_text(json.dumps(out, indent=2)+'\n')
    print(json.dumps(out, indent=2))


if __name__ == '__main__':
    main()
