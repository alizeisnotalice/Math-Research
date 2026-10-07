#!/usr/bin/env python3
"""Exact postprocessing only: fixed single-gridpoint packet share bands."""
import gzip
import hashlib
import json
import math
from fractions import Fraction as F
from pathlib import Path
from time import perf_counter

BASE = Path(__file__).resolve().parent
PREFIX = 'packet_share_band_critical_audit_20261007'


def digest(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def band(p):
    assert 0 < p <= 1
    j, theta = 0, F(1, 2)
    while p <= theta:
        j += 1
        theta /= 2
    assert theta < p <= 2 * theta
    return j, theta


def main():
    start = perf_counter()
    regpath = BASE / (PREFIX + '_registration.json')
    reg = json.loads(regpath.read_text())
    for name, expected in reg['inputs_sha256'].items():
        assert digest(BASE / name) == expected, name
    prior = json.loads((BASE / 'packet_capture_critical_audit_20261007_results.json').read_text())
    rows = [r for r in prior['records'] if r['eta'] == '1/2']
    assert len(rows) == 24
    result = []
    rounds = []
    for n in (4, 16, 64):
        saved = json.loads(gzip.decompress((BASE / f'critical_source_tail_arrival_20261007_n{n}_exact_profiles.json.gz').read_bytes()))
        inp = saved['input']
        weights = list(map(F, inp['weights']))
        qs = inp['qs']
        atom_mass = [w / q**n for w, q in zip(weights, qs)]
        assert sum(weights) == 1
        assert sum(q**n * mass for q, mass in zip(qs, atom_mass)) == 1
        indices = []
        nr = sorted((r for r in rows if r['n'] == n), key=lambda r: r['point_index'])
        assert [r['point_index'] for r in nr] == list(range(8))
        for r in nr:
            old = saved['oracle_records'][r['point_index']]
            assert old['receiver_x'] == r['receiver_x']
            assert old['winner_R'] == r['saved_winner_R']
            assert old['M'] == r['saved_full_M']
            assert old['critical_exceedance'] == r['saved_E']
            assert saved['critical_tau'] == r['saved_tau']
            m, R, M, tau = map(F, [r['full_m'], r['saved_winner_R'], r['saved_full_M'], r['saved_tau']])
            assert M == m / R**n and (M > tau) == r['saved_E']
            Nk = [math.prod(counts) for counts in r['packet_capture_counts']]
            assert sum(count * mass for count, mass in zip(Nk, atom_mass)) == m > 0
            assert [F(x) for x in r['packet_mi']] == [count * mass for count, mass in zip(Nk, atom_mass)]
            components, grouped = [], {}
            total_share, budget = F(0), F(0)
            for k, (q, w, mass, count) in enumerate(zip(qs, weights, atom_mass, Nk)):
                component = {'component_index': k, 'q': q, 'weight': str(w), 'fixed_global_packet_count': q**n,
                             'single_packet_Mi': str(mass), 'captured_count_Nk': count,
                             'saved_per_axis_capture_counts': r['packet_capture_counts'][k],
                             'captured_mass': str(count * mass), 'captured_packet_fullness': '1' if count else None}
                if count:
                    p = mass / m
                    j, theta = band(p)
                    share, fee = count * p, count * theta
                    assert fee < share <= 2 * fee
                    component.update({'single_packet_share_p': str(p), 'j': j, 'theta_j': str(theta),
                                      'actual_share_Nk_p': str(share), 'theta_count': str(fee),
                                      'p_equals_upper_endpoint': p == 2 * theta})
                    g = grouped.setdefault(j, {'theta': theta, 'count': 0, 'share': F(0), 'components': []})
                    g['count'] += count
                    g['share'] += share
                    g['components'].append(k)
                    total_share += share
                    budget += fee
                else:
                    component.update({'single_packet_share_p': None, 'j': None, 'theta_j': None,
                                      'actual_share_Nk_p': '0', 'theta_count': '0', 'p_equals_upper_endpoint': None})
                components.append(component)
            bands = []
            for j, g in sorted(grouped.items()):
                fee = g['theta'] * g['count']
                assert fee < g['share'] <= 2 * fee
                bands.append({'j': j, 'theta_j': str(g['theta']), 'captured_packet_count': g['count'],
                              'component_indices': g['components'], 'actual_share_sum': str(g['share']),
                              'theta_count': str(fee)})
            assert total_share == 1 and F(1, 2) <= budget < 1
            assert sum(F(g['theta_count']) for g in bands) == budget
            indices.append(len(result))
            result.append({'n': n, 'point_index': r['point_index'], 'source_seed': r['source_seed'],
                           'probe_type': r['probe_type'], 'receiver_x': r['receiver_x'], 'winner_R': str(R),
                           'M': str(M), 'tau': str(tau), 'E': r['saved_E'], 'full_m': str(m),
                           'components': components, 'merged_bands': bands, 'total_captured_packet_count': sum(Nk),
                           'sum_Nk_Mi': str(m), 'sum_Nk_p': str(total_share),
                           'row_theta_count_sum': str(budget), 'row_lower_slack': str(budget - F(1, 2)),
                           'row_upper_slack': str(1 - budget), 'E_gated_share_sum': str(total_share if r['saved_E'] else 0),
                           'E_gated_theta_count_sum': str(budget if r['saved_E'] else 0), 'all_exact_checks_passed': True})
        rr = [result[i] for i in indices]
        er = [r for r in rr if r['E']]
        jj = [g['j'] for r in rr for g in r['merged_bands']]
        values = [F(r['row_theta_count_sum']) for r in rr]
        ev = [F(r['row_theta_count_sum']) for r in er]
        rounds.append({'n': n, 'record_indices': indices, 'receiver_count': len(rr), 'E_receiver_count': len(er),
                       'occupied_j_min': min(jj), 'occupied_j_max': max(jj),
                       'distinct_j': sorted(set(jj)), 'occupied_band_count_per_row': [len(r['merged_bands']) for r in rr],
                       'row_budget_min': str(min(values)), 'row_budget_max': str(max(values)),
                       'row_budget_min_float': float(min(values)), 'row_budget_max_float': float(max(values)),
                       'E_row_budget_min_float': float(min(ev)) if ev else None,
                       'E_row_budget_max_float': float(max(ev)) if ev else None,
                       'captured_component_incidence_count': sum(bool(c['captured_count_Nk']) for r in rr for c in r['components']),
                       'observed_exact_upper_endpoint_count': sum(c['p_equals_upper_endpoint'] is True for r in rr for c in r['components']),
                       'all_exact_checks_passed': True})
    out = {'status': 'complete', 'registration_sha256': digest(regpath), 'script_sha256': digest(Path(__file__)),
           'inputs_sha256': reg['inputs_sha256'], 'rounds': rounds, 'records': result,
           'oracle_rerun': False, 'new_samples': False, 'exact_fraction_arithmetic': True,
           'global_packet_mass_W': '1', 'new_source_normalization': False,
           'finite_deterministic_row_audit_only': True, 'elapsed_seconds': perf_counter() - start,
           'conclusion': 'All complete-source single-atom packet rows have sum Np=1 and 1/2<=sum theta N<1. No source integral or cross-band square inference.'}
    target = BASE / (PREFIX + '_results.json')
    with target.open('x') as f:
        json.dump(out, f, ensure_ascii=False, indent=2)
        f.write('\n')
    print(json.dumps({'status': out['status'], 'rounds': rounds, 'results_sha256': digest(target),
                      'elapsed_seconds': out['elapsed_seconds']}, ensure_ascii=False, indent=2))


if __name__ == '__main__':
    main()
