#!/usr/bin/env python3
"""Pure postprocessing of the only B015 batch: no native calls/imports."""
import sys
sys.dont_write_bytecode = True
import hashlib
import itertools
import json
import os
from decimal import Decimal, localcontext
from fractions import Fraction as Q
from pathlib import Path
from time import perf_counter
import prepare as prep

HERE = Path(__file__).resolve().parent
T, M64, M625 = Q(6, 125), Q(3, 256000), Q(3, 250000)

def interval(obj):
    # Independently implemented exact endpoint decoding (do not import search).
    assert set(obj) == {'lower_mantissa', 'upper_mantissa', 'common_binary_exponent'}
    a, b, e = int(obj['lower_mantissa']), int(obj['upper_mantissa']), int(obj['common_binary_exponent'])
    assert a <= b and abs(e) <= 10000
    den = 1 << max(0, -e)
    mul = 1 << max(0, e)
    return Q(a*mul, den), Q(b*mul, den)

def decimal(q):
    with localcontext() as ctx:
        ctx.prec = 65
        return str(Decimal(q.numerator)/Decimal(q.denominator))

def s2_distance_squared(nums, den):
    x = [Q(n, den) for n in nums]
    vs = [(1, *s) for s in itertools.product((-1, 1), repeat=5)]
    norm2 = sum(t*t for t in x)
    candidates = []
    for i, v in enumerate(vs):
        xv = sum(t*a for t, a in zip(x, v))
        for j in range(i):
            w = vs[j]
            xw = sum(t*a for t, a in zip(x, w))
            g = sum(a*b for a, b in zip(v, w))
            candidates.append((norm2-(6*xv*xv-2*g*xv*xw+6*xw*xw)/(36-g*g), j, i))
    assert len(candidates) == 496
    d2, j, i = min(candidates)
    assert d2 > 0
    return {'distance_squared_exact': str(d2), 'closest_pair_indices': [j, i],
            'pairs_checked': 496, 'strictly_off_S2': True}

def main():
    tic = perf_counter()
    reg = json.loads((HERE/'preparation.json').read_text())
    release = json.loads((HERE/'release.json').read_text())
    state = json.loads((HERE/'checkpoint.json').read_text())
    assert prep.sha(HERE/'release.json') == '01d229c70878f3a3ce4d9a4a4236c0a8d358f2a661250b25cec5a552c62e6d25'
    for name, digest in release['prepared_file_sha256'].items():
        assert prep.sha(HERE/name) == digest
    assert prep.sha(release['C014_final_report_path']) == release['C014_final_report_sha256']
    for row in reg['inputs']:
        assert prep.sha(row['path']) == row['sha256']
    assert state['status'] == 'POINT_OBSTRUCTION_TO_STRICT_64_STOP_IMMEDIATELY'
    assert state['native_attempts'] == len(state['native_results']) == 2
    assert state['callback_requests'] == 2 and len(state['seed_runs']) == 1
    assert state['counterexample_native_index'] == 2
    assert state['elapsed_seconds'] < 1200 and state['all_spawned_native_processes_waited']
    assert not state['global_minimum_certified']
    _, rays = prep.fan_data()
    forms = prep.all_zero_forms(rays)
    results = []
    keys = set()
    for rec, proposal in zip(state['native_results'], state['proposal_records']):
        assert rec['returncode'] == 0 and rec['status'] == 'POINT_INTERVALS_VALIDATED'
        assert prep.sha(rec['stdout_path']) == rec['stdout_sha256']
        assert prep.sha(rec['stderr_path']) == rec['stderr_sha256']
        assert Path(rec['stderr_path']).read_text() == ''
        raw = json.loads(Path(rec['stdout_path']).read_text())
        assert raw['status'] == 'generic_nonzero_supported'
        assert raw['precision_bits'] == 512
        assert raw['full_group_terms'] == 967680 and raw['computed_sign_pairs'] == 483840
        assert raw['zero_node_terms'] == 0 and raw['confluent_terms'] == 0
        assert raw['unique_node_values'] == 16931
        ns, den = rec['numerators'], rec['denominator']
        assert len(ns) == 6 and len(raw['gradient']) == 6
        assert proposal['evaluated_numerators'] == ns and proposal['evaluated_denominator'] == den
        key = tuple(Q(n, den) for n in ns)
        assert key not in keys
        keys.add(key)
        pf = prep.preflight(ns, den, forms, rays)
        assert pf['feasible'] and pf['int64_safe'] and pf['zero_form_count'] == 0
        center = [Q(n, reg['seeds'][0]['denominator']) for n in reg['seeds'][0]['numerators']]
        assert max(abs(Q(n, den)-c) for n, c in zip(ns[1:], center[1:])) <= Q(1, 20)
        k, p = interval(raw['K']), interval(raw['p_T'])
        grad = [interval(v) for v in raw['gradient']]
        eu = interval(raw['Euler_x_grad_plus5K_plusTp'])
        assert k[0] > 0 and p[0] > 0 and eu[0] <= 0 <= eu[1]
        assert k == tuple(Q(v) for v in rec['check']['K_interval'])
        elo, ehi = 5*k[0]+T*p[0], 5*k[1]+T*p[1]
        for n, pair in zip(ns, grad):
            assert n > 0
            elo += Q(n, den)*pair[0]
            ehi += Q(n, den)*pair[1]
        assert elo <= 0 <= ehi
        for pair in [k, p, *grad]:
            assert pair[1]-pair[0] < max(Q(1), abs(sum(pair)/2))*Q(2)**-256
        # Simple decimal rational enclosures, proved by Fraction comparison.
        display_den = 10**16
        l = k[0].numerator*display_den//k[0].denominator
        u = -((-k[1].numerator*display_den)//k[1].denominator)
        assert Q(l, display_den) < k[0] <= k[1] < Q(u, display_den)
        generic_radius = min(Q(abs(sum(a*n for a, n in zip(f, ns))), 2*den*sum(map(abs, f))) for f in forms)
        assert generic_radius > 0
        q = (k[0]/T, k[1]/T)
        qprime = ((T*p[0]-k[1])/(T*T), (T*p[1]-k[0])/(T*T))
        free_directions = []
        for i in range(1, 6):
            a, b = grad[i]
            sign = 1 if b < 0 else (-1 if a > 0 else 0)
            assert sign != 0  # actual saved results strictly separate these signs
            xx = Q(ns[i], den)
            lower = max(Q(1, 10**8), center[i]-Q(1, 20))
            upper = min(1-Q(1, 10**8), center[i]+Q(1, 20))
            if sign == 1:
                gap_allowance = Q(ns[i-1]-ns[i], den)-Q(1, 10**8)
                local_allowance = upper-xx
                dd = (a, b)
            else:
                gap_allowance = (Q(ns[i]-ns[i+1], den) if i < 5 else xx)-Q(1, 10**8)
                local_allowance = xx-lower
                dd = (-b, -a)
            maximum = min(gap_allowance, local_allowance)
            assert maximum >= 0
            free_directions.append({'coordinate': i+1, 'signed_unit_direction': sign,
                'gradient_interval_exact': [str(a), str(b)], 'gradient_display_midpoint': decimal((a+b)/2),
                'directional_derivative_interval_exact': [str(v) for v in dd],
                'strict_downhill_derivative': dd[1] < 0,
                'axis_step_supremum_remaining_in_registered_local_domain': str(maximum),
                'step_supremum_display': decimal(maximum),
                'note': 'derivative only; no finite-step decrease or new point evaluated'})
        row = {'native_index': rec['index'], 'numerators': ns, 'denominator': den,
               'K_lower_exact': str(k[0]), 'K_upper_exact': str(k[1]),
               'K_display_midpoint': decimal(sum(k)/2),
               'coarse_rational_enclosure': [str(Q(l, display_den)), str(Q(u, display_den))],
               'K_width_exact': str(k[1]-k[0]), 'K_width_display': decimal(k[1]-k[0]),
               'threshold64': str(M64), 'threshold625': str(M625),
               'strictly_below_threshold64': k[1] < M64,
               'strictly_below_threshold625': k[1] < M625,
               'threshold64_minus_K_upper_exact': str(M64-k[1]),
               'threshold64_minus_K_upper_display': decimal(M64-k[1]),
               'point_cost_interval_exact': [str(T/(64*k[1])), str(T/(64*k[0]))],
               'point_cost_display_lower': decimal(T/(64*k[1])),
               'q_K_over_T_interval_exact': [str(v) for v in q],
               'q_K_over_T_midpoint_display': decimal(sum(q)/2),
               'dq_dT_interval_exact': [str(v) for v in qprime],
               'dq_dT_midpoint_display': decimal(sum(qprime)/2),
               'dq_dT_strict_sign': 'negative' if qprime[1] < 0 else ('positive' if qprime[0] > 0 else 'unresolved'),
               'free_coordinate_downhill_diagnostics': free_directions,
               'preflight': pf, 'open_generic_linf_radius': str(generic_radius),
               'S2_geometry': s2_distance_squared(ns, den),
               'Euler_exact_interval_contains_zero': True,
               'source_sha256': rec['source_sha256'], 'binary_sha256': rec['binary_sha256'],
               'raw_output_sha256': rec['stdout_sha256']}
        results.append(row)
    assert not results[0]['strictly_below_threshold64']
    assert results[1]['strictly_below_threshold64'] and results[1]['strictly_below_threshold625']
    # Reconstruct the one actual zero-node repair without any kernel replay.
    proposal = state['proposal_records'][1]
    assert proposal['initial_preflight']['zero_forms'] == [[1, -2, 1, 0, 0, 0]]
    assert len(proposal['repairs']) == 1 and proposal['repairs'][0]['k'] == 1
    repair = next(prep.gap_repairs(proposal['initial_numerators'], proposal['initial_denominator']))
    for key, value in repair.items():
        assert proposal['repairs'][0][key] == value
    assert repair['numerators'] == state['native_results'][1]['numerators']
    assert sum(repair['gap_delta_numerators']) == 0
    assert min(repair['candidate_gap_numerators']) >= 10**6
    out = {'status': 'POSTPROCESS_PASS_FIXED_T48_POINT_OBSTRUCTION',
           'pid': os.getpid(), 'elapsed_seconds': perf_counter()-tic,
           'additional_native_calls': 0, 'native_calls_in_registered_batch': 2,
           'first_64_obstruction_triggered_stop': True, 'second_seed_never_started': True,
           'results': results, 'actual_repair_reconstructed_exactly': True,
           'global_minimum_certified': False, 'C6_bound_improved': False,
           'proof_dependency': 'C014 native interval formula; parent reviewed, independent F engine audit remains separate',
           'source_hashes': [{'path': str(HERE/f), 'sha256': prep.sha(HERE/f)}
                             for f in ('preparation.json', 'release.json', 'checkpoint.json', 'check_results.py')]}
    prep.dump(HERE/'checked_results.json', out)
    print(json.dumps({k:out[k] for k in ('status','pid','elapsed_seconds','additional_native_calls')}, indent=2))
    for row in results:
        print(json.dumps({k:row[k] for k in ('native_index','K_display_midpoint','coarse_rational_enclosure',
               'K_width_display','threshold64_minus_K_upper_display','point_cost_display_lower','open_generic_linf_radius',
               'q_K_over_T_midpoint_display','dq_dT_midpoint_display','dq_dT_strict_sign')}, indent=2))
        print(json.dumps(row['free_coordinate_downhill_diagnostics'],indent=2))

if __name__ == '__main__':
    main()
