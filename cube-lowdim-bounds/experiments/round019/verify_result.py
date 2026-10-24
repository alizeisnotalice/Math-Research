#!/usr/bin/env python3
"""Read-only saved-record verification. Never imports search or evaluates geometry."""
import argparse
import hashlib
import json
import math
from fractions import Fraction
from pathlib import Path
from time import monotonic

ROOT = Path(__file__).resolve().parent
SEARCH_HASH = '8ebed200d5c1e3ae29799df3440f23c30cbd4eedbefa4671a3cd3dc3787b0851'
PHI_HASH = 'baef432fb4a37cc4c8199775c226b8782259bb774eaf92bde6af3d63397ad286'
BOUNDS = [(-.2,.2),(-.2,.2),(-.15,.15),(-.15,.15),(-.3,.3)]

def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def verify(trace=None):
    start = monotonic()
    r = json.loads((ROOT/'result.json').read_text())
    p = json.loads((ROOT/'phi_reference.json').read_text())
    assert digest(ROOT/'search.py') == r['search_sha256'] == SEARCH_HASH
    assert digest(ROOT/'phi_reference.py') == r['phi_reference']['source_sha256'] == PHI_HASH
    assert r['status'] == 'finite_numerical_run_closed' and r['failure'] is None
    assert r['search_stop_reason'] == 'optimizer_returned'
    assert r['bounds'] == [list(b) for b in BOUNDS]
    assert r['seed'] == 2026102419
    lo, hi = Fraction(p['lower']), Fraction(p['upper'])
    assert lo < hi and hi-lo == Fraction(p['width']) and 0 < hi-lo < Fraction(1,10**49)
    assert r['phi_reference']['lower'] == p['lower'] and r['phi_reference']['upper'] == p['upper']
    assert r['phi_reference']['float_display'] == float(hi)
    h,v,w = (r['Lin_parameters'][k] for k in ('h','V','w'))
    u = (2+math.sqrt(22))/3
    assert h == (1+u)/2 and v == h+1 and w == u*u-1
    mass = 2*(1+w)
    assert r['complete_cell_mass'] == mass
    expected = {'controls':2,'DE':640,'audits':2}
    for kind, n in expected.items():
        assert r['counts'][kind] == dict(started=n,completed=n,interrupted=0,failed=0)
    assert r['total_configurations_started'] == r['maximum_scheduled_configurations'] == 644
    assert r['hard_total_started_limit'] == 1000
    o = r['optimizer']
    assert (o['dimension'],o['popsize'],o['maxiter'],o['scipy_nfev'],o['scipy_nit']) == (5,8,15,640,15)
    assert o['polish'] is False and o['workers'] == 1 and o['scipy_success'] is False
    assert o['completed_generation_callbacks'] == 15 and o['x0_in_population_and_counted'] is True
    t = r['timing']
    assert 0 <= t['prepare_seconds'] <= t['search_end_seconds'] < t['prepare_controls_search_limit'] == 45
    assert t['search_end_seconds'] <= t['numerical_stage_seconds'] < t['total_numerical_limit'] == 60
    controls = r['controls']
    assert len(controls) == 2 and controls[0]['z'] == [0.]*5
    assert controls[1]['z'] == [.1,-.1,.1,-.1,.2]
    max_mass_error = max_score_error = 0.

    def objective_record(row):
        nonlocal max_mass_error,max_score_error
        assert row['status'] == 'completed' and len(row['z']) == 5
        assert all(a <= z <= b for z,(a,b) in zip(row['z'],BOUNDS))
        assert row['kind'] in ('controls','DE')
        assert 0 <= row['started_at_seconds'] <= row['completed_at_seconds'] <= t['search_end_seconds']
        assert row['period_mass'] == mass
        weights = [1.,w*(1+row['z'][4]),1.,w*(1-row['z'][4])]
        err = abs(math.fsum(weights)-mass)
        assert min(weights)>0 and err <= 1e-12*mass
        max_mass_error = max(max_mass_error,err)
        err = abs(row['score']-row['union_area']/row['period_mass'])
        assert err == 0
        max_score_error = max(max_score_error,err)
        g = row['geometry']
        assert g['periods'] == [2*h,2*v]
        assert g['patch_atoms_before_aggregation'] == g['patch_positions_after_aggregation'] == 36
        nx,ny = g['unique_x_coordinates'],g['unique_y_coordinates']
        assert g['endpoint_rectangles'] == nx*(nx+1)//2*(ny*(ny+1)//2)
        assert 0 < g['positive_unwrapped_boxes'] <= g['distinct_nonempty_closures'] <= g['endpoint_rectangles']
        assert 0 <= row['union_area'] <= 4*h*v+1e-9

    for row in controls+[r[k] for k in ('best_completed','best_nonzero_completed','best_DE_completed')]:
        objective_record(row)
    assert r['best_completed'] == controls[0]
    checks = r['independent_geometry_checks']
    assert len(checks) == 2
    sources = [controls[0],r['best_nonzero_completed']]
    for n,(c,s) in enumerate(zip(checks,sources),643):
        assert c['configuration'] == n and c['status'] == 'completed' and c['passed'] is True
        assert c['source_configuration'] == s['configuration'] and c['z'] == s['z']
        assert c['recorded_score'] == s['score']
        assert c['endpoint_cell_score'] == c['endpoint_cell_area']/c['independently_rebuilt_mass']
        assert c['absolute_score_difference'] == abs(c['endpoint_cell_score']-s['score']) <= c['tolerance'] == 1e-8
        assert c['mass_difference'] == abs(c['independently_rebuilt_mass']-s['period_mass'])
        assert c['mass_difference'] <= c['mass_tolerance'] == 1e-12*s['period_mass']
        assert t['search_end_seconds'] <= c['started_at_seconds'] <= c['completed_at_seconds'] <= t['numerical_stage_seconds']
    gain = r['best_completed']['score']-float(hi)
    assert gain == r['best_numerical_minus_Phi_upper_display']
    assert r['best_nonzero_completed']['score']-float(hi) == r['best_nonzero_numerical_minus_Phi_upper_display']
    assert r['numerical_improvement_candidate'] == (gain > r['numeric_gain_tolerance']) is False
    assert r['zero_control_possible_larger_witness_gain'] is False
    trace_report = None
    if trace is not None:
        assert digest(trace) == r['full_trajectory']['sha256']
        events = [json.loads(line) for line in trace.read_text().splitlines()]
        assert len(events) == 1288
        starts, ends = {},{}
        previous_time = -1.
        for e in events:
            n = e['configuration']
            target = starts if e['status']=='started' else ends
            assert n not in target
            assert e['status'] in ('started','completed')
            target[n] = e
            stamp = e['started_at_seconds'] if e['status']=='started' else e['completed_at_seconds']
            assert stamp >= previous_time
            previous_time = stamp
        assert sorted(starts) == sorted(ends) == list(range(1,645))
        for n in starts:
            assert all(ends[n][k]==val for k,val in starts[n].items() if k!='status')
        objectives = [ends[n] for n in range(1,643)]
        for e in objectives: objective_record(e)
        assert [e['kind'] for e in objectives[:2]] == ['controls']*2
        assert all(e['kind']=='DE' for e in objectives[2:])
        assert objectives[:2] == controls and [ends[643],ends[644]] == checks
        assert max(objectives,key=lambda e:e['score']) == r['best_completed']
        assert max((e for e in objectives if any(e['z'])),key=lambda e:e['score']) == r['best_nonzero_completed']
        assert max(objectives[2:],key=lambda e:e['score']) == r['best_DE_completed']
        de = objectives[2:]
        trace_report = {'events':len(events),'configurations':len(ends),'objective_records_checked':len(objectives),
                        'all_started_have_exactly_one_completed_terminal':True,'sha256':digest(trace),'bytes':trace.stat().st_size,
                        'DE_variable_min':[min(e['z'][j] for e in de) for j in range(5)],
                        'DE_variable_max':[max(e['z'][j] for e in de) for j in range(5)],
                        'DE_distinct_values':[len({e['z'][j] for e in de}) for j in range(5)]}
    return {'status':'saved_record_verification_passed','mode':'full_saved_trajectory' if trace is not None else 'archive_summary_records_only',
            'scope':'Saved arithmetic/ledger/hash/reference/check-record verification only; no geometry, union, or optimizer evaluated.',
            'geometry_evaluations':0,'optimizer_calls':0,'search_sha256':SEARCH_HASH,'phi_source_sha256':PHI_HASH,
            'result_sha256':digest(ROOT/'result.json'),'verifier_sha256':digest(Path(__file__)),
            'counts':r['counts'],'total_configurations':644,'all_trajectory_records_checked':trace is not None,
            'trace':trace_report,'maximum_seed_mass_identity_error':max_mass_error,'maximum_saved_area_over_mass_error':max_score_error,
            'saved_independent_geometry_records_checked':2,'independent_geometry_recomputed':False,
            'best_nonzero_score':r['best_nonzero_completed']['score'],'best_nonzero_minus_Phi_upper_display':r['best_nonzero_numerical_minus_Phi_upper_display'],
            'numerical_improvement_candidate':r['numerical_improvement_candidate'],'original_numerical_seconds':t['numerical_stage_seconds'],
            'verification_elapsed_seconds':monotonic()-start}

if __name__=='__main__':
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--trace',type=Path,help='Optional complete local JSONL; absent means summary only.')
    ap.add_argument('--output',type=Path,default=ROOT/'verification.json')
    args = ap.parse_args()
    result = verify(args.trace)
    args.output.write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps({'status':result['status'],'mode':result['mode'],'geometry_evaluations':0,'elapsed_seconds':result['verification_elapsed_seconds']}))
