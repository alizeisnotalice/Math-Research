#!/usr/bin/env python3
"""Saved-certificate arithmetic/provenance only. Never imports an old program.

No geometry, union, old runner, MATLAB, network, or subprocess execution.
Python's arbitrary-precision integers and Fraction provide exact evidence.
"""
from pathlib import Path
from fractions import Fraction as Q
from decimal import Decimal, localcontext
from datetime import datetime, timezone
import hashlib
import json
import os
import platform
import sys
import time

HERE = Path(__file__).resolve().parent
N6 = Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube-dim-6/cube-dim-6-bounds')
EVIDENCE = Path('/Users/zhengzhihao/Documents/Codex/2026-10-27/n-n-5-1-1-verify/work/evidence-20261025/cube-finite-dim-bounds')
CHAIN = [
    'audit/round058/a2_next_event_review',
    'audit/round064/a2_one_more_event_review',
    'audit/round066/post_s2_local_event_review',
    'audit/round068/adjacent_hidden_chain_review',
    'audit/round070/first_overlap_cover_review',
    'audit/round073/post_s7_mask_review',
    'audit/round075/post_s8_overlap_review',
]


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def read(p):
    return json.loads(p.read_text())


def write(name, value):
    (HERE / name).write_text(json.dumps(value, ensure_ascii=False, indent=2) + '\n')


def dec(q):
    with localcontext() as ctx:
        ctx.prec = 70
        return str(Decimal(q.numerator) / Decimal(q.denominator))


def mapped(old):
    marker = '/cube-dim-6-bounds/'
    assert marker in old, old
    return N6 / old.split(marker, 1)[1]


def main():
    start = time.monotonic()
    cpu = time.process_time()
    sources = {}
    expected_checks = []
    historical_mutable_ledger = []

    def consume(p, purpose, expected=None, origin=None):
        p = p.resolve()
        assert p.is_file(), f'missing source: {p}'
        h = sha(p)
        entry = sources.setdefault(str(p), {'path': str(p), 'bytes': p.stat().st_size,
                                           'sha256': h, 'purposes': [], 'historical_references': []})
        assert entry['sha256'] == h, f'source changed within check: {p}'
        if purpose not in entry['purposes']:
            entry['purposes'].append(purpose)
        if origin and origin not in entry['historical_references']:
            entry['historical_references'].append(origin)
        if expected:
            assert h == expected, f'hash mismatch: {p}: {h} != {expected}'
            expected_checks.append({'path': str(p), 'expected': expected, 'actual': h})
        assert len(sources) <= 250
        return p

    # All geometric continuation reports are read as dependencies. Their saved
    # exact arithmetic and geometries are not replayed by this program.
    for stem in CHAIN:
        report = consume(N6 / stem / 'report.md', 'inherited geometric proof, read only')
        assert report.read_text().strip()
        result = consume(N6 / stem / 'review_results.json', 'inherited chain result')
        d = read(result)
        assert d['status'].startswith('pass_')
        for old, expected in d.get('input_sha256', {}).items():
            p = mapped(old)
            if p == N6 / 'results/bounds_ledger.json':
                # R058 captured a then-current mutable ledger. The current
                # ledger is separately captured below; byte identity is neither
                # expected nor falsely asserted across historical snapshots.
                historical_mutable_ledger.append({'historical_reference': old,
                    'historical_sha256': expected, 'rule': 'historical snapshot only; not a current-source hash requirement'})
                continue
            consume(p, 'hash-bound chain dependency', expected, old)
        manifest = N6 / stem / 'manifest.json'
        if manifest.exists():
            consume(manifest, 'chain artifact seal')
            m = read(manifest)
            for x in m.get('artifacts', []):
                consume(manifest.parent / x.get('path', x.get('name')),
                        'sealed chain artifact', x['sha256'])

    for rel in ['lower/round075/post_s8_overlap/manifest.json',
                'lower/round073/post_s7_mask_diagnosis/manifest.json']:
        p = consume(N6 / rel, 'author artifact seal')
        m = read(p)
        for x in m['artifacts']:
            consume(p.parent / x['path'], 'sealed author artifact', x['sha256'])
        for x in m['frozen_dependencies']:
            assert x['matches'] and x['frozen_sha256'] == x['current_sha256']
            consume(mapped(x['path']), 'frozen author dependency', x['frozen_sha256'], x['path'])

    for rel in ['audit/lower008/parametric_union/interface_proof.md',
                'audit/lower008/parametric_union/candidate_results.json',
                'audit/round008/lower_low_resource_review/report.md',
                'lower/round075/primary/registration.json',
                'lower/round075/primary/run_receipt.json',
                'audit/round091/cost_source_inventory/ERRATUM.md',
                'audit/round091/cost_source_inventory/erratum_comparison.json']:
        consume(N6 / rel, 'proof or source trace')
    consume(EVIDENCE / 'proofs/calibration_and_transfer.md', 'prior transfer proof, independently restated here')
    r4p = consume(EVIDENCE / 'baseline/original_certificates.json', 'R4 original interval')
    ledgerp = consume(N6 / 'results/bounds_ledger.json', 'current accepted lower slot and comparison intervals')
    dp = consume(N6 / 'lower/round075/post_s8_overlap/results.json', 'D075 exact candidate',
                 'dec62274abcf8204cf698cd7665dd6c754db3fee87a55beb5bfa1ece676e2817')
    ep = consume(N6 / 'audit/round075/post_s8_overlap_review/review_results.json', 'E075 independent certificate',
                 '03ee9356878432d22cf88c42889bddbc19a7bd021bb99da654249b088d0f0590')
    pp = consume(N6 / 'lower/round075/primary/results.json', 'primary saved-only arithmetic check',
                 'a1a1612ac7c0d0f7787e40e1dc457e46847cf762a5b09620b7c05fc37b3438f5')
    rp = consume(N6 / 'audit/round058/a2_next_event_review/review_results.json', 'saved positive root isolations')
    l, d, e, p, roots = map(read, (ledgerp, dp, ep, pp, rp))
    lower = l['lower']
    assert lower['id'] == 'R6_rectangular_exact_event_round075'
    assert lower['certificate_sha256'] == sha(dp)
    coefficients = e['new_ideal_ratio']['coefficients']
    assert len(coefficients) == 36
    assert set(coefficients) == {f'{i},{j}' for i in range(6) for j in range(6)}
    assert coefficients == d['candidate']['coefficients'] == lower['exact_representation']['coefficients']
    assert roots['root_denominator'] == e['root_denominator'] == 10**100
    ul, uh = map(Q, roots['root_enclosures']['u'])
    vl, vh = map(Q, roots['root_enclosures']['v'])
    assert 0 < ul < uh and ul**6 < 2 < uh**6
    assert 0 < vl < vh and vl**6 < 3 < vh**6
    assert uh-ul == vh-vl == Q(1, 10**100)
    lo = hi = Q(0)
    for key, value in coefficients.items():
        i, j = map(int, key.split(','))
        c = Q(value)
        ml, mh = ul**i * vl**j, uh**i * vh**j
        lo += c*(ml if c >= 0 else mh)
        hi += c*(mh if c >= 0 else ml)
    assert [lo, hi] == list(map(Q, e['new_ideal_ratio']['exact_interval']))
    assert [lo, hi] == list(map(Q, p['ideal_interval']))
    assert [lo, hi] == list(map(Q, lower['independent_tighter_ideal_interval']))
    dlo, dhi = map(Q, d['candidate']['exact_interval'])
    assert [dlo, dhi] == list(map(Q, lower['exact_interval']))
    assert dlo <= lo <= hi <= dhi

    # Saved structure assertions are sanity checks of supplied evidence, not
    # renewed geometry or a full union computation.
    assert e['all_simultaneous_ties'] == [[31, 2]]
    assert e['new_s8_right_germ_verified'] and e['all_127_minmax_choices_unchanged']
    assert e['complete_effective_groups'] == [[110,111,112,114,115,117], [0]]
    assert len(e['all_seven_tie_labels_effective']) == 7
    assert len(e['six_disjoint_projection_pairs']) == 6
    assert all(Q(x['gap']['exact_interval'][0]) > 0 for x in e['six_disjoint_projection_pairs'])
    assert len(e['new_complete_order_gaps']) == 36 and len(e['new_effective_gaps']) == 21
    assert e['next_ties'] == [[7,2]] and e['s9_left_state_equals_s8_right']
    assert e['new_endpoint_all_127_boxes_valid']
    assert e['resource']['full_union_calls'] == 0
    f = e['finite_L1']
    eps, theta = Q(f['epsilon']), Q(f['threshold'])
    k = f['receiver_index_radius']
    assert eps == Q(1,10**22) and Q(f['bump_halfwidth']) == eps/2
    assert theta == 1-Q(1,10**16) and k == 10**16
    assert f['source_index_radius'] == k+2 and f['max_relative_source_index'] == 2
    assert Q(f['minimum_new_box_width']['exact_interval'][0]) > 2*eps
    assert f['all_box_widths_below_one']
    loss = 64*127*12*eps
    inward = lo-loss
    boundary = Q(2*k+1,2*k+5)**6
    finite = theta*inward*boundary
    target = Q(2073015624830,10**12)
    assert loss == Q(f['union_loss_upper']) and inward == Q(f['inward_ratio_lower'])
    assert boundary == Q(f['boundary_factor'])
    assert finite == Q(f['weak_ratio_lower']) == Q(p['finite_weak_ratio_lower'])
    assert target == Q(f['strict_rational_target']) == Q(p['strict_target'])
    assert finite-target == Q(f['margin_over_target']) == Q(p['finite_margin'])
    old_upper = Q(e['complete_current073_upper'])
    assert finite > target > old_upper
    assert finite-old_upper == Q(f['margin_over_complete073_upper'])

    frozen = read(N6/'audit/round075/post_s8_overlap_review/frozen_073_baseline.json')
    assert old_upper == Q(frozen['comparison_upper']) == max(Q(v[1]) for v in frozen['complete_ideal_intervals'].values())
    r4 = next(x['lower_bound_interval'] for x in read(r4p) if x['dimension'] == 4)
    n5 = l['dynamic_inherited_n5']['interval']
    assert r4 == l['inherited_lower']['interval']
    erratum = read(N6/'audit/round091/cost_source_inventory/erratum_comparison.json')
    assert n5 == erratum['n5_complete_construction_ratio_interval']
    r4lo, r4hi = map(Q,r4)
    n5lo, n5hi = map(Q,n5)
    assert r4lo <= r4hi < n5lo <= n5hi < target < finite < lo <= hi
    comparison = {'R4_interval':r4, 'n5_R003_construction_interval':n5,
        'n5_status':'saved accepted n5 construction interval from n6 ledger and erratum; not an upper bound on C5',
        'R4_upper_lt_n5_lower':True, 'n5_upper_lt_R075_explicit_finite_target':True,
        'n5_upper_lt_R075_full_ideal_lower':True,
        'exact_positive_margins':{'n5_lower_minus_R4_upper':str(n5lo-r4hi),
            'target_minus_n5_upper':str(target-n5hi), 'R075_ideal_lower_minus_n5_upper':str(lo-n5hi)},
        'display_only':{'R4':[dec(r4lo),dec(r4hi)],'n5':[dec(n5lo),dec(n5hi)]}}
    baseline = {'dimension':7, 'status':'inherited_exact_lower_via_proved_dimension_monotonicity',
        'model':{'measure':'Lebesgue','cubes':'centered axis-parallel','radius':'one common r per cube','inputs':'general L1','scales':'all positive radii'},
        'new_seven_dimensional_construction':False,
        'bound_relation':'C7 >= C6 >= R075; C7 > 2073015624830/1000000000000',
        'source_lower_id':lower['id'],
        'R075':{'expression':'sum_{i=0}^5 sum_{j=0}^5 c[i,j]*u^i*v^j',
            'positive_roots':{'u':'u^6=2,u>0','v':'v^6=3,v>0'},
            'coefficients':coefficients, 'author_full_interval':d['candidate']['exact_interval'],
            'independent_tighter_interval':e['new_ideal_ratio']['exact_interval'],
            'root_enclosures':roots['root_enclosures'],
            'spacing':lower['exact_representation']['spacing'],
            'ideal_display_lower':dec(lo), 'ideal_display_upper':dec(hi),
            'source_sha256':{'author':sha(dp),'independent':sha(ep),'primary':sha(pp),'ledger':sha(ledgerp)}},
        'explicit_finite_L1':{'source_dimension':6,'parameters':{'epsilon':str(eps),'eta':str(eps/2),'theta':str(theta),'K':k},
            'union_loss_upper':str(loss),'inward_ratio_lower':str(inward), 'boundary_factor':str(boundary),
            'weak_ratio_lower':str(finite),'strict_target':'2073015624830/1000000000000',
            'margin_over_target':str(finite-target),'display_weak_ratio_lower':dec(finite),
            'C7_transfer':'long seventh-coordinate truncations converge from below; a sufficiently long finite truncation also exceeds strict_target'},
        'comparison':comparison,
        'proof_dependency_chain':['accepted E008 affine union','E058','E064','E066','E068','E070','E073','E075','primary075'],
        'verification_scope':'saved hashes, coefficient identity, exact positive-root brackets, outward interval arithmetic, saved finite-L1 formula and exact interval comparisons; inherited geometry, no geometry/union rerun'}
    write('baseline.json',baseline)
    for source in sources.values():
        assert sha(Path(source['path'])) == source['sha256'], 'source changed during check'
    write('source_manifest.json',{'status':'PASS','sources':list(sources.values()),
        'source_file_count':len(sources),'expected_hash_check_count':len(expected_checks),
        'expected_hash_checks':expected_checks,
        'historical_mutable_ledger_snapshots_not_asserted_current':historical_mutable_ledger,
        'path_remap_rule':'historical .../cube-dim-6-bounds/<rel> mapped read-only to current designated n6 root; every immutable expected hash checked',
        'source_modifications':0})
    receipt = {'status':'PASS','utc':datetime.now(timezone.utc).isoformat(),'pid':os.getpid(),
        'python':platform.python_version(),'wall_seconds':time.monotonic()-start,'cpu_seconds':time.process_time()-cpu,
        'command':[sys.executable,'-B',str(Path(__file__).resolve())],
        'check_sha256':sha(Path(__file__).resolve()),'registration_sha256':sha(HERE/'registration.json'),
        'counts':{'unique_source_files':len(sources),'immutable_expected_hash_checks':len(expected_checks),
            'positive_root_brackets':2,'R075_coefficients':36,'finite_formula_cases':1,'interval_comparisons':3,
            'geometry_recomputations':0,'union_calls':0,'old_runner_calls':0,'MATLAB_calls':0,'subagents':0},
        'wall_limit_seconds':None,'arithmetic':'exact Python integers/Fraction; Decimal display only',
        'outputs':{name:sha(HERE/name) for name in ['baseline.json','source_manifest.json']},
        'scope_note':'Historical execution logs are source documentation only, never this run evidence.'}
    write('receipt.json',receipt)
    print(json.dumps({'status':'PASS','source_files':len(sources),'immutable_hash_checks':len(expected_checks),
        'ideal_lower':dec(lo),'finite_lower':dec(finite),'finite_minus_target':dec(finite-target),
        'wall_seconds':receipt['wall_seconds'],'geometry_calls':0,'union_calls':0},ensure_ascii=False,indent=2))


if __name__ == '__main__':
    main()
