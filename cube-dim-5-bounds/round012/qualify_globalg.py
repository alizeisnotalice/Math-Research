#!/usr/bin/env python3
"""Bind the primary g acceptance to independently audited, frozen evidence.

No heat kernel or special function is evaluated. The analytic responsibility
is recorded separately; source hashes identify inputs, not mathematical proof.
"""
from pathlib import Path
from fractions import Fraction as F
import hashlib
import json
import os
import time

HERE = Path(__file__).resolve().parent
PROJECT = HERE.parent
EXTERNAL = PROJECT.parent.parent / 'Math-Research-cube5-lower-audit-cont' / 'cube-dim-5-bounds' / 'audit_kernel081'


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b''):
            h.update(chunk)
    return h.hexdigest()


def read(path):
    return json.loads(Path(path).read_text())


def main():
    start = time.monotonic()
    output = HERE / 'globalg_acceptance.json'
    assert not output.exists(), 'completed acceptance must not be overwritten'
    reg = read(HERE / 'globalg_acceptance_registration.json')
    locked = reg['inputs']
    paths = {
        'numeric081_sha256': EXTERNAL / 'numeric_result.json',
        'qualification081_sha256': EXTERNAL / 'qualification_result.json',
        'v2_checker084_sha256': PROJECT / 'upper_portable084/check_portable_tree_v2.py',
        'v2_result084_sha256': PROJECT / 'upper_portable084/package_001/independent_check_v2_replay004.json',
        'theory083_report_sha256': PROJECT / 'audit_globalg083/report.md',
    }
    for key, path in paths.items():
        assert sha(path) == locked[key], (key, path)
    joined = read(paths['qualification081_sha256'])
    assert joined['status'] == 'all_values_structure_and_closed_cube_g_qualified_given_locked075_029'
    assert joined['final_state_sha256'] == locked['final070_sha256']
    assert joined['numeric_points'] == 1318 and joined['independent_components'] == 6590
    assert joined['H_epsilon_homogeneity_used'] is False
    assert joined['same_alpha1_L1_all_positive_scale_qualification'].startswith('locked029')
    # Bind all proof, structural, numeric and actual final-state sources.
    for path, digest in joined['input_locks'].items():
        assert sha(path) == digest, path
    numeric = read(paths['numeric081_sha256'])
    assert numeric['status'] == 'numeric_complete_all1318_H_plus4grad'
    assert numeric['verified_count'] == len(numeric['points']) == 1318
    assert numeric['shared_production_code'] is False
    assert numeric['shared_Arb_or_log_library'] is False
    manifest = read(PROJECT / 'upper_cover070/manifest.json')
    last = manifest['batches'][-1]
    overlay_path = Path(last['directory']) / 'state.json'
    assert sha(overlay_path) == last['files']['state.json']
    cache = read(overlay_path)['new_cache']
    assert len(cache) == 1318
    points = set()
    components = 0
    for index, item in enumerate(numeric['points']):
        assert item['index'] == index and item['all5_strict_containment']
        point = tuple(item['point'])
        assert point not in points
        points.add(point)
        independent = item['independent_values_512_dyadic']
        original = cache[','.join(point)]['values']
        assert len(independent) == len(original) == 5
        for (low, high), (a, b) in zip(independent, original):
            # This primary check strengthens the numerical audit's <= test.
            assert F(a) < F(low) <= F(high) < F(b)
            components += 1
        assert item['geometry']['all_cones'] == 3504
        assert item['geometry']['distinct_projections'] == 280
    portable = read(paths['v2_result084_sha256'])
    assert portable['status'] == 'PASS_portable_closed_root_certificate_package'
    assert portable['root_ranges_contiguous_and_exact']
    assert portable['canonical_root_catalog_count'] == portable['roots'] == 52360
    assert portable['unique_leaf_nodes_actual_set'] == portable['unique_leaf_anchors_actual_set'] == 97960
    assert portable['unique_tree_nodes_actual_set'] == 143560
    assert portable['fixed_log_independent_checks'] == 97960
    bundle = PROJECT / 'upper_portable084/package_001'
    assert sha(bundle / 'manifest.json') == portable['bundle_manifest_sha256']
    chunks = read(bundle / 'manifest.json')['chunks']
    assert len(chunks) == 26
    for chunk in chunks:
        assert sha(bundle / chunk['file']) == chunk['sha256']
    old_I = F(71555899, 10**12)
    g = F(1, 10000)
    assert 3125 * old_I**4 > 32**4 * g**5
    result = {
        'status': 'primary_accepts_complete_alpha1_cube_g_component',
        'g': str(g),
        'canonical_G': 'G(x)>g*||x||infinity^-4 for x!=0; G(0)=infinity',
        'new_independent_components_strictly_inside_production_balls': components,
        'closed_roots': 52360,
        'closed_leaves': 97960,
        'portable_v2_G2_repair_accepted': True,
        'all_26_portable_chunks_hash_bound': True,
        'analytical_review': 'main_review_pending.md and locked075/029/083',
        'shared_inputs_limit': '017 fan and inherited old H intervals; new1318 use distinct integer partial-fraction path',
        'complete_I_cap': str(old_I),
        'paired_cost_strictly_greater_than32': True,
        'formal_C5_upper_bound': '32',
        'formal_C5_bounds_changed': False,
        'post_dual_rounds': 0,
        'new_kernel_or_special_function_calls': 0,
        'pid': os.getpid(),
        'seconds': time.monotonic() - start,
        'source_sha256': sha(__file__),
        'registration_sha256': sha(HERE / 'globalg_acceptance_registration.json'),
        'evidence_sha256': {key: sha(path) for key, path in paths.items()},
    }
    output.write_text(json.dumps(result, ensure_ascii=False, indent=2) + '\n')
    print(json.dumps(result, ensure_ascii=False))


if __name__ == '__main__':
    main()
