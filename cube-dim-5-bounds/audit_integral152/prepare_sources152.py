"""Read-only finite identity checks, no new factor/mean/old power."""
from collections import Counter
from datetime import datetime, timezone
from fractions import Fraction
from math import factorial
from pathlib import Path
import gzip
import hashlib
import json
import os
import time

HERE = Path(__file__).resolve().parent
MAIN = Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube-dim-5/cube-dim-5-bounds')
AUTHOR = Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube5-curvature128-sol/cube-dim-5-bounds/upper_integral150')


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read(path):
    path = Path(path)
    raw = path.read_bytes()
    return json.loads(gzip.decompress(raw) if path.suffix == '.gz' else raw)


def write(name, obj):
    with (HERE / name).open('x') as handle:
        json.dump(obj, handle, sort_keys=True, indent=2)
        handle.write('\n')


def main():
    start = time.monotonic()
    pre = read(AUTHOR / 'preflight_results.json')
    manifest = read(AUTHOR / 'manifest.json')
    locks = {str(AUTHOR / 'manifest.json'): sha(AUTHOR / 'manifest.json')}
    for name, row in manifest['files'].items():
        locks[str(AUTHOR / name)] = row['sha256']
    for obj in list(pre['source_locks'].values()) + list(pre['selected_packet_locks'].values()):
        assert locks.get(obj['path'], obj['sha256']) == obj['sha256']
        locks[obj['path']] = obj['sha256']
    extra = ['upper_alpha1_029/report.md', 'upper_pathmean111/report.md',
             'upper_pathmean111/pathmean.py', 'upper_pathmean111/secant.py',
             'upper_pathmean111/manifest.json', 'upper_path112/report.md',
             'upper_path112/production_v2.py', 'upper_path112/secant.py',
             'upper_path112/pathmean.py', 'audit_path121/manifest.json',
             'audit_path121/report.md', 'audit_path121/whole_results.json']
    for name in extra:
        locks[str(MAIN / name)] = sha(MAIN / name)
    write('source_registration.json', {
        'task': 'audit_integral152_readonly_source_identity', 'pid': os.getpid(),
        'utc': datetime.now(timezone.utc).isoformat(), 'code_sha256': sha(__file__),
        'finite_selected_packets': 64, 'old_source_nodes': 97972,
        'new_factor_mean_H_R_or_old_power_calls': 0, 'source_locks': locks})
    for path, h in locks.items():
        assert sha(path) == h, path
    assert locks[str(AUTHOR / 'manifest.json')] == '1e404882a7a47459eeddac53cfd2985036e8fb3a409293db1931c561e0a85127'
    old = MAIN / 'upper_path112'
    rows, certs = read(old / 'nodes.json.gz'), read(old / 'node_certificates.json.gz')
    tree, geos = read(old / 'numeric_hierarchy.json.gz'), read(old / 'geometries.json.gz')
    reg, selected = read(old / 'scalar_registration.json'), read(AUTHOR / 'selected64.json')
    rm, cm = {r['key']: r for r in rows}, {r['key']: r for r in certs}
    ns = tree['nodes']
    assert len(rm) == len(cm) == len(ns) == 97972
    assert len(tree['frontier']) == 94580 and len(selected) == 64
    assert pre['time_endpoints'] == reg['time_endpoints']
    times = reg['time_endpoints']
    assert len(times) == 73 and times[0] == '1/512'
    assert all(Fraction(a) < Fraction(b) for a, b in zip(times, times[1:]))
    inherited = read(MAIN / 'audit_path121/manifest.json')
    assert inherited['status'] == 'complete_pass'
    assert inherited['source112_manifest_sha256'] == locks[str(old / 'manifest.json')]
    assert inherited['nodes'] == 97972 and inherited['volume'] == '1'
    assert tree['budget']['UH'] == pre['saved_complete_UH']
    eligible = [k for k in tree['frontier']
                if cm[k]['new_selfU'] < cm[k]['old_selfU']
                and rm[k]['old_U_i'] >= 1 << 144
                and rm[k]['old_L_i'] < rm[k]['old_U_i']]
    expected = sorted(eligible, key=lambda k: (
        -(ns[k]['selfU'] - ns[k]['selfL']), ns[k]['level'], ns[k]['indices']))[:64]
    assert [x['key'] for x in selected] == expected
    for rank, item in enumerate(selected):
        key = item['key']
        r, c, n = rm[key], cm[key], ns[key]
        assert item['rank'] == rank and item['source_node'] == r and item['source112_certificate'] == c
        assert key in tree['roots'] and key in tree['frontier']
        assert r['level'] == n['level'] == 5 and r['geometry'] == 'g002'
        assert n['parent'] is None and not n['children']
        assert n['effU'] == n['selfU'] == c['new_selfU'] == item['saved_selfU'] < c['old_selfU']
        assert n['effL'] == n['selfL'] == r['old_selfL'] == item['saved_selfL']
        assert c['old_selfU'] == r['old_selfU'] and c['old_M_i'] == r['old_M_i']
        assert 0 <= r['old_L_i'] <= c['new_M_i'] <= r['old_U_i']
        assert r['old_U_i'] >= 1 << 144 and r['old_L_i'] < r['old_U_i']
        weight = factorial(4)
        for mult in Counter(r['indices']).values():
            weight //= factorial(mult)
        assert r['weight'] == n['weight'] == weight
        assert r['indices'] == n['indices'] and Fraction(r['point'][0]) == 1
        offsets = [[Fraction(i, 1 << r['level']) - Fraction(p),
                    Fraction(i + 1, 1 << r['level']) - Fraction(p)]
                   for i, p in zip(r['indices'], r['point'][1:])]
        geometry = geos[r['geometry']]
        assert offsets == [[Fraction(x) for x in ab] for ab in geometry['offsets']]
        assert offsets == [[Fraction(x) for x in ab] for ab in geometry['profile']['offsets']]
        assert max(abs(x) for ab in offsets for x in ab) == Fraction(geometry['profile']['h'])
        assert max(abs((a+b)/2) for a, b in offsets) == Fraction(geometry['profile']['mu_infinity'])
        packet = read(r['packet_path'])
        assert packet['key'] == key and packet['point'] == r['point']
        assert [e['t'] for e in packet['endpoints']] == times
        assert sha(r['packet_path']) == r['packet_sha256'] == c['packet_sha256']
        for e in packet['endpoints']:
            for bound in ('lower', 'upper'):
                denominator = Fraction(e['H'][bound]).denominator
                assert denominator & (denominator - 1) == 0
    for path, h in locks.items():
        assert sha(path) == h, path
    write('input_locks.json', dict(sorted(locks.items())))
    write('source_check.json', {
        'status': 'selected64_original_active_chord_root_geometry_packet_identity_bound',
        'pid': os.getpid(), 'seconds': time.monotonic() - start,
        'selected': 64, 'original_nodes': 97972, 'eligible_saved_frontier': len(eligible),
        'same_time_endpoints': 73, 'geometry_ids': ['g002'], 'all_levels': [5],
        'all_unexpanded_roots': True, 'all_strict112_chord_arm_active': True,
        'all_U_at_least_1over65536': True, 'all_orbit_and_volume_identity_checked': True,
        'all_offsets_from_real_indices_and_anchor': True, 'dyadic_H_inputs': True,
        'source_locks_before_after_unchanged': len(locks),
        'new_factor_mean_H_R_or_old_power_calls': 0})
    print(json.dumps({'pid': os.getpid(), 'selected': 64, 'source_locks': len(locks),
                      'new_factor_mean_H_R_or_old_power_calls': 0}))


if __name__ == '__main__':
    main()
