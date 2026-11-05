"""Portable rescaling audit; adapted from the independent review script.

Default: verify historical C4 coverage/provenance, rerun C1/C2, independently
integrate both kernels, and transfer C4 by exact positive scaling. No C4 rerun.
Source/adaptation details: provenance/UPSTREAM.json. Outputs are not proof inputs.
"""
from pathlib import Path
from fractions import Fraction as F
from hashlib import sha256
from datetime import datetime, timezone
import argparse
import json
import sys
import time

if not __debug__:
    raise SystemExit("Run without -O: mathematical assertions must be enabled.")

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / 'vendor'))
from flint import arb, fmpq, ctx
from hl_core import Spec
from hl_check import check_c12, c4_cells, provenance, prov_digest

ctx.prec = 192


def sha(path):
    return sha256(path.read_bytes()).hexdigest()


def aq(q):
    q = F(q)
    return arb(fmpq(q.numerator, q.denominator))


def exact_end(x, which):
    return F(str(getattr(x, which)().fmpq()))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, default=ROOT / 'runs' / 'latest-verification.json')
    args = parser.parse_args()
    t0 = time.monotonic()

    # Package integrity, separate from the old 13-file freeze.
    manifest_path = ROOT / 'package-manifest.json'
    package_files = json.loads(manifest_path.read_text())['files']
    for rel, digest in package_files.items():
        p = (ROOT / rel).resolve()
        assert p.is_relative_to(ROOT)
        assert sha(p) == digest, rel
    assert args.out.resolve() not in {(ROOT / p).resolve() for p in package_files}
    assert args.out.resolve() != manifest_path

    bindings = json.loads((ROOT / 'provenance/subset-bindings.json').read_text())
    freeze_path = ROOT / 'provenance/original-rescaling-freeze.json'
    assert sha(freeze_path) == bindings['original_13_file_freeze_sha256']
    assert sha(freeze_path) == '02e3e0b08d7d0ad1b9dced57f1678b21d783a704fd633d2ac1e7a8c6fbae312b'
    old_freeze = json.loads(freeze_path.read_text())
    for item in bindings['copied_frozen_subset']:
        assert old_freeze['files'][item['original_relative_path']] == item['sha256']
        assert sha(ROOT / item['package_relative_path']) == item['sha256']
    for name, digest in bindings['vendor_source_hashes'].items():
        assert sha(ROOT / 'vendor' / name) == digest

    oldpath = ROOT / 'certificates/original.json'
    newpath = ROOT / 'certificates/scaled.json'
    old = json.loads(oldpath.read_text())
    new = json.loads(newpath.read_text())
    kappa = F(10000993, 10000000)
    claim = F(3590867754339, 10**12)
    assert sha(oldpath) == bindings['original_certificate_sha256']
    assert F(new['c']) * kappa == F(old['c'])
    assert set(new['a']) == set(old['a'])
    assert len(old['a']) == 1181
    assert all(F(new['a'][key]) * kappa == F(value) for key, value in old['a'].items())
    assert set(old) - {'c', 'a', 'bound', 'note'} == set(new) - {'c', 'a', 'bound', 'note'}
    for key in set(old) - {'c', 'a', 'bound', 'note'}:
        assert new[key] == old[key]
    assert F(new['bound']) == claim

    oldspec, newspec = Spec(oldpath), Spec(newpath)
    recorded = json.loads((ROOT / 'records/original-full-report.json').read_text())
    old_prov = recorded['provenance']
    assert recorded['spec_sha256'] == sha(oldpath)
    assert recorded['args']['only'] == 'c1,c2,c3,c4'
    assert recorded['args']['cells'] == ''
    assert all(recorded[key]['ok'] for key in ('c1', 'c2', 'c3', 'c4'))
    assert recorded['c4']['prov_consistent']
    # Old source content is bound. Old absolute paths and binary digests are
    # historical data, not requirements imposed on this machine.
    assert old_prov['hl_check.py'] == sha(ROOT / 'vendor/hl_check.py')
    assert old_prov['hl_core.py'] == sha(ROOT / 'vendor/hl_core.py')
    results = recorded['c4']['results']
    expected = set(c4_cells(oldspec))
    assert len(results) == len(expected) == recorded['c4']['cells'] == 1368
    assert {(r['p'], r['q']) for r in results} == expected
    assert all(r['ok'] and r['nfail'] == 0 and not r['fails'] for r in results)
    assert all(r['spec_sha'] == sha(oldpath) and r['prov'] == prov_digest(old_prov) for r in results)
    assert sum(r['leaves'] for r in results) == recorded['c4']['leaves'] == 443754
    assert set(c4_cells(newspec)) == expected

    saved_scaled = json.loads((ROOT / 'records/scaled-c123-report.json').read_text())
    assert saved_scaled['spec_sha256'] == sha(newpath)
    assert saved_scaled['args']['only'] == 'c1,c2,c3'
    assert all(saved_scaled[key]['ok'] for key in ('c1', 'c2', 'c3'))
    assert F(saved_scaled['claim']) == claim

    actual_prov = provenance()
    assert actual_prov['hl_check.py'] == sha(ROOT / 'vendor/hl_check.py')
    assert actual_prov['hl_core.py'] == sha(ROOT / 'vendor/hl_core.py')
    # Deliberately no recorded_provenance == actual_provenance assertion.
    checks = {}
    for name, spec in [('original', oldspec), ('scaled', newspec)]:
        r = check_c12(spec, maxdepth=10)
        assert not r['c1_fail'] and not r['c2_fail']
        assert r['c1_cells'] == 1368 and r['c2_cells'] == 272
        lo = exact_end(r['c2_min'], 'lower')
        assert lo > (kappa if name == 'original' else F(1))
        # Direct orbit expansion and integration, not check_c3 or printed balls.
        n = int(spec.N)
        total = sum((spec.coef(i, j) for i in range(-n, n+1) for j in range(-n, n+1)), F(0))
        alpha = aq(spec.alpha)
        mass = aq(spec.c) * alpha * ((2-alpha) * aq(spec.R).log()).exp() / (2-alpha)
        mass += aq(spec.h**2 * total) / 2
        checks[name] = {
            'c1_cells': r['c1_cells'], 'c1_leaves': r['c1_leaves'],
            'c2_cells': r['c2_cells'], 'c2_leaves': r['c2_leaves'],
            'c2_lower_exact': str(lo), 'c2_lower_display': r['c2_min'].str(40),
            'coefficient_orbit_sum_exact': str(total),
            'half_mass_lower_exact': str(exact_end(mass, 'lower')),
            'half_mass_upper_exact': str(exact_end(mass, 'upper')),
            'half_mass_display': mass.str(45)}
        if name == 'scaled':
            assert exact_end(mass, 'upper') < claim
            assert exact_end(mass, 'lower') > F('3.590867754338446537')
            assert exact_end(mass, 'upper') < F('3.590867754338446538')
            checks[name]['claim_gap_lower_exact'] = str(claim-exact_end(mass, 'upper'))
        else:
            assert exact_end(mass, 'upper') < F(359123, 100000)

    result = {
        'status': 'PASS', 'utc': datetime.now(timezone.utc).isoformat(),
        'wall_seconds': time.monotonic()-t0, 'precision_bits': ctx.prec,
        'kappa_exact': str(kappa), 'claim_exact': str(claim),
        'original_sha256': sha(oldpath), 'scaled_sha256': sha(newpath),
        'original_report_sha256': sha(ROOT / 'records/original-full-report.json'),
        'package_manifest_sha256': sha(manifest_path), 'package_files_verified': len(package_files),
        'original_freeze_subset_files_verified': len(bindings['copied_frozen_subset']),
        'original_full_13_file_freeze_reverified': False,
        'c4_historical_record_cells_checked': len(expected),
        'c4_historical_record_accepted_leaves': 443754,
        'c4_numeric_leaf_bounds_rerun': False,
        'c4_summary_record_is_not_a_per_leaf_certificate': True,
        'c4_new_kernel': 'Inherited by exact positive scaling, not numerically rerun.',
        'checks': checks, 'historical_provenance': old_prov,
        'actual_environment_provenance': actual_prov,
        'historical_and_actual_provenance_required_equal': False,
        'lean_build_executed': False}
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
