#!/usr/bin/env python3
"""Check archived hashes and exact saved comparisons; not a full proof checker."""
from pathlib import Path
from fractions import Fraction as Q
import argparse
import hashlib
import json

ROOT = Path(__file__).resolve().parent


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--round', type=int)
    args = parser.parse_args()
    manifests = sorted(ROOT.glob('audit/round*/primary/release_manifest.json'))
    manifest_path = (ROOT / f'audit/round{args.round:03d}/primary/release_manifest.json'
                     if args.round else manifests[-1])
    manifest = json.loads(manifest_path.read_text())
    for relative, digest in manifest['files'].items():
        path = ROOT / relative
        assert path.is_file(), relative
        assert hashlib.sha256(path.read_bytes()).hexdigest() == digest, relative
    ledger = json.loads((ROOT / manifest['ledger_snapshot']).read_text())
    lower, upper = ledger['lower'], ledger['upper']
    assert Q(lower['finite_L1_exact_lower']) > Q(ledger['opening']['lower']['construction_interval'][1])
    assert Q(lower['finite_L1_exact_lower']) > Q(lower['simple_strict_lower'])
    assert Q(lower['construction_interval'][0]) <= Q(lower['construction_interval'][1])
    assert Q(upper['F_interval'][0]) <= Q(upper['F_interval'][1]) < Q(upper['strict_rational_upper'])
    assert Q(upper['strict_rational_upper']) < Q(ledger['opening']['upper']['exact'])
    author = json.loads((ROOT / lower['certificate']).read_text())
    if 'trials' in author:
        row = next(t for t in author['trials'] if t.get('winner'))
    else:
        row = author
    assert row['union_interval'] == lower['construction_interval']
    assert row['finite_L1']['lower'] == lower['finite_L1_exact_lower']
    audited = json.loads((ROOT / lower['audit']).read_text())
    assert audited['status'] == 'PASS'
    assert audited['union_interval'] == lower['construction_interval']
    assert audited['finite_L1']['lower'] == lower['finite_L1_exact_lower']
    print(json.dumps({'status': 'PASS', 'round': manifest['round'],
                      'files': len(manifest['files']),
                      'scope': 'hashes and exact saved comparisons; full analytic proofs and union algorithms remain separate'}, indent=2))


if __name__ == '__main__':
    main()
