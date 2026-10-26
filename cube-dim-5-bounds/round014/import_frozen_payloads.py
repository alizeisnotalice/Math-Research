"""Byte-copy frozen whitelisted evidence, without executing its code.

Extends round013/import_sources.py to list manifests using the key `file`.
All bytes and destinations are checked before any copy. External caches stay
outside the repository; importing an artifact does not qualify its mathematics.
"""
import hashlib
import json
import os
from pathlib import Path
import sys
import time

PROJECT = Path(__file__).resolve().parent.parent


def digest(data):
    return hashlib.sha256(data).hexdigest()


def main():
    registration, receipt = map(Path, sys.argv[1:])
    assert not receipt.exists()
    started = time.monotonic()
    reg = json.loads(registration.read_bytes())
    payloads = []
    sources = []
    for spec in reg['sources']:
        directory = Path(spec['directory'])
        manifest_path = directory / spec['manifest']
        raw = manifest_path.read_bytes()
        assert digest(raw) == spec['sha256']
        manifest = json.loads(raw)
        selected = manifest[spec.get('payload_field', 'files')]
        if isinstance(selected, dict):
            rows = [{'path': key, 'sha256': value if isinstance(value, str)
                     else value['sha256']} for key, value in selected.items()]
        else:
            rows = selected
        rows = [*rows, {'path': spec['manifest'], 'sha256': spec['sha256']}]
        for row in rows:
            path = Path(row.get('path', row.get('file')))
            path = path if path.is_absolute() else directory / path
            relative = path.relative_to(directory)
            assert '..' not in relative.parts
            data = path.read_bytes()
            assert digest(data) == row['sha256']
            assert 'bytes' not in row or len(data) == row['bytes']
            assert len(data) < 100_000_000
            target = PROJECT / spec['target'] / relative
            assert target.is_relative_to(PROJECT)
            assert not target.exists() or target.read_bytes() == data
            payloads.append((target, data, str(path), row['sha256']))
        sources.append({'manifest': str(manifest_path),
                        'sha256': spec['sha256'], 'payload_count': len(rows)})
    for target, data, _, expected in payloads:
        if not target.exists():
            target.parent.mkdir(parents=True, exist_ok=True)
            target.write_bytes(data)
        assert digest(target.read_bytes()) == expected
    result = {'status': 'all_frozen_payloads_byte_verified', 'pid': os.getpid(),
              'seconds': time.monotonic() - started,
              'registration_sha256': digest(registration.read_bytes()),
              'source_manifests': sources,
              'new_mathematical_evaluations': 0,
              'total_bytes': sum(len(data) for _, data, _, _ in payloads),
              'files': [{'path': str(path.relative_to(PROJECT)),
                         'bytes': len(data), 'source': source, 'sha256': expected}
                        for path, data, source, expected in payloads]}
    receipt.write_text(json.dumps(result, ensure_ascii=False, indent=2) + '\n')
    print(json.dumps({key: result[key] for key in
                     ('status', 'pid', 'seconds', 'total_bytes')}))
    print('files', len(payloads))


if __name__ == '__main__':
    main()
