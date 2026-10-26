"""Copy only frozen source manifests' explicitly selected repository payloads.

No mathematical experiment is run. Existing destinations must be byte identical.
Absolute external cache locks are retained as provenance, never copied as payload.
"""
from pathlib import Path
import argparse
import hashlib
import json
import os
import time

HERE = Path(__file__).resolve().parent
PROJECT = HERE.parent


def digest(data):
    return hashlib.sha256(data).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('registration')
    parser.add_argument('receipt')
    args = parser.parse_args()
    started = time.monotonic()
    registration = Path(args.registration)
    reg = json.loads(registration.read_text())
    payloads = []
    manifests = []
    for source in reg['sources']:
        directory = Path(source['directory'])
        manifest_path = directory / source['manifest']
        raw = manifest_path.read_bytes()
        assert digest(raw) == source['sha256'], str(manifest_path)
        manifest = json.loads(raw)
        selected = manifest[source['payload_field']]
        if isinstance(selected, dict):
            rows = [{'path': key, 'sha256': value if isinstance(value, str)
                     else value['sha256']} for key, value in selected.items()]
        else:
            rows = selected
        rows = [*rows, {'path': source['manifest'], 'sha256': source['sha256']}]
        for row in rows:
            path = Path(row['path'])
            path = path if path.is_absolute() else directory / path
            relative = path.relative_to(directory)
            assert '..' not in relative.parts
            data = path.read_bytes()
            assert digest(data) == row['sha256'], str(path)
            if 'bytes' in row:
                assert len(data) == row['bytes'], str(path)
            target = PROJECT / source['target'] / relative
            assert target.is_relative_to(PROJECT) and len(data) < 100_000_000
            payloads.append((target, data, str(path), row['sha256']))
        manifests.append({'source': str(manifest_path), 'sha256': source['sha256'],
                          'selected_files_including_manifest': len(rows)})
    # All sources verified before the first copy.
    for target, data, source, sha in payloads:
        if target.exists():
            assert target.read_bytes() == data, str(target)
        else:
            target.parent.mkdir(parents=True, exist_ok=True)
            target.write_bytes(data)
        assert digest(target.read_bytes()) == sha
    result = {'status': 'all_selected_frozen_source_payloads_byte_verified',
              'pid': os.getpid(), 'seconds': time.monotonic()-started,
              'registration_sha256': digest(registration.read_bytes()),
              'source_manifests': manifests,
              'files': [{'path': str(t.relative_to(PROJECT)), 'bytes': len(b),
                         'source': s, 'sha256': h} for t, b, s, h in payloads],
              'total_bytes': sum(len(b) for _, b, _, _ in payloads),
              'new_mathematical_evaluations': 0}
    receipt = Path(args.receipt)
    assert not receipt.exists()
    receipt.write_text(json.dumps(result, ensure_ascii=False, indent=2)+'\n')
    print(json.dumps({k: result[k] for k in ('status', 'pid', 'seconds', 'total_bytes')}))
    print('files', len(payloads))


if __name__ == '__main__':
    main()
