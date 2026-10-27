#!/usr/bin/env python3
"""Read-only byte integrity check. This is not a mathematical proof check."""
import hashlib
import json
from pathlib import Path


def main():
    project = Path(__file__).resolve().parents[1]
    manifest = json.loads((project/'rounds/001-manifest.json').read_text())
    for item in manifest['files']:
        path = project/item['path']
        assert path.is_file(), item['path']
        data = path.read_bytes()
        assert len(data) == item['bytes'], item['path']
        assert hashlib.sha256(data).hexdigest() == item['sha256'], item['path']
    print('PASS:', len(manifest['files']),
          'immutable round001 files; integrity only, not a proof check')


if __name__ == '__main__':
    main()
