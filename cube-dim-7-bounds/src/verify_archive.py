#!/usr/bin/env python3
"""Read-only integrity of an immutable round archive; not a proof check."""
import argparse
import hashlib
import json
from pathlib import Path


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--round', type=int, required=True)
    args = parser.parse_args()
    assert args.round > 0
    project = Path(__file__).resolve().parents[1]
    manifest = json.loads((project/f'rounds/{args.round:03d}-manifest.json').read_text())
    for item in manifest['files']:
        path = (project/item['path']).resolve()
        assert path.is_relative_to(project) and path.is_file(), item['path']
        data = path.read_bytes()
        assert len(data) == item['bytes'], item['path']
        assert hashlib.sha256(data).hexdigest() == item['sha256'], item['path']
    print('PASS:', len(manifest['files']), 'round', args.round,
          'immutable file hashes; integrity only, not a proof check')


if __name__ == '__main__':
    main()
