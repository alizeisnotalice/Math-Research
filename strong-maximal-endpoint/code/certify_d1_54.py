#!/usr/bin/env python3
"""Portable launcher for the immutable r054 D1 core.

This wrapper performs path resolution and hash/race auditing only.  All
mathematical interval operations live in cert_d1_core54.py, whose SHA-256 is
recorded in the emitted JSON.
"""
from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import os
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
CORE = HERE / 'cert_d1_core54.py'

def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()

def root_at_or_above(start: Path) -> Path | None:
    start=start.resolve()
    for candidate in (start, *start.parents):
        if (candidate/'rounds/r052/fair52.json').is_file() and (candidate/'rounds/r030/receiver_seed.json').is_file():
            return candidate
    return None

def output_arg(rest: list[str], default: Path) -> Path:
    for i, item in enumerate(rest):
        if item == '--out' and i+1 < len(rest): return Path(rest[i+1]).resolve()
    return default

def main() -> None:
    parser=argparse.ArgumentParser(add_help=False)
    parser.add_argument('--project-root', type=Path)
    parser.add_argument('-h','--help', action='store_true')
    ours, rest=parser.parse_known_args()
    if ours.help:
        print('usage: certify_d1_54.py [--project-root PROJECT] [core options: --m M --workers W --out FILE]')
        return
    root=(ours.project_root.resolve() if ours.project_root else (root_at_or_above(Path.cwd()) or root_at_or_above(HERE)))
    if root is None:
        env=os.environ.get('R054_PROJECT_ROOT')
        root=root_at_or_above(Path(env)) if env else None
    if root is None: raise SystemExit('cannot locate project root; pass --project-root containing rounds/r052 and rounds/r030')
    fair=root/'rounds/r052/fair52.json'; seed=root/'rounds/r030/receiver_seed.json'
    start={'core_sha256':sha(CORE),'fair_sha256':sha(fair),'seed_sha256':sha(seed),'project_root':str(root)}
    out=output_arg(rest,Path('/tmp/r054_d1_certificate.json'))
    if '--out' not in rest: rest += ['--out',str(out)]
    if '--snapshot' not in rest: rest += ['--snapshot',str(CORE)]
    spec=importlib.util.spec_from_file_location('r054_core',CORE)
    if spec is None or spec.loader is None: raise SystemExit('cannot import immutable core')
    core=importlib.util.module_from_spec(spec)
    # The core uses a fork pool whose workers must be able to resolve the
    # defining module of lower_row during result serialization.
    sys.modules[spec.name]=core
    spec.loader.exec_module(core)
    core.REPO=root; core.FAIR=fair; core.SEED=seed; core.OUTDIR=out.parent
    argv_save=sys.argv; sys.argv=[str(CORE),*rest]
    core_exit=None
    try: core.main()
    except SystemExit as exc:
        # A deliberately coarse smoke grid writes a valid negative diagnostic
        # and then exits nonzero.  Still complete the wrapper's I/O/hash audit.
        core_exit=exc
    finally: sys.argv=argv_save
    end={'core_sha256':sha(CORE),'fair_sha256':sha(fair),'seed_sha256':sha(seed)}
    if {k:start[k] for k in end} != end: raise SystemExit('input changed during calculation; discard output and rerun')
    record=json.loads(out.read_text())
    if (record.get('source_sha256')!=start['core_sha256'] or
        record.get('fair_sha256')!=start['fair_sha256'] or
        record.get('seed_sha256')!=start['seed_sha256']):
        raise SystemExit('core result hashes do not match startup inputs')
    record['launcher_startup_hashes']=start
    record['launcher_final_hashes']=end
    record['launcher_scope']='Path resolution and hash/race audit only; no mathematical operations.'
    out.write_text(json.dumps(record,indent=2)+'\n')
    if core_exit is not None: raise core_exit

if __name__=='__main__': main()
