"""Verify frozen evidence in temporary copies. No tail MC or large cache."""
from pathlib import Path
import argparse
import hashlib
import json
import math
import os
import shutil
import subprocess
import sys
import tempfile
import time
import numpy as np

HERE = Path(__file__).resolve().parent


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def compare(a, b):
    if isinstance(a, dict):
        assert a.keys() == b.keys()
        for k in a:
            if k != 'seconds':
                compare(a[k], b[k])
    elif isinstance(a, list):
        assert len(a) == len(b)
        for x, y in zip(a, b):
            compare(x, y)
    elif isinstance(a, float):
        assert math.isclose(a, b, rel_tol=1e-11, abs_tol=1e-12), (a, b)
    else:
        assert a == b, (a, b)


def main():
    started = time.time()
    parser = argparse.ArgumentParser()
    parser.add_argument('--output')
    args = parser.parse_args()
    provenance = json.loads((HERE/'provenance.json').read_text())
    for file, digest in provenance['payload_SHA256'].items():
        assert sha(HERE/file) == digest, file
    # Registered components remain byte-for-byte intact.
    h = HERE/'round038_highdim'
    plan = json.loads((h/'HIGH_DIM_PLAN.json').read_text())
    for key, name in [('script_SHA256', 'probe_highdim.py'),
                      ('source_SHA256', 'source_tables.py'),
                      ('builder_SHA256', 'build_cache.py')]:
        assert sha(h/name) == plan[key]
    cond = HERE/'round038_conditional'
    saved = json.loads((cond/'RESULTS.json').read_text())
    assert saved['source_SHA256'] == sha(cond/'probe.py')
    assert saved['plan_SHA256'] == sha(cond/'PLAN.md')
    old = HERE.parent/'037/cell/probe.py'
    pb_reg = json.loads((HERE/'round038_pb_bound/registration.json').read_text())
    assert sha(old) == pb_reg['cube_source_SHA256']
    assert sha(HERE/'round038_pb_bound/probe.py') == pb_reg['script_SHA256']
    runs = []
    with tempfile.TemporaryDirectory(prefix='cube-round038-check-') as tmp:
        root = Path(tmp)
        for name in ['round038_conditional', 'round038_highdim', 'round038_pb_bound']:
            shutil.copytree(HERE/name, root/name)
        def run(relative, output_name):
            p = root/relative
            (p.parent/output_name).unlink()
            proc = subprocess.run([sys.executable, str(p)], cwd=root,
                                  env=dict(os.environ, PYTHONDONTWRITEBYTECODE='1'),
                                  capture_output=True, text=True)
            if proc.returncode:
                raise RuntimeError(proc.stdout+proc.stderr)
            data = json.loads((p.parent/output_name).read_text())
            runs.append({'script': relative, 'exit_code': 0})
            return data
        exact = run('round038_conditional/probe.py', 'RESULTS.json')
        compare(saved['cases'], exact['cases'])
        # Restore original bytes so postprocessing hash identity is preserved.
        shutil.copy2(cond/'RESULTS.json', root/'round038_conditional/RESULTS.json')
        layered = run('round038_conditional/layered_postprocess.py', 'LAYERED_RESULTS.json')
        original_layered = json.loads((cond/'LAYERED_RESULTS.json').read_text())
        compare(original_layered['cases'], layered['cases'])
        shutil.copy2(cond/'LAYERED_RESULTS.json', root/'round038_conditional/LAYERED_RESULTS.json')
        be = run('round038_pb_bound/check_bounds.py', 'checks.json')
        local = run('round038_pb_bound/highdim_small_checks.py', 'highdim_small_checks.json')
        audit = run('round038_conditional/read_saved_stats.py', 'HIGH_DIM_SAVED_STATS_AUDIT.json')
    out = {'scope': 'frozen hashes; exact small-source replay; saved DP/statistics and local checks only',
           'tail_MC_resampled': False, 'large_cache_built': False,
           'python': sys.version, 'numpy': np.__version__,
           'payload_hash_count': len(provenance['payload_SHA256']),
           'exact_conditional_cases': len(exact['cases']),
           'exact_conditional_replay_matches': True, 'layered_replay_matches': True,
           'saved_CI_checks': audit['checked_one_sided_CIs'],
           'max_CI_difference': audit['max_CI_absolute_recompute_error'],
           'max_ESS_difference': audit['max_ESS_absolute_recompute_error'],
           'BE_saved_DP_checks': be['BE_checks'], 'local_checks': local,
           'runs': runs, 'seconds': time.time()-started,
           'full_floating_point_certificate': False}
    if args.output:
        target = Path(args.output)
        if target.exists():
            raise RuntimeError('Refusing replacement of audit output')
        target.write_text(json.dumps(out, indent=2)+'\n')
    print(json.dumps(out, indent=2))


if __name__ == '__main__':
    main()
