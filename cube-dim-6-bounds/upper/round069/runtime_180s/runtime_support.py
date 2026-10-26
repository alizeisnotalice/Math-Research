"""Fixed runtime policy and release boundary; no body/oracle imports."""
import hashlib
import json
import os
from pathlib import Path
import signal
import time

P = Path(__file__).resolve().parent
C068 = P.parents[1] / 'round068' / 'cover_batch001'
PRODUCTION = dict(soft=170.0, group=178.0, parent=179.0, watchdog=179.5,
                  supervisor=179.8, worker_fallback=178.0, reap=1.0)
FIXTURE = dict(soft=0.35, group=0.85, parent=1.25, watchdog=1.45,
               supervisor=1.6, worker_fallback=1.1, reap=0.12)
NAMES = ('normal_child', 'soft_timeout', 'cleanup_hang', 'descendant_group',
         'receipt_write_hang', 'gate_and_deadline_isolation')


def require(ok, message):
    if not ok:
        raise ValueError(message)


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read(path):
    return json.loads(Path(path).read_text())


def write_new(path, value):
    with Path(path).open('x') as f:
        json.dump(value, f, sort_keys=True, indent=2)
        f.write('\n')
        f.flush()
        os.fsync(f.fileno())


def arm_exit(seconds):
    """Default SIGALRM terminates even inside error handling/receipt writes."""
    signal.signal(signal.SIGALRM, signal.SIG_DFL)
    signal.setitimer(signal.ITIMER_REAL, max(0.001, seconds))


def validate_operation(args):
    require(args and args[0] in ('native', 'control'), 'Only native/control operation')
    if args[0] == 'native':
        import re
        require(len(args) == 2 and re.fullmatch(r'r[01]*\.a[0-4]\.n[0-9]{6}', args[1]),
                'Exactly one registered reservation token required')
    else:
        require(len(args) >= 2, 'Exactly one explicit control operation required')
        command = args[1]
        count = {'status': 2, 'authorize': 2, 'reserve': 2, 'attach': 4,
                 'verify': 3, 'import': 3, 'halt': 3}.get(command)
        require(count is not None and len(args) == count, 'No reset, timeout override or batch loop')
    return list(args)


def verify_pins():
    pins = read(P / 'input_manifest.json')
    for name, digest in pins['source_sha256'].items():
        require(sha(P / name) == digest, 'Runtime source changed: ' + name)
    for path, digest in pins['readonly_sha256'].items():
        require(sha(path) == digest, 'Frozen input changed: ' + path)
    return pins


def production_gate(args):
    validate_operation(args)
    verify_pins()
    release = read(P / 'runtime_release.json')  # deliberately absent in preparation
    require(release.get('schema') == 'C069-runtime-release-v1' and
            release.get('owner') == 'main' and release.get('production_authorized') is True,
            'No main runtime release')
    require(release.get('deadline_policy') == PRODUCTION and release.get('max_new_nodes') == 8 and
            release.get('max_native_reservations') == 8, 'Runtime/call policy mismatch')
    require(release.get('preparation_manifest_sha256') == sha(P / 'preparation_manifest.json'),
            'Runtime release not bound to this seal')
    manifest = read(P / 'preparation_manifest.json')
    for name, digest in manifest['artifact_sha256'].items():
        require(sha(P / name) == digest, 'Sealed runtime artifact changed: ' + name)
    acceptance = release['runtime_acceptance']
    require(sha(acceptance['path']) == acceptance['sha256'], 'Runtime acceptance hash mismatch')
    verdict = read(acceptance['path'])
    require(verdict.get('status') == 'PASS_RUNTIME_180S_ACCEPTED_BY_MAIN' and
            verdict.get('runtime_preparation_manifest_sha256') == release['preparation_manifest_sha256'],
            'Runtime audit/main acceptance missing')
    require(bool(verdict.get('independent_review_artifacts')), 'Missing independent runtime review')
    for item in verdict['independent_review_artifacts']:
        require(sha(item['path']) == item['sha256'], 'Independent runtime review changed')
    require(release.get('C068_main_release_sha256') == sha(C068 / 'main_release.json'),
            'C068 main release binding mismatch')
    # Importing the original context here checks its independent core/math release too.
    import sys
    sys.path.insert(0, str(C068))
    from batch_context import Context
    Context().verified_release(C068 / 'main_release.json', release['C068_main_release_sha256'])
    return sha(P / 'runtime_release.json')


def require_worker(guard):
    require(guard['mode'] == 'production' and guard['worker_pid'] == os.getpid() and
            guard['worker_pgid'] == os.getpgrp() == os.getpid() and
            guard['supervisor_pid'] == os.getppid(), 'Only the guarded original group leader may run tasks')
    require(time.monotonic() < guard['started_monotonic'] + PRODUCTION['soft'], 'Work deadline reached')
    require(production_gate(guard['operation']) == guard['runtime_release_sha256'], 'Release changed after launch')


def bounded_reap(child, deadline):
    try:
        code = child.wait(timeout=max(0.0001, deadline-time.monotonic()))
        return dict(reaped=True, returncode=code)
    except TimeoutError:
        return dict(reaped=False, returncode=None)
    except __import__('subprocess').TimeoutExpired:
        return dict(reaped=False, returncode=None)
