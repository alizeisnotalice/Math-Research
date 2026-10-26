"""One-shot, six registered runtime-only configurations, no automatic retry."""
import signal
import time
START = time.monotonic()
signal.signal(signal.SIGALRM, signal.SIG_DFL)
signal.alarm(170)
import ast
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import shutil
import stat
import subprocess
import sys
import traceback
sys.dont_write_bytecode = True

P = Path(__file__).resolve().parent
LOCK = json.loads((P / 'input_manifest.json').read_text())
AUTHOR = Path(LOCK['author_root'])
REG = json.loads((P / 'registration.json').read_text())
ROWS = []
OWNED = []


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def save(path, data):
    with path.open('x') as stream:
        json.dump(data, stream, indent=2, sort_keys=True)
        stream.write('\n')


def check_pins():
    assert sha(P / 'registration.json') == LOCK['registration_sha256']
    assert sha(AUTHOR / 'preparation_manifest.json') == LOCK['author_seal_sha256']
    for name, digest in LOCK['author_artifact_sha256'].items():
        assert sha(AUTHOR / name) == digest, name
    for name, digest in LOCK['readonly_sha256'].items():
        assert sha(Path(name)) == digest, name
    assert not (AUTHOR / 'runtime_release.json').exists(), 'Production release appeared'


def ps(pid):
    result = subprocess.run(['/bin/ps', '-o', 'pid=,ppid=,pgid=,stat=', '-p', str(pid)],
                            capture_output=True, text=True, timeout=0.5)
    text = result.stdout.strip()
    fields = text.split()
    return dict(pid=pid, ps_exit=result.returncode, raw=text,
                live=bool(fields) and not fields[-1].startswith('Z'))


def stop_direct(proc):
    if proc.poll() is None:
        proc.kill()
    try:
        return proc.wait(timeout=0.4)
    except subprocess.TimeoutExpired:
        return None


def wait_file(path, deadline):
    while not path.exists():
        assert time.monotonic() < deadline, 'Missing bounded fixture marker: ' + str(path)
        time.sleep(0.003)
    return json.loads(path.read_text())


def make_copy(case):
    target = P / '_scratch' / case['id']
    target.mkdir(parents=True, exist_ok=False)
    pins = json.loads((AUTHOR / 'input_manifest.json').read_text())
    for name, digest in LOCK['author_source_sha256'].items():
        shutil.copyfile(AUTHOR / name, target / name)
        assert sha(target / name) == digest
    # Only the fixture module and its local pin are deliberately replaced.
    shutil.copyfile(P / 'independent_stub.py', target / 'fixture_stubs.py')
    pins['source_sha256']['fixture_stubs.py'] = sha(target / 'fixture_stubs.py')
    save(target / 'input_manifest.json', pins)
    (target / 'independent_case.txt').write_text(case['id'] + '\n')
    return target, {name: sha(target / name) for name in LOCK['author_source_sha256']}


def collect(folder):
    files, specials = {}, {}
    if not folder.exists():
        return dict(files=files, special_files=specials)
    for item in sorted(folder.iterdir()):
        mode = item.lstat().st_mode
        if stat.S_ISREG(mode):
            raw = item.read_bytes()
            text = raw.decode('utf-8')
            files[item.name] = dict(sha256=hashlib.sha256(raw).hexdigest(), bytes=len(raw),
                content=json.loads(text) if item.suffix == '.json' else text)
        else:
            specials[item.name] = 'FIFO' if stat.S_ISFIFO(mode) else 'DIRECTORY' if stat.S_ISDIR(mode) else str(mode)
    return dict(files=files, special_files=specials)


def observed_pids(folder):
    values = []
    if (folder / 'handles.json').exists():
        data = json.loads((folder / 'handles.json').read_text())
        values.extend(data[k] for k in ('supervisor_pid', 'worker_pid', 'watchdog_pid'))
    for name in ('child.json', 'grandchild.json'):
        if (folder / name).exists():
            values.append(json.loads((folder / name).read_text())['pid'])
    return sorted(set(values))


def owned_emergency_cleanup(proc, folder):
    # Only original numeric handles belonging to this single newly started test.
    # Worker leader is still owned by the supervisor while this is needed.
    if proc.poll() is None and (folder / 'handles.json').exists():
        handles = json.loads((folder / 'handles.json').read_text())
        try:
            os.killpg(handles['worker_pgid'], signal.SIGKILL)
        except ProcessLookupError:
            pass
    if proc.poll() is None:
        try:
            os.killpg(proc.pid, signal.SIGKILL)
        except ProcessLookupError:
            pass
    stop_direct(proc)


def run_fault(case):
    target, hashes = make_copy(case)
    folder = target / 'fixtures' / case['entry_name']
    before = time.monotonic()
    sentinel = None
    out = (target / 'supervisor.stdout').open('xb')
    err = (target / 'supervisor.stderr').open('xb')
    proc = subprocess.Popen([sys.executable, '-B', str(target / 'envelope.py'),
                             'fixture', case['entry_name']],
                            stdout=out, stderr=err, start_new_session=True)
    OWNED.append(proc.pid)
    observation = dict(id=case['id'], started_utc=datetime.now(timezone.utc).isoformat(),
                       source_copy_sha256=hashes, original_supervisor_pid=proc.pid,
                       reviewer_emergency_cleanup=False, reviewer_injected_signal=None)
    try:
        if case['id'].startswith('R2_'):
            sentinel = subprocess.Popen([sys.executable, '-B', '-c',
                'import signal,time; signal.signal(signal.SIGALRM,signal.SIG_DFL); signal.alarm(4); time.sleep(3)'],
                stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, start_new_session=True)
            OWNED.append(sentinel.pid)
            observation['sentinel_pid'] = sentinel.pid
        if case['id'].startswith('R5_'):
            handles = wait_file(folder / 'handles.json', before + 2)
            entered = wait_file(folder / 'cleanup_block_entered.json', before + 2)
            assert entered['retained_reader'] and entered['filled_bytes'] > 0
            os.kill(handles['watchdog_pid'], signal.SIGKILL)
            observation['reviewer_injected_signal'] = dict(pid=handles['watchdog_pid'],
                signal='SIGKILL', cause='registered owned watchdog fault', elapsed=time.monotonic()-before)
        code = proc.wait(timeout=max(0.01, before + 4 - time.monotonic()))
        observation.update(supervisor_exit_code=code, driver_wall_seconds=time.monotonic()-before)
        if sentinel:
            observation['sentinel_before_reviewer_cleanup'] = ps(sentinel.pid)
            assert observation['sentinel_before_reviewer_cleanup']['live'], 'Unrelated sentinel was killed'
        # Bounded observation window permits init to collect dead grandchildren.
        states = [ps(pid) for pid in observed_pids(folder)]
        deadline = time.monotonic() + 0.35
        while any(item['live'] for item in states) and time.monotonic() < deadline:
            time.sleep(0.01)
            states = [ps(pid) for pid in observed_pids(folder)]
        observation['owned_process_observations'] = states
        evidence = collect(folder)
        observation.update(evidence)
        observation['supervisor_stdout'] = (target / 'supervisor.stdout').read_text()
        observation['supervisor_stderr'] = (target / 'supervisor.stderr').read_text()
        assert not any(item['live'] for item in states), 'Owned process still live'
        assert observation['driver_wall_seconds'] < 2.1, 'Fixture exceeded bounded expected duration'
        files = evidence['files']
        before_receipt = files['before_receipt.json']['content']
        events = before_receipt['events']
        kinds = [item['kind'] for item in events]
        if case['id'].startswith(('R1_', 'R5_')):
            marker = files['cleanup_block_entered.json']['content']
            assert marker['retained_reader'] and marker['filled_bytes'] > 0
            assert 'cleanup_returned.json' not in files
            assert before_receipt['task_result'] is None and code == 1
            assert before_receipt['status'] == 'STOP_NO_RETRY'
            if case['id'].startswith('R1_'):
                assert 'GROUP_SIGKILL_SENT' in kinds, 'Independent watchdog did not fire'
                kill = next(item for item in events if item['kind'] == 'GROUP_SIGKILL_SENT')
                assert 0.85 <= kill['elapsed'] < 1.1, kill
                assert before_receipt['group_kill_acknowledged'] is True
            else:
                assert 'GROUP_SIGKILL_SENT' not in kinds
                assert before_receipt['group_sigkill_sent_by_supervisor'] is True
                assert before_receipt['group_kill_acknowledged'] is False
                assert before_receipt['reaps']['worker']['reaped'] is False
        elif case['id'].startswith('R2_'):
            handles = files['handles.json']['content']
            for name in ('child.json', 'grandchild.json'):
                assert files[name]['content']['pgid'] == handles['worker_pgid']
            assert files['early_return.json']['content']['child_poll'] is None
            assert observation['sentinel_pid'] != handles['worker_pgid']
            assert code == 0 and before_receipt['task_result']['logical_exit_code'] == 0
            assert before_receipt['group_sigkill_sent_by_supervisor'] is True
        elif case['id'].startswith('R3_'):
            assert evidence['special_files'].get('blocked_receipt.fifo') == 'FIFO'
            assert 'receipt.json' not in files and 'completion.json' not in files
            assert code == -signal.SIGKILL
            start = files['start.json']['content']['started_monotonic']
            # Exit SIGKILL before the supervisor's SIGALRM fallback; the source's
            # only parent SIGKILL comes from the independent watchdog at 1.25.
            inferred = time.monotonic() - start
            observation['post_observation_elapsed_from_supervisor_entry'] = inferred
            assert 1.20 < observation['driver_wall_seconds'] < 1.6
        elif case['id'].startswith('R4_'):
            assert evidence['special_files'].get('receipt.json') == 'DIRECTORY'
            assert 'completion.json' not in files and code != 0
            assert any(name in observation['supervisor_stderr'] for name in ('FileExistsError', 'IsADirectoryError'))
            assert before_receipt['group_sigkill_sent_by_supervisor'] is True
        observation['status'] = 'PASS'
    except BaseException as error:
        observation['status'] = 'FAIL_STOP'
        observation['error'] = repr(error)
        observation['traceback'] = traceback.format_exc()
        observation['reviewer_emergency_cleanup'] = True
        owned_emergency_cleanup(proc, folder)
        observation.update(collect(folder))
        observation['supervisor_stdout'] = (target / 'supervisor.stdout').read_text()
        observation['supervisor_stderr'] = (target / 'supervisor.stderr').read_text()
        raise
    finally:
        if sentinel:
            observation['sentinel_reviewer_cleanup_exit'] = stop_direct(sentinel)
            observation['sentinel_after_reviewer_cleanup'] = ps(sentinel.pid)
        out.close()
        err.close()
        save(P / (case['id'] + '.json'), observation)
        ROWS.append(dict(id=case['id'], status=observation['status'],
            supervisor_exit_code=observation.get('supervisor_exit_code'),
            driver_wall_seconds=observation.get('driver_wall_seconds'),
            observation_sha256=sha(P / (case['id'] + '.json'))))


def readonly_snapshot():
    root = AUTHOR.parents[1] / 'round068' / 'cover_batch001'
    paths = [root / name for name in ('state', 'attempts', 'operation_receipts')]
    paths.append(AUTHOR / 'operations')
    return {str(base): {str(x.relative_to(base)): sha(x) for x in base.rglob('*') if x.is_file()}
            if base.exists() else None for base in paths}


def run_gate(case):
    started = time.monotonic()
    before = readonly_snapshot()
    probes = []
    commands = [
        ['envelope.py', 'production', 'native', 'r.a0.n000001'],
        ['envelope.py', 'production', 'control', 'status'],
        ['native_task.py'], ['control_task.py']]
    for args in commands:
        result = subprocess.run([sys.executable, '-B', str(AUTHOR / args[0]), *args[1:]],
                                capture_output=True, text=True, timeout=3)
        assert result.returncode != 0
        assert ('runtime_release.json' in result.stderr if args[0] == 'envelope.py'
                else 'standalone execution is forbidden' in result.stderr)
        probes.append(dict(args=args, exit_code=result.returncode, stdout=result.stdout, stderr=result.stderr))
    assert readonly_snapshot() == before, 'Rejected gate wrote operation/state/attempt artifacts'
    files = ('envelope.py', 'watchdog.py', 'worker.py', 'runtime_support.py', 'native_task.py', 'control_task.py')
    parsed = {name: ast.parse((AUTHOR / name).read_text()) for name in files}
    waits, alarms, pops = [], [], []
    for name, tree in parsed.items():
        for node in ast.walk(tree):
            if not isinstance(node, ast.Call) or not isinstance(node.func, ast.Attribute):
                continue
            if node.func.attr == 'wait':
                assert any(k.arg == 'timeout' for k in node.keywords), (name, node.lineno)
                waits.append(dict(file=name, line=node.lineno, bounded=True))
            if node.func.attr in ('alarm', 'setitimer'):
                assert not any(isinstance(arg, ast.Constant) and arg.value == 0 for arg in node.args)
                alarms.append(dict(file=name, line=node.lineno))
            if node.func.attr == 'Popen':
                pops.append(dict(file=name, line=node.lineno,
                    start_new_session=next((ast.literal_eval(k.value) for k in node.keywords
                                            if k.arg == 'start_new_session'), False)))
    assert next(p['start_new_session'] for p in pops if p['file'] == 'native_task.py') is False
    envelope_pops = [p['start_new_session'] for p in pops if p['file'] == 'envelope.py']
    assert envelope_pops == [True, False]
    # Parse fixed policy directly, without executing any copied production task.
    policies = {}
    for node in parsed['runtime_support.py'].body:
        if isinstance(node, ast.Assign) and isinstance(node.value, ast.Call):
            name = node.targets[0].id if isinstance(node.targets[0], ast.Name) else None
            if name in ('PRODUCTION', 'FIXTURE'):
                policies[name] = {k.arg: ast.literal_eval(k.value) for k in node.value.keywords}
    assert policies['PRODUCTION'] == dict(soft=170., group=178., parent=179., watchdog=179.5,
                                       supervisor=179.8, worker_fallback=178., reap=1.)
    observation = dict(id=case['id'], status='PASS', probes=probes,
        unchanged_state_attempt_operation_snapshot=before, fixed_policies=policies,
        explicit_wait_timeout_calls=waits, alarm_calls_without_cancellation=alarms,
        process_launch_group_flags=pops, elapsed_seconds=time.monotonic()-started,
        native_calls=0, body_kernel_calls=0, Arb_calls=0)
    save(P / (case['id'] + '.json'), observation)
    ROWS.append(dict(id=case['id'], status='PASS', elapsed_seconds=observation['elapsed_seconds'],
                    observation_sha256=sha(P / (case['id'] + '.json'))))


def main():
    check_pins()
    save(P / 'run_started.json', dict(pid=os.getpid(), started_utc=datetime.now(timezone.utc).isoformat(),
        registration_sha256=sha(P / 'registration.json'), harness_sha256=sha(Path(__file__)),
        stub_sha256=sha(P / 'independent_stub.py'), python=sys.executable, automatic_retries=0))
    code = 0
    failure = None
    try:
        for case in REG['configurations']:
            if case['entry_name'] is None:
                run_gate(case)
            else:
                run_fault(case)
            print(case['id'], 'PASS', flush=True)
        check_pins()
    except BaseException as error:
        failure = dict(error=repr(error), traceback=traceback.format_exc())
        code = 1
    result = dict(schema='A069-independent-runtime-review-results-v1',
        status='PASS_SIX_REGISTERED_RUNTIME_CONFIGURATIONS' if code == 0 else 'FAIL_STOP_NO_RETRY',
        registered_configurations=6, started_configurations=len(ROWS),
        passed_configurations=sum(row['status'] == 'PASS' for row in ROWS),
        results=ROWS, failure=failure, elapsed_seconds=time.monotonic()-START,
        body_kernel_calls=0, native_calls=0, Arb_calls=0, old_normal_soft_reruns=0,
        production_gate='CLOSED_NO_RELEASE_CREATED', automatic_retries=0,
        author_source_and_readonly_pins_unchanged=True if code == 0 else None,
        all_descendants_reaped_claim=False, actual_body_native_timeout_path_tested=False)
    save(P / 'results.json', result)
    print(json.dumps(result, sort_keys=True), flush=True)
    return code


if __name__ == '__main__':
    sys.exit(main())
