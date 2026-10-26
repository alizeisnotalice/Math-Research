"""One already-reserved C068 native attempt, hard 180 s, no state reset/retry.

Usage after reviewed main release: python -B run_native_once.py TOKEN
The unique main scheduler must persist RESERVE_NEXT_ONCE before this command.
"""
import sys
sys.dont_write_bytecode = True
import json
import os
from pathlib import Path
import signal
import subprocess
import time
from datetime import datetime, timezone
import core_v2 as V
from batch_context import Context, P, read, sha, write_new, require
import math_contract as M


def utc():
    return datetime.now(timezone.utc).isoformat()


def main():
    require(len(sys.argv) == 2, 'Usage: run_native_once.py TOKEN')
    tic = time.monotonic()
    def timeout(signum, frame):
        raise TimeoutError('180-second whole-program wall limit')
    signal.signal(signal.SIGALRM, timeout)
    signal.setitimer(signal.ITIMER_REAL, 180)
    context = Context()
    latest = V.recover(P / 'state', context)
    state, token = latest['state'], sys.argv[1]
    require(state['environment'] == 'MAIN_RELEASED' and state['release'] is not None, 'No reviewed main production release')
    context.verified_release(P / 'main_release.json', state['release']['release_sha256'])
    require(state['halted'] is None and state['inflight'] is not None and state['inflight']['token'] == token,
            'No matching persisted original reservation')
    node = state['nodes'][state['inflight']['node_id']]
    attempt = node['attempts'][-1]
    require(attempt['state'] == 'RESERVED_BEFORE_ANY_LAUNCH' and attempt['execution_handle'] is None,
            'An observed/resolved attempt cannot launch again')
    expected = attempt['input']
    (P / 'attempts').mkdir(exist_ok=True)
    folder = context.attempt_directory(token)
    folder.mkdir()  # exclusive launch intent: existing attempt is never replayed
    write_new(folder/'reserved_input.json', dict(token=token, input=expected, input_sha256=attempt['input_sha256'],
                                               reserved_snapshot_sha256=latest['sha256'], owned_wrapper_pid=os.getpid(),
                                               timeout_seconds=180, synthetic_fixture=False))
    try:
        command = M.native_command(context, expected)
    except BaseException as error:
        signal.setitimer(signal.ITIMER_REAL, 0)
        write_new(folder/'failure.json', dict(token=token, input_sha256=attempt['input_sha256'],
                  reason='VALIDATION_ERROR_STOP', error=repr(error), native_launch_observed=False,
                  retained_artifacts=[dict(path=str(folder/'reserved_input.json'),sha256=sha(folder/'reserved_input.json'))],
                  synthetic_fixture=False, automatic_retries=0, timeout_seconds=180, wall_seconds=time.monotonic()-tic))
        raise
    intent = dict(schema='C068-native-launch-intent-v1', token=token, input=expected,
                  input_sha256=attempt['input_sha256'], reserved_snapshot_sha256=latest['sha256'],
                  release_sha256=state['release']['release_sha256'], command=command,
                  source_sha256=context.external_sha(context.paths['source']),
                  binary_sha256=context.external_sha(context.paths['binary']),
                  fan_sha256=context.external_sha(context.paths['fan']), runtime_sha256=context.runtime_pins,
                  owned_wrapper_pid=os.getpid(), started_utc=utc(), timeout_seconds=180, synthetic_fixture=False)
    write_new(folder / 'intent.json', intent)
    native_path, stderr_path = folder / 'native.json', folder / 'stderr.txt'
    child = None
    receipt = None
    try:
        with native_path.open('xb') as out, stderr_path.open('xb') as err:
            child = subprocess.Popen(command, cwd=P, stdout=out, stderr=err, start_new_session=True)
            process_record = dict(token=token, input_sha256=attempt['input_sha256'],
                                  owned_wrapper_pid=os.getpid(), owned_native_pid=child.pid,
                                  owned_process_group_id=child.pid, launched_utc=utc(),
                                  timeout_seconds=180, status='RUNNING_ORIGINAL_NATIVE')
            write_new(folder / 'original_process.json', process_record)
            code = child.wait(timeout=max(0.001, 180 - (time.monotonic() - tic)))
        context.verify_pins()
        receipt = dict(schema='C068-native-receipt-v1', token=token, input_sha256=attempt['input_sha256'],
                       command=command, **{key:expected[key] for key in ('model','model_definition','window','center')},
                       source_sha256=intent['source_sha256'], binary_sha256=intent['binary_sha256'],
                       fan_sha256=intent['fan_sha256'], runtime_sha256=context.runtime_pins,
                       intent_sha256=sha(folder / 'intent.json'), release_sha256=state['release']['release_sha256'],
                       original_process_sha256=sha(folder / 'original_process.json'),
                       native_exit_code=code, stderr_empty=stderr_path.stat().st_size == 0,
                       native_sha256=sha(native_path), stderr_sha256=sha(stderr_path),
                       owned_wrapper_pid=os.getpid(), owned_native_pid=child.pid,
                       started_utc=intent['started_utc'], ended_utc=utc(), wall_seconds=time.monotonic()-tic,
                       timeout_seconds=180, synthetic_fixture=False, automatic_retries=0)
        write_new(folder / 'native_receipt.json', receipt)
        require(code == 0 and receipt['stderr_empty'], 'Native failed: preserve receipt and halt without retry')
        print(json.dumps(dict(token=token, native_pid=child.pid, exit_code=code,
                              native_receipt=str(folder / 'native_receipt.json'), next='MAIN_ATTACH_ORIGINAL_HANDLE_THEN_VERIFY')))
    except BaseException as error:
        signal.setitimer(signal.ITIMER_REAL, 0)
        if child is not None and child.poll() is None:
            try:
                os.killpg(child.pid, signal.SIGKILL)
            except ProcessLookupError:
                pass
            child.wait()
        retained = [dict(path=str(path), sha256=sha(path)) for path in
                    (folder/'reserved_input.json', folder/'intent.json', folder/'original_process.json', native_path, stderr_path, folder/'native_receipt.json') if path.exists()]
        failure = dict(token=token, input_sha256=attempt['input_sha256'],
                       reason='PROGRAM_TIMEOUT_STOP' if isinstance(error,(TimeoutError,subprocess.TimeoutExpired)) else 'NATIVE_ERROR_STOP',
                       error=repr(error), retained_artifacts=retained, native_launch_observed=child is not None,
                       native_exit_code=child.returncode if child is not None else None,
                       only_owned_group_reaped=child is not None, ended_utc=utc(), wall_seconds=time.monotonic()-tic,
                       timeout_seconds=180, synthetic_fixture=False, automatic_retries=0)
        write_new(folder / 'failure.json', failure)
        raise
    finally:
        signal.setitimer(signal.ITIMER_REAL, 0)


if __name__ == '__main__':
    main()
