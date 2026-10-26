"""One explicit operation, never a batch loop; supervisor outside workload PG."""
import signal
import time
START = time.monotonic()
signal.signal(signal.SIGALRM, signal.SIG_DFL)
signal.setitimer(signal.ITIMER_REAL, 179.8)
import json
import os
from pathlib import Path
import select
import subprocess
import sys
sys.dont_write_bytecode = True
from runtime_support import (P, PRODUCTION, FIXTURE, NAMES, arm_exit, require,
                             write_new, production_gate, verify_pins, bounded_reap)


def main():
    require(len(sys.argv) >= 3 and sys.argv[1] in ('fixture', 'production'),
            'Usage: envelope.py fixture NAME | production native TOKEN | production control OP ...')
    mode, args = sys.argv[1], sys.argv[2:]
    if mode == 'fixture':
        require(len(args) == 1 and args[0] in NAMES, 'Only four fixed revision configurations')
        require(not (P/'runtime_release.json').exists(), 'Fixtures disabled once a runtime release exists')
        verify_pins()
        policy = FIXTURE
        spec = args[0]
        folder = P/'fixtures'/spec
        release_sha = None
    else:
        policy = PRODUCTION
        release_sha = production_gate(args)
        spec = json.dumps(args)
        folder = P/'operations'/f'{time.time_ns()}-{os.getpid()}'
    arm_exit(START + policy['supervisor'] - time.monotonic())
    folder.parent.mkdir(exist_ok=True)
    folder.mkdir()  # no overwrite or restart of a fixture/operation
    write_new(folder/'start.json', dict(schema='C069-runtime-start-v1', mode=mode, args=args,
              supervisor_pid=os.getpid(), supervisor_pgid=os.getpgrp(), started_monotonic=START,
              deadline_policy=policy, runtime_release_sha256=release_sha,
              native_calls=0 if mode=='fixture' else None, automatic_retries=0))
    gate_r, gate_w = os.pipe()
    result_r, result_w = os.pipe()
    control_r, control_w = os.pipe()
    event_r, event_w = os.pipe()
    worker = watchdog = None
    events = []
    buffers = {result_r:b'', event_r:b''}
    soft_sent = False
    group_killed = False
    group_ack = False
    cleanup_error = None
    original_worker_observation = None
    logical = None
    timed_out = False
    reaps = {}
    deadline = START + policy['group'] + 0.04
    stdout = (folder/'worker.stdout').open('xb')
    stderr = (folder/'worker.stderr').open('xb')
    dogout = (folder/'watchdog.stdout').open('xb')
    dogerr = (folder/'watchdog.stderr').open('xb')
    def drain(wait):
        nonlocal logical, group_ack
        ready, _, _ = select.select([result_r, event_r], [], [], wait)
        for fd in ready:
            raw = os.read(fd, 16384)
            buffers[fd] += raw
            while b'\n' in buffers[fd]:
                line, buffers[fd] = buffers[fd].split(b'\n', 1)
                value = json.loads(line)
                events.append(value)
                if fd == result_r and value['kind'] == 'TASK_RESULT':
                    logical = value
                if fd == event_r and value.get('pid') == watchdog.pid and value['kind'] in (
                        'GROUP_KILL_ACKNOWLEDGED', 'GROUP_SIGKILL_SENT', 'GROUP_NOT_FOUND'):
                    group_ack = True
    try:
        worker = subprocess.Popen([sys.executable, '-B', str(P/'worker.py'), mode, str(START), str(os.getpid()),
                                   str(gate_r), str(result_w), spec], pass_fds=(gate_r,result_w),
                                  start_new_session=True, stdout=stdout, stderr=stderr)
        # Native tasks never create another session; the original leader is retained.
        watchdog = subprocess.Popen([sys.executable, '-B', str(P/'watchdog.py'), mode, str(START), str(os.getpid()),
                                     str(worker.pid), str(control_r), str(event_w)],
                                    pass_fds=(control_r,event_w), stdout=dogout, stderr=dogerr)
        for fd in (gate_r, result_w, control_r, event_w):
            os.close(fd)
        os.set_blocking(result_r, False)
        os.set_blocking(event_r, False)
        os.set_blocking(control_w, False)
        write_new(folder/'handles.json', dict(supervisor_pid=os.getpid(), worker_pid=worker.pid,
                  worker_pgid=worker.pid, watchdog_pid=watchdog.pid, inherited_native_pgid=worker.pid))
        while not any(e['kind']=='WATCHDOG_ARMED' for e in events):
            require(time.monotonic() < START+policy['soft'], 'Watchdog did not arm before work deadline')
            drain(0.005)
        os.write(gate_w, b'G')
        os.close(gate_w)
        while logical is None and not group_ack and time.monotonic() < deadline:
            if not soft_sent and time.monotonic() >= START+policy['soft']:
                os.kill(worker.pid, signal.SIGUSR1)
                soft_sent = True
                events.append(dict(kind='SOFT_WORK_CUTOFF_SENT', elapsed=time.monotonic()-START))
            drain(0.005)
        timed_out = logical is None
    except BaseException as error:
        events.append(dict(kind='SUPERVISOR_ERROR', error=repr(error), elapsed=time.monotonic()-START))
        timed_out = True
    # No total-deadline cancellation: cleanup and all writes remain under both fuses.
    cleanup_end = min(time.monotonic()+policy['reap'], START+policy['parent']-0.04)
    if worker is not None:
        # The original leader has not been wait/poll/reaped and cannot be PID-reused.
        # A recorded independent termination ACK makes a second kill unnecessary.
        group_absent = False
        if not group_ack:
            try:
                if mode=='fixture' and spec=='receipt_error_gate_unacked_eperm':
                    raise PermissionError(1, 'Registered fixture-only unacknowledged EPERM')
                os.killpg(worker.pid, signal.SIGKILL)
                group_killed = True
            except ProcessLookupError:
                group_absent = True
            except OSError as error:
                cleanup_error = dict(error=repr(error), errno=error.errno, worker_pid=worker.pid,
                                     worker_pgid=worker.pid, identity_anchor_unreaped=True,
                                     termination_confirmed=False, elapsed=time.monotonic()-START)
                events.append(dict(kind='GROUP_KILL_ERROR_UNRESOLVED', **cleanup_error))
                timed_out = True
        if watchdog is not None and (group_killed or group_absent):
            try:
                os.write(control_w, b'K')
            except (BrokenPipeError, BlockingIOError) as error:
                events.append(dict(kind='WATCHDOG_ACK_REQUEST_FAILED', error=repr(error)))
        # Wait boundedly for dog acknowledgment before releasing the PG identity anchor.
        ack_end = min(time.monotonic()+policy['reap']/2, cleanup_end)
        while not group_ack and time.monotonic() < ack_end:
            drain(0.002)
        if cleanup_error is not None:
            # Observe only the original unreaped PID, never a name search or new kill.
            try:
                observation = subprocess.run(['/bin/ps','-p',str(worker.pid),'-o','pid=,ppid=,pgid=,stat='],
                                             capture_output=True, text=True,
                                             timeout=max(0.001,min(0.1,cleanup_end-time.monotonic())))
                original_worker_observation = dict(pid=worker.pid, ps_exit_code=observation.returncode,
                    stdout=observation.stdout, stderr=observation.stderr, identity_anchor_unreaped=True)
            except (OSError,subprocess.TimeoutExpired) as error:
                original_worker_observation = dict(pid=worker.pid,error=repr(error),state='UNKNOWN')
        reap_end = cleanup_end
        if group_ack or watchdog is None:
            reaps['worker'] = bounded_reap(worker, reap_end)
        else:
            reaps['worker'] = dict(reaped=False, returncode=None, reason='No watchdog ACK; PG identity not released')
    receipt = dict(schema='C069-runtime-receipt-v1', mode=mode, args=args,
                   status='TASK_RESULT_OBSERVED' if logical is not None and not timed_out else 'STOP_NO_RETRY',
                   task_result=logical, events=events, group_sigkill_sent_by_supervisor=group_killed,
                   group_kill_acknowledged=group_ack, reaps=reaps, deadline_policy=policy,
                   cleanup_error=cleanup_error, original_worker_observation=original_worker_observation,
                   supervisor_pid=os.getpid(), worker_pid=worker.pid if worker else None,
                   watchdog_pid=watchdog.pid if watchdog else None, elapsed_before_receipt=time.monotonic()-START,
                   all_descendants_reaped_claim=False, automatic_retries=0,
                   native_calls=0 if mode=='fixture' else None)
    receipt_written = False
    try:
        write_new(folder/'before_receipt.json', receipt)
        if mode=='production' and args[0]=='native' and (
                timed_out or logical is None or logical['logical_exit_code']!=0):
            from runtime_support import record_native_stop
            record_native_stop(args[1], folder/'before_receipt.json', receipt)
        if mode=='fixture' and spec=='receipt_retained_reader_hang':
            from fixture_stubs import block_full_pipe
            block_full_pipe(folder, 'receipt', START)
        if mode=='fixture' and spec=='receipt_error_gate_unacked_eperm':
            (folder/'receipt.json').mkdir()  # registered EISDIR receipt-write injection
        write_new(folder/'receipt.json', receipt)
        receipt_written = True
    except BaseException as error:
        timed_out = True
        write_new(folder/'receipt_failure.json', dict(status='STOP_RECEIPT_WRITE_ERROR_NO_RETRY',
                  error=repr(error), elapsed=time.monotonic()-START, group_ack=group_ack,
                  cleanup_error=cleanup_error, reaps=reaps, final_receipt_written=False,
                  all_descendants_reaped_claim=False, automatic_retries=0))
    if watchdog is not None:
        if group_ack:
            try:
                os.write(control_w, b'F')
            except (BrokenPipeError, BlockingIOError):
                pass
        else:
            # Do not send a false termination ACK/disarm. EOF tells the original
            # watchdog to attempt its own group cleanup; held PID remains unreaped.
            os.close(control_w)
        reaps['watchdog'] = bounded_reap(watchdog, min(cleanup_end, START+policy['supervisor']-0.02))
        if reaps['watchdog']['reaped'] and worker is not None and not reaps['worker']['reaped']:
            reaps['worker'] = bounded_reap(worker, cleanup_end)
    write_new(folder/'completion.json', dict(reaps=reaps, elapsed=time.monotonic()-START,
              final_receipt_written=receipt_written, cleanup_error=cleanup_error,
              all_descendants_reaped_claim=False, automatic_retries=0))
    print(json.dumps(dict(folder=str(folder), receipt=str(folder/'receipt.json'), task_result=logical,
                          elapsed=time.monotonic()-START)))
    return 0 if logical is not None and logical['logical_exit_code']==0 and not timed_out else 1


if __name__ == '__main__':
    sys.exit(main())
