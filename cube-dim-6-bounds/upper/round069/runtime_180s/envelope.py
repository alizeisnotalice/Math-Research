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
        require(len(args) == 1 and args[0] in NAMES, 'Only six fixed short stub configurations')
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
                if value['kind'] == 'TASK_RESULT':
                    logical = value
                if value['kind'] in ('GROUP_KILL_ACKNOWLEDGED', 'GROUP_SIGKILL_SENT', 'GROUP_NOT_FOUND'):
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
    if worker is not None:
        # The original leader has not been wait/poll/reaped and cannot be PID-reused.
        try:
            os.killpg(worker.pid, signal.SIGKILL)
            group_killed = True
        except ProcessLookupError:
            pass
        if watchdog is not None:
            try:
                os.write(control_w, b'K')
            except (BrokenPipeError, BlockingIOError):
                pass
        # Wait boundedly for dog acknowledgment before releasing the PG identity anchor.
        ack_end = min(time.monotonic()+policy['reap']/2, START+policy['parent']-0.05)
        while not group_ack and time.monotonic() < ack_end:
            drain(0.002)
        reap_end = min(time.monotonic()+policy['reap']/2, START+policy['parent']-0.04)
        if group_ack or watchdog is None:
            reaps['worker'] = bounded_reap(worker, reap_end)
        else:
            reaps['worker'] = dict(reaped=False, returncode=None, reason='No watchdog ACK; PG identity not released')
    receipt = dict(schema='C069-runtime-receipt-v1', mode=mode, args=args,
                   status='TASK_RESULT_OBSERVED' if logical is not None and not timed_out else 'STOP_NO_RETRY',
                   task_result=logical, events=events, group_sigkill_sent_by_supervisor=group_killed,
                   group_kill_acknowledged=group_ack, reaps=reaps, deadline_policy=policy,
                   supervisor_pid=os.getpid(), worker_pid=worker.pid if worker else None,
                   watchdog_pid=watchdog.pid if watchdog else None, elapsed_before_receipt=time.monotonic()-START,
                   all_descendants_reaped_claim=False, automatic_retries=0,
                   native_calls=0 if mode=='fixture' else None)
    write_new(folder/'before_receipt.json', receipt)
    if mode=='fixture' and spec=='receipt_write_hang':
        os.mkfifo(folder/'blocked_receipt.fifo')
        # Real blocking open: independent dog must kill this supervisor during the write phase.
        with (folder/'blocked_receipt.fifo').open('w') as stream:
            stream.write('must never complete without a reader')
    write_new(folder/'receipt.json', receipt)
    if watchdog is not None:
        try:
            os.write(control_w, b'F')
        except (BrokenPipeError, BlockingIOError):
            pass
        reaps['watchdog'] = bounded_reap(watchdog, min(time.monotonic()+policy['reap']/2,
                                                     START+policy['supervisor']-0.02))
    write_new(folder/'completion.json', dict(reaps=reaps, elapsed=time.monotonic()-START,
              final_receipt_written=True, all_descendants_reaped_claim=False, automatic_retries=0))
    print(json.dumps(dict(folder=str(folder), receipt=str(folder/'receipt.json'), task_result=logical,
                          elapsed=time.monotonic()-START)))
    return 0 if logical is not None and logical['logical_exit_code']==0 and not timed_out else 1


if __name__ == '__main__':
    sys.exit(main())
