"""Original workload group leader; native descendants inherit this group."""
import signal
import time
signal.signal(signal.SIGALRM, signal.SIG_DFL)
signal.setitimer(signal.ITIMER_REAL, 178)
import json
import os
import sys
from runtime_support import P, PRODUCTION, FIXTURE, require, production_gate


def main():
    mode, start, parent, gate, result, spec = sys.argv[1:]
    require(mode in ('fixture', 'production'), 'Invalid worker mode')
    policy = FIXTURE if mode == 'fixture' else PRODUCTION
    start, parent, gate, result = float(start), int(parent), int(gate), int(result)
    require(os.getpid() == os.getpgrp() and os.getppid() == parent, 'Worker must be original own group leader')
    def kill_own_group(signum, frame):
        os.killpg(os.getpgrp(), signal.SIGKILL)
    signal.signal(signal.SIGALRM, kill_own_group)
    signal.setitimer(signal.ITIMER_REAL, max(0.001, start + policy['worker_fallback']-time.monotonic()))
    def soft(signum, frame):
        raise TimeoutError('Fixed soft work cutoff; preserve original attempt, never retry')
    signal.signal(signal.SIGUSR1, soft)
    # A worker cannot import an oracle until the independently armed watchdog ACK.
    require(os.read(gate, 1) == b'G', 'No armed watchdog start permission')
    os.close(gate)
    os.set_blocking(result, False)
    def emit(value):
        os.write(result, (json.dumps(value, separators=(',', ':'))+'\n').encode())
    emit(dict(kind='WORKER_STARTED', pid=os.getpid(), pgid=os.getpgrp(), supervisor_pid=parent))
    try:
        if mode == 'fixture':
            from fixture_stubs import run
            data = run(spec, start, policy)
        else:
            operation = json.loads(spec)
            release_sha = production_gate(operation)
            guard = dict(mode=mode, started_monotonic=start, worker_pid=os.getpid(), worker_pgid=os.getpgrp(),
                         supervisor_pid=parent, runtime_release_sha256=release_sha, operation=operation)
            if operation[0] == 'native':
                import native_task
                sys.argv = ['native_task.py'] + operation[1:]
                native_task.main(guard)
            else:
                import control_task
                sys.argv = ['control_task.py'] + operation[1:]
                control_task.main(guard)
            data = dict(operation=operation, runtime_release_sha256=release_sha)
        emit(dict(kind='TASK_RESULT', logical_exit_code=0, detail=data))
    except BaseException as error:
        emit(dict(kind='TASK_RESULT', logical_exit_code=1, error=repr(error)))
    # Retain the original unreaped leader as the PG identity anchor until cleanup.
    signal.signal(signal.SIGUSR1, signal.SIG_IGN)
    while True:
        signal.pause()


if __name__ == '__main__':
    main()
