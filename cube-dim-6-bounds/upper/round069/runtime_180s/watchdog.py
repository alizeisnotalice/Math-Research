"""Independent, fixed-deadline process-group and supervisor fuse.

No potentially blocking disk write precedes either deadline action. The pipe
comes from the original supervisor; no PID is discovered by name or scan.
"""
import signal
import time
START = time.monotonic()
signal.signal(signal.SIGALRM, signal.SIG_DFL)
signal.setitimer(signal.ITIMER_REAL, 179.5)
import json
import os
import select
import sys
from runtime_support import PRODUCTION, FIXTURE, arm_exit, require


def main():
    mode, start, parent, pgid, control, events = sys.argv[1:]
    require(mode in ('fixture', 'production'), 'Invalid watchdog mode')
    policy = FIXTURE if mode == 'fixture' else PRODUCTION
    start, parent, pgid, control, events = float(start), int(parent), int(pgid), int(control), int(events)
    arm_exit(start + policy['watchdog'] - time.monotonic())
    require(parent == os.getppid() and pgid > 1 and pgid != os.getpgrp(), 'Not an owned external worker group')
    os.set_blocking(control, False)
    os.set_blocking(events, False)
    def emit(kind, **extra):
        try:
            os.write(events, (json.dumps(dict(kind=kind, pid=os.getpid(), elapsed=time.monotonic()-start, **extra))+'\n').encode())
        except (BlockingIOError, BrokenPipeError):
            pass
    emit('WATCHDOG_ARMED', supervisor_pid=parent, worker_pgid=pgid)
    group_done = False
    parent_done = False
    while True:
        now = time.monotonic()
        if not group_done and (now >= start + policy['group'] or os.getppid() != parent):
            try:
                os.killpg(pgid, signal.SIGKILL)
                emit('GROUP_SIGKILL_SENT', pgid=pgid)
            except ProcessLookupError:
                emit('GROUP_NOT_FOUND', pgid=pgid)
            group_done = True
        if not parent_done and now >= start + policy['parent']:
            # A direct child cannot mistake a recycled PID for its still-live parent.
            if os.getppid() == parent:
                os.kill(parent, signal.SIGKILL)
                emit('SUPERVISOR_SIGKILL_SENT', supervisor_pid=parent)
            parent_done = True
            return
        ready, _, _ = select.select([control], [], [], 0.01)
        if ready:
            raw = os.read(control, 32)
            if b'K' in raw:
                # Supervisor already killed the still-unreaped group leader.
                group_done = True
                emit('GROUP_KILL_ACKNOWLEDGED')
            if b'F' in raw:
                require(group_done, 'Disarm before owned group cleanup forbidden')
                emit('WATCHDOG_DISARMED_AFTER_RECEIPT')
                return
            if not raw:
                if not group_done:
                    try:
                        os.killpg(pgid, signal.SIGKILL)
                        emit('GROUP_SIGKILL_ON_PARENT_PIPE_CLOSE', pgid=pgid)
                    except ProcessLookupError:
                        pass
                return


if __name__ == '__main__':
    main()
