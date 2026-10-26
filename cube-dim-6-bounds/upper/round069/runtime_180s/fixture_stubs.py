"""Six named runtime-only configurations. No body, jet, flint or oracle imports."""
import signal
import time
import os
from pathlib import Path
import subprocess
import sys
from runtime_support import P, FIXTURE, PRODUCTION, NAMES, require, write_new, validate_operation, bounded_reap, production_gate


def launch(name, role):
    return subprocess.Popen([sys.executable, '-B', str(P/'fixture_stubs.py'), name, role])


def child(name, role):
    folder = P/'fixtures'/name
    write_new(folder/f'{role}.json', dict(pid=os.getpid(), ppid=os.getppid(), pgid=os.getpgrp(), role=role))
    if role=='tree_child':
        descendant = launch(name, 'grandchild')
        # Intentionally not reaped by this fault fixture; group death is not a reap claim.
        write_new(folder/'grandchild_handle.json', dict(pid=descendant.pid, pgid=os.getpgrp()))
    if name=='normal_child':
        time.sleep(0.015)
        return
    while True:
        time.sleep(0.03)


def run(name, start, policy):
    require(name in NAMES and policy == FIXTURE, 'Only fixed isolated fixture mode')
    folder = P/'fixtures'/name
    if name=='gate_and_deadline_isolation':
        bad = [['native', 'r.a0.n000001', '--timeout', '999'], ['control', 'status', '999'],
               ['control', 'reset'], ['native', 'not-a-reservation']]
        rejected = []
        for args in bad:
            try:
                validate_operation(args)
            except ValueError:
                rejected.append(args)
        require(len(rejected)==4 and PRODUCTION['soft']==170 and PRODUCTION['group']==178 and
                PRODUCTION['supervisor']<180 and PRODUCTION['reap']==1, 'Deadline/argument isolation failed')
        require(not (P/'runtime_release.json').exists(), 'Preparation production gate unexpectedly open')
        try:
            production_gate(['control', 'status'])
        except FileNotFoundError as error:
            require(Path(error.filename)==P/'runtime_release.json', 'Wrong gate rejection')
        else:
            raise ValueError('Unauthorized production gate accepted')
        require('core_v2' not in sys.modules and 'math_contract' not in sys.modules,
                'Fixture crossed into mathematical task imports')
        return dict(rejected_operations=rejected, production_release_present=False,
                    actual_production_gate_rejected=True, production_deadlines=PRODUCTION, no_body_import=True)
    if name=='receipt_write_hang':
        return dict(work_completed=True, next='SUPERVISOR_REAL_FIFO_OPEN_BLOCK')
    role = 'tree_child' if name=='descendant_group' else 'child'
    process = launch(name, role)  # no setsid: inherits worker PG
    write_new(folder/'stub_handle.json', dict(pid=process.pid, worker_pid=os.getpid(), pgid=os.getpgrp()))
    try:
        code = process.wait(timeout=max(0.001, start+policy['soft']-time.monotonic()+0.2))
        return dict(child_returncode=code, child_reaped=True, worker_pgid=os.getpgrp())
    except (TimeoutError, subprocess.TimeoutExpired) as error:
        write_new(folder/'soft_observed.json', dict(error=repr(error), elapsed=time.monotonic()-start))
        if name in ('cleanup_hang', 'descendant_group'):
            os.mkfifo(folder/'blocked_cleanup.fifo')
            with (folder/'blocked_cleanup.fifo').open('w') as stream:
                stream.write('cleanup intentionally blocked')
        process.kill()  # original direct child only
        outcome = bounded_reap(process, time.monotonic()+policy['reap'])
        write_new(folder/'stub_cleanup.json', dict(**outcome, kill_sent=True, elapsed=time.monotonic()-start))
        raise


if __name__=='__main__':
    signal.signal(signal.SIGALRM, signal.SIG_DFL)
    signal.setitimer(signal.ITIMER_REAL, 1.55)
    require(len(sys.argv)==3 and sys.argv[1] in NAMES, 'Fixed stub invocation only')
    child(*sys.argv[1:])
