"""Four fixed runtime revision scenarios; never imports body/math modules."""
import json
import os
import signal
import subprocess
import sys
import time
from pathlib import Path
from runtime_support import P, FIXTURE, PRODUCTION, NAMES, require, write_new, production_gate, validate_operation


def block_full_pipe(folder, phase, start):
    reader, writer = os.pipe()
    os.set_blocking(writer, False)
    count = 0
    while True:
        try:
            count += os.write(writer, b'x'*4096)
        except BlockingIOError:
            break
    os.set_blocking(writer, True)
    write_new(folder/f'{phase}_full_pipe.json', dict(reader_fd=reader, writer_fd=writer,
              retained_reader=True, drained=False, filled_bytes=count, pid=os.getpid(),
              pgid=os.getpgrp(), elapsed=time.monotonic()-start))
    os.write(writer, b'x')  # reader is retained and never drained; this must block
    write_new(folder/f'{phase}_unexpected_return.json', dict(error='Full pipe write unexpectedly returned'))
    raise RuntimeError('Full pipe did not remain blocked')


def run(name, start, policy):
    require(name in NAMES and policy==FIXTURE, 'Fixed fixture mode only')
    folder=P/'fixtures'/name
    if name=='retained_reader_cleanup':
        try:
            while True:
                signal.pause()
        except TimeoutError:
            block_full_pipe(folder,'cleanup',start)
    if name=='early_return_descendants':
        child=subprocess.Popen([sys.executable,'-B',str(P/'fixture_stubs.py'),'child',name])
        write_new(folder/'spawn_handle.json',dict(child_pid=child.pid, worker_pid=os.getpid(), worker_pgid=os.getpgrp()))
        while not (folder/'grandchild.json').exists() and time.monotonic()<start+policy['soft']-0.03:
            time.sleep(0.002)
        require((folder/'grandchild.json').exists(),'Actual grandchild did not reach readiness')
        return dict(early_return_with_descendants=True,child_pid=child.pid,worker_pgid=os.getpgrp())
    if name=='receipt_retained_reader_hang':
        return dict(work_completed=True,next='SUPERVISOR_RETAINED_READER_FULL_PIPE')
    if name=='receipt_error_gate_unacked_eperm':
        rejected=[]
        for args in (['native','r.a0.n000001','--timeout','999'],['control','reset'],['control','status','999']):
            try:
                validate_operation(args)
            except ValueError:
                rejected.append(args)
        require(len(rejected)==3,'Argument override accepted')
        try:
            production_gate(['control','status'])
        except FileNotFoundError as error:
            require(Path(error.filename)==P/'runtime_release.json','Wrong production gate failure')
        else:
            raise ValueError('Absent release accepted')
        require('core_v2' not in sys.modules and 'math_contract' not in sys.modules,'Gate crossed math boundary')
        return dict(actual_absent_release_rejected=True,invalid_arguments_rejected=rejected,
                    production_policy=PRODUCTION,no_math_import=True,
                    next='INJECT_UNACKED_EPERM_AND_EISDIR_RECEIPT_WRITE_ERROR')


def auxiliary(role,name):
    signal.signal(signal.SIGALRM,signal.SIG_DFL)
    signal.setitimer(signal.ITIMER_REAL,1.7)
    if role=='sentinel':
        path=P/'sentinel_started.json'
    else:
        require(name=='early_return_descendants','Only registered descendants')
        path=P/'fixtures'/name/f'{role}.json'
    write_new(path,dict(pid=os.getpid(),ppid=os.getppid(),pgid=os.getpgrp(),role=role))
    if role=='child':
        subprocess.Popen([sys.executable,'-B',str(P/'fixture_stubs.py'),'grandchild',name])
    while True:
        time.sleep(0.03)


if __name__=='__main__':
    require(len(sys.argv)==3 and sys.argv[1] in ('child','grandchild','sentinel'),'Registered auxiliary only')
    auxiliary(*sys.argv[1:])
