"""A069 new fault injections; never a body kernel or mathematical fixture."""
import json
import os
from pathlib import Path
import signal
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parent


def put(path, value):
    with path.open('x') as stream:
        json.dump(value, stream, sort_keys=True)
        stream.write('\n')


def descendant(role, folder):
    # Short independent test fuse; the tested group kill must precede it.
    signal.signal(signal.SIGALRM, signal.SIG_DFL)
    signal.alarm(4)
    put(folder / (role + '.json'), dict(pid=os.getpid(), pgid=os.getpgrp(),
                                      ppid=os.getppid(), role=role))
    if role == 'child':
        subprocess.Popen([sys.executable, '-B', __file__, 'grandchild', str(folder)],
                         start_new_session=False)
    while True:
        time.sleep(0.02)


def run(spec, start, policy):
    case = (ROOT / 'independent_case.txt').read_text().strip()
    folder = ROOT / 'fixtures' / spec
    put(folder / 'injection_entered.json', dict(case=case, pid=os.getpid(),
        pgid=os.getpgrp(), elapsed=time.monotonic() - start))
    if case in ('R1_retained_reader_full_pipe', 'R5_watchdog_death_during_cleanup'):
        try:
            time.sleep(10)  # interrupted by the unmodified worker's soft handler
        except TimeoutError:
            # Keep BOTH endpoints live. There is a reader, so this cannot turn
            # into the author's BrokenPipeError. Nobody ever drains that reader.
            reader, writer = os.pipe()
            os.set_blocking(writer, False)
            filled = 0
            while True:
                try:
                    filled += os.write(writer, b'x' * 16384)
                    assert filled < 8 * 1024 * 1024, 'Unexpected pipe capacity'
                except BlockingIOError:
                    break
            os.set_blocking(writer, True)
            put(folder / 'cleanup_block_entered.json', dict(pid=os.getpid(),
                reader_fd=reader, writer_fd=writer, retained_reader=True,
                filled_bytes=filled, elapsed=time.monotonic() - start))
            os.write(writer, b'blocked-after-full' * 1024)
            put(folder / 'cleanup_returned.json', dict(unexpected=True))
            raise AssertionError('Full pipe unexpectedly drained')
    elif case == 'R2_early_return_descendants_and_sentinel':
        child = subprocess.Popen([sys.executable, '-B', __file__, 'child', str(folder)],
                                 start_new_session=False)
        while not (folder / 'grandchild.json').exists():
            if time.monotonic() >= start + policy['soft'] - 0.03:
                raise TimeoutError('Descendant marker not ready within registered fixture')
            time.sleep(0.002)
        put(folder / 'early_return.json', dict(child_pid=child.pid,
            child_poll=child.poll(), elapsed=time.monotonic() - start))
        return dict(case=case, leaves_live_descendants=True)
    elif case == 'R3_supervisor_receipt_fifo':
        return dict(case=case, original_envelope_fifo_branch=True)
    elif case == 'R4_receipt_path_exception':
        (folder / 'receipt.json').mkdir()
        return dict(case=case, receipt_path_is_directory=True)
    raise AssertionError('Unregistered independent case')


if __name__ == '__main__':
    descendant(sys.argv[1], Path(sys.argv[2]))
