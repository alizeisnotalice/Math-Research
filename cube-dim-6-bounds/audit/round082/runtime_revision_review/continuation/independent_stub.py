"""Three independently registered combined runtime faults; no math imports."""
import json
import os
from pathlib import Path
import signal
import sys
import time

P = Path(__file__).resolve().parent


def save(path, value):
    with path.open('x') as stream:
        json.dump(value, stream, sort_keys=True)
        stream.write('\n')


def block_full_pipe(folder, phase, start):
    reader, writer = os.pipe()
    os.set_blocking(writer, False)
    count = 0
    while True:
        try:
            count += os.write(writer, b'p' * 8192)
            assert count < 8 * 1024 * 1024
        except BlockingIOError:
            break
    os.set_blocking(writer, True)
    save(folder / (phase + '_blocked.json'), dict(pid=os.getpid(), pgid=os.getpgrp(),
        elapsed=time.monotonic()-start, retained_reader=True, filled_bytes=count,
        reader_fd=reader, writer_fd=writer))
    os.write(writer, b'blocked')
    save(folder / (phase + '_unexpected_return.json'), dict(error='Unexpected drain'))
    raise AssertionError('Full pipe drained')


def run(name, start, policy):
    case = (P / 'independent_case.txt').read_text().strip()
    folder = P / 'fixtures' / name
    save(folder / 'entered.json', dict(case=case, pid=os.getpid(), elapsed=time.monotonic()-start))
    if case == 'A3_dead_watchdog_plus_supervisor_receipt_hang':
        return dict(case=case, next='Original supervisor receipt hook')
    if case == 'A1_wrong_pipe_ack_plus_cleanup_hang':
        handles = json.loads((folder / 'handles.json').read_text())
        forged = dict(kind='GROUP_SIGKILL_SENT', pid=handles['watchdog_pid'],
                      pgid=handles['worker_pgid'], from_untrusted_result_pipe=True)
        # The original worker's result descriptor, never the watchdog event pipe.
        os.write(int(sys.argv[5]), (json.dumps(forged)+'\n').encode())
        save(folder / 'forged_message.json', forged)
    elif case == 'A2_dead_watchdog_plus_cleanup_and_receipt_exception':
        (folder / 'receipt.json').mkdir()
    else:
        raise AssertionError('Unregistered independent combination')
    try:
        while True:
            signal.pause()
    except TimeoutError:
        block_full_pipe(folder, 'cleanup', start)
