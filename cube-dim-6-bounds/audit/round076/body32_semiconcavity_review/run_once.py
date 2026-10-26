#!/usr/bin/env python3
"""Single attempt: <=165s work, <=3s kill/reap, hard wall alarm179s."""
import datetime
import hashlib
import json
import os
from pathlib import Path
import resource
import signal
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parent
start = time.monotonic()
signal.alarm(179)
receipt_path = ROOT/'execution_receipt.json'
fd = os.open(receipt_path, os.O_CREAT | os.O_EXCL | os.O_WRONLY, 0o644)
os.close(fd)
receipt = {'schema': 'B076-single-attempt-receipt-v1', 'attempt': 1,
           'guard_pid': os.getpid(), 'started_utc': datetime.datetime.now(datetime.timezone.utc).isoformat(),
           'hard_wall_seconds': 179, 'program_total_ceiling_seconds': 180,
           'work_timeout_seconds': 165, 'cleanup_timeout_seconds': 3,
           'guard_cpu_limit_seconds': 2, 'worker_cpu_limit_seconds': 27,
           'automatic_retries': 0, 'status': 'PRECHECK'}


def save():
    receipt['elapsed_wall_seconds'] = time.monotonic()-start
    receipt_path.write_text(json.dumps(receipt, indent=2)+'\n')


def child_limit():
    resource.setrlimit(resource.RLIMIT_CPU, (27, 27))


save()
try:
    reg_bytes = (ROOT/'registration.json').read_bytes()
    reg = json.loads(reg_bytes)
    if reg['state'] != 'SEALED_READY':
        raise ValueError('Source seal is not ready')
    for name, expected in reg['frozen_local_sha256'].items():
        if hashlib.sha256((ROOT/name).read_bytes()).hexdigest() != expected:
            raise ValueError('Local hash mismatch: '+name)
    for pin in reg['readonly_source_pins']:
        if hashlib.sha256(Path(pin['path']).read_bytes()).hexdigest() != pin['sha256']:
            raise ValueError('Readonly source hash mismatch: '+pin['path'])
    receipt['registration_sha256'] = hashlib.sha256(reg_bytes).hexdigest()
    command = [sys.executable, '-B', str(ROOT/'check_parallel_quotient.py')]
    receipt['command'] = command
    # Set the guard's low hard CPU limit only AFTER spawning. This avoids
    # passing a hard limit of2 to a worker that must set a limit of27.
    proc = subprocess.Popen(command, cwd=ROOT, start_new_session=True,
                            stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True,
                            preexec_fn=child_limit)
    resource.setrlimit(resource.RLIMIT_CPU, (2, 2))
    receipt.update(child_pid=proc.pid, owned_process_group=proc.pid, status='RUNNING')
    save()
    try:
        out, err = proc.communicate(timeout=165)
        receipt['status'] = 'PASS' if proc.returncode == 0 else 'FAILED_FIRST_STOP'
    except subprocess.TimeoutExpired:
        os.killpg(proc.pid, signal.SIGKILL)
        out, err = proc.communicate(timeout=3)
        receipt['status'] = 'TIMEOUT_FIRST_STOP'
    receipt.update(exit_code=proc.returncode, stdout=out, stderr=err)
except Exception as exc:
    receipt.update(status='PRECHECK_OR_ENVIRONMENT_FAILED_FIRST_STOP', error=repr(exc))
receipt['guard_cpu_seconds'] = time.process_time()
usage = resource.getrusage(resource.RUSAGE_CHILDREN)
receipt['children_cpu_seconds'] = usage.ru_utime+usage.ru_stime
receipt['completed_utc'] = datetime.datetime.now(datetime.timezone.utc).isoformat()
save()
print(json.dumps(receipt))
raise SystemExit(0 if receipt['status'] == 'PASS' else 1)
