#!/usr/bin/env python3
"""Single-use owned-process guard; 165s work plus <=10s cleanup, hard180s."""
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
signal.alarm(180)
resource.setrlimit(resource.RLIMIT_CPU, (2, 2))
start = time.monotonic()
if (ROOT / "execution_receipt.json").exists() or (ROOT / "results.json").exists():
    raise SystemExit("Single-use guard: existing receipt/results; no automatic rerun allowed")
registration_bytes = (ROOT / "registration.json").read_bytes()
registration = json.loads(registration_bytes)
for name, expected in registration["frozen_local_sha256"].items():
    actual = hashlib.sha256((ROOT / name).read_bytes()).hexdigest()
    if actual != expected:
        raise SystemExit("Source/input hash mismatch: " + name)
for item in registration["readonly_source_pins"]:
    actual = hashlib.sha256(Path(item["path"]).read_bytes()).hexdigest()
    if actual != item["sha256"]:
        raise SystemExit("Readonly source mismatch: " + item["path"])
receipt = {
    "schema": "A076-single-execution-receipt-v1",
    "status": "STARTING", "guard_pid": os.getpid(),
    "started_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
    "registration_sha256": hashlib.sha256(registration_bytes).hexdigest(),
    "command": [sys.executable, str(ROOT / "check_fixed.py")],
    "total_hard_wall_seconds_including_cleanup": 180,
    "work_wait_seconds": 165, "cleanup_wait_seconds": 10,
    "child_cpu_limit_seconds": 27, "guard_cpu_limit_seconds": 2,
    "attempt": 1, "automatic_retries": 0,
    "process_group_ownership": "new child process group; only this group may be killed"
}


def save():
    receipt["elapsed_wall_seconds"] = time.monotonic() - start
    (ROOT / "execution_receipt.json").write_text(json.dumps(receipt, indent=2) + "\n")


save()
proc = subprocess.Popen(receipt["command"], cwd=ROOT, stdout=subprocess.PIPE,
                        stderr=subprocess.PIPE, text=True, start_new_session=True,
                        env={**os.environ, "PYTHONDONTWRITEBYTECODE": "1"})
receipt["child_pid"] = proc.pid
receipt["owned_child_pgid"] = proc.pid
receipt["status"] = "RUNNING"
save()
try:
    stdout, stderr = proc.communicate(timeout=165)
    receipt["status"] = "COMPLETE" if proc.returncode == 0 else "FAILED_FIRST_STOP"
except subprocess.TimeoutExpired:
    receipt["status"] = "TIMEOUT_KILL_OWNED_GROUP_NO_RETRY"
    save()
    os.killpg(proc.pid, signal.SIGKILL)
    stdout, stderr = proc.communicate(timeout=10)
receipt["returncode"] = proc.returncode
receipt["stdout"] = stdout
receipt["stderr"] = stderr
receipt["rusage_children_user_seconds"] = resource.getrusage(resource.RUSAGE_CHILDREN).ru_utime
receipt["rusage_children_system_seconds"] = resource.getrusage(resource.RUSAGE_CHILDREN).ru_stime
receipt["guard_cpu_seconds"] = time.process_time()
receipt["completed_utc"] = datetime.datetime.now(datetime.timezone.utc).isoformat()
save()
print(json.dumps(receipt))
raise SystemExit(0 if receipt["status"] == "COMPLETE" else 1)
