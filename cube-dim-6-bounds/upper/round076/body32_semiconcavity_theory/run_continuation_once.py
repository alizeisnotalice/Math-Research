#!/usr/bin/env python3
"""Explicitly released continuation after the recorded zero-science startup failure."""
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
# Preserve a hard ceiling that the child may lower to27; guard soft budget2.
resource.setrlimit(resource.RLIMIT_CPU, (2, 29))
start = time.monotonic()
if (ROOT / "continuation_receipt.json").exists() or (ROOT / "results.json").exists():
    raise SystemExit("Single-use continuation: receipt/results already exist; no automatic rerun")
supplement_bytes = (ROOT / "supplement_registration.json").read_bytes()
supplement = json.loads(supplement_bytes)
for name, expected in supplement["pinned_local_sha256"].items():
    if hashlib.sha256((ROOT / name).read_bytes()).hexdigest() != expected:
        raise SystemExit("Supplement hash mismatch: " + name)
registration_bytes = (ROOT / "registration.json").read_bytes()
registration = json.loads(registration_bytes)
for name, expected in registration["frozen_local_sha256"].items():
    if hashlib.sha256((ROOT / name).read_bytes()).hexdigest() != expected:
        raise SystemExit("Original source/input hash mismatch: " + name)
for item in registration["readonly_source_pins"]:
    if hashlib.sha256(Path(item["path"]).read_bytes()).hexdigest() != item["sha256"]:
        raise SystemExit("Readonly source mismatch: " + item["path"])
receipt = {
    "schema": "A076-explicit-continuation-receipt-v1",
    "status": "STARTING", "guard_pid": os.getpid(),
    "started_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
    "supplement_registration_sha256": hashlib.sha256(supplement_bytes).hexdigest(),
    "original_registration_sha256": hashlib.sha256(registration_bytes).hexdigest(),
    "prior_failed_receipt_sha256": supplement["pinned_local_sha256"]["execution_receipt.json"],
    "command": [sys.executable, str(ROOT / "check_fixed.py")],
    "total_hard_wall_seconds_including_cleanup": 180,
    "work_wait_seconds": 165, "cleanup_wait_seconds": 10,
    "child_cpu_soft_hard_seconds": [27, 27],
    "guard_cpu_soft_hard_seconds": [2, 29],
    "overall_invocation_number": 2, "explicit_continuation_number": 1,
    "prior_startup_failure_invocations": 1, "prior_scientific_configurations": 0,
    "automatic_retries": 0,
    "process_group_ownership": "new child process group; only this group may be killed"
}


def save():
    receipt["elapsed_wall_seconds"] = time.monotonic() - start
    (ROOT / "continuation_receipt.json").write_text(json.dumps(receipt, indent=2) + "\n")


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
