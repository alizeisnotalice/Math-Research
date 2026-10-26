#!/usr/bin/env python3
"""One local released batch; 180s total including its own process cleanup."""
import datetime,hashlib,json,os,resource,signal,subprocess,sys,time
from pathlib import Path
ROOT=Path(__file__).resolve().parent
signal.alarm(180);resource.setrlimit(resource.RLIMIT_CPU,(2,29));start=time.monotonic()
if (ROOT/"execution_receipt.json").exists() or (ROOT/"results.json").exists():raise SystemExit("Existing receipt/results: no automatic retry")
raw=(ROOT/"registration.json").read_bytes();reg=json.loads(raw)
for name,h in reg["frozen_local_sha256"].items():
    if hashlib.sha256((ROOT/name).read_bytes()).hexdigest()!=h:raise SystemExit("Frozen local hash mismatch: "+name)
for pin in reg["readonly_source_pins"]:
    if hashlib.sha256(Path(pin["path"]).read_bytes()).hexdigest()!=pin["sha256"]:raise SystemExit("Readonly source mismatch: "+pin["path"])
receipt={"schema":"A093-single-execution-receipt-v1","status":"STARTING","guard_pid":os.getpid(),"started_utc":datetime.datetime.now(datetime.timezone.utc).isoformat(),
 "registration_sha256":hashlib.sha256(raw).hexdigest(),"command":[sys.executable,str(ROOT/"check_fixed.py")],"total_wall_limit_including_cleanup_seconds":180,
 "work_wait_seconds":165,"cleanup_wait_seconds":10,"guard_cpu_soft_hard_seconds":[2,29],"child_cpu_soft_hard_seconds":[27,27],"attempt":1,"automatic_retries":0}
def save():
    receipt["elapsed_wall_seconds"]=time.monotonic()-start;(ROOT/"execution_receipt.json").write_text(json.dumps(receipt,indent=2)+"\n")
save()
p=subprocess.Popen(receipt["command"],cwd=ROOT,stdout=subprocess.PIPE,stderr=subprocess.PIPE,text=True,start_new_session=True,env={**os.environ,"PYTHONDONTWRITEBYTECODE":"1"})
receipt.update(child_pid=p.pid,owned_child_pgid=p.pid,status="RUNNING");save()
try:
    stdout,stderr=p.communicate(timeout=165);receipt["status"]="COMPLETE" if p.returncode==0 else "FAILED_FIRST_STOP"
except subprocess.TimeoutExpired:
    receipt["status"]="TIMEOUT_KILL_OWNED_GROUP_NO_RETRY";save();os.killpg(p.pid,signal.SIGKILL);stdout,stderr=p.communicate(timeout=10)
receipt.update(returncode=p.returncode,stdout=stdout,stderr=stderr,child_user_cpu_seconds=resource.getrusage(resource.RUSAGE_CHILDREN).ru_utime,
 child_system_cpu_seconds=resource.getrusage(resource.RUSAGE_CHILDREN).ru_stime,guard_cpu_seconds=time.process_time(),completed_utc=datetime.datetime.now(datetime.timezone.utc).isoformat())
save();print(json.dumps(receipt));raise SystemExit(0 if receipt["status"]=="COMPLETE" else 1)
