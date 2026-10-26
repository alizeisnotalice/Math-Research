"""Own process group only; work170s + cleanup10s <=180s, no retry."""
from pathlib import Path
import subprocess, sys, os, signal, time, json
H = Path(__file__).resolve().parent
receipt = H / 'process_receipt.json'
if receipt.exists():
    raise RuntimeError('Single invocation only: receipt already exists')
start = time.monotonic()
p = subprocess.Popen([sys.executable, '-B', str(H / 'check_saved.py')], cwd=H,
                     stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                     text=True, start_new_session=True)
r = {'status': 'RUNNING', 'guard_pid': os.getpid(), 'owned_pid': p.pid,
     'owned_process_group': p.pid, 'wall_limit_including_cleanup_seconds': 180,
     'scientific_work_timeout_seconds': 170, 'cleanup_reserve_seconds': 10,
     'automatic_retries': 0, 'kill_scope': 'owned child process group only'}
receipt.write_text(json.dumps(r, indent=2) + '\n')
try:
    output, _ = p.communicate(timeout=max(0, 170 - (time.monotonic() - start)))
except subprocess.TimeoutExpired:
    r['timeout'] = True
    os.killpg(p.pid, signal.SIGKILL)
    output, _ = p.communicate(timeout=max(0.001, 180 - (time.monotonic() - start)))
r.update(status='EXITED_AND_REAPED', exit_code=p.returncode,
         wall_seconds_including_cleanup=time.monotonic() - start, child_reaped=True)
receipt.write_text(json.dumps(r, indent=2) + '\n')
print(output, end='')
print(json.dumps(r))
raise SystemExit(p.returncode)
