#!/usr/bin/env python3
"""Launch exactly one owned scientific child with a hard180-second wall limit."""
import sys
sys.dont_write_bytecode=True
from pathlib import Path
import subprocess,os,signal,time,json,datetime
H=Path(__file__).resolve().parent
assert len(sys.argv)==2 and sys.argv[1]in ['producer','checker']
mode=sys.argv[1];path=H/(mode+'_process_receipt.json')
assert not path.exists(),'no automatic or silent rerun'
cmd=[sys.executable,'-B',str(H/('evaluate.py'if mode=='producer'else'check_saved.py'))]
if mode=='producer':cmd+=['--run']
receipt=dict(mode=mode,status='STARTING',wall_timeout_seconds=180,started_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),command=cmd,automatic_retries=0)
def save():path.write_text(json.dumps(receipt,indent=2)+'\n')
save();tic=time.monotonic();p=None
try:
    p=subprocess.Popen(cmd,cwd=H,stdout=subprocess.PIPE,stderr=subprocess.PIPE,text=True,start_new_session=True)
    receipt.update(status='RUNNING',owned_pid=p.pid,owned_process_group=p.pid);save()
    try:
        out,err=p.communicate(timeout=180)
        receipt.update(status='COMPLETED'if p.returncode==0 else'FAILED',exit_code=p.returncode)
    except subprocess.TimeoutExpired:
        os.killpg(p.pid,signal.SIGKILL)
        out,err=p.communicate()
        receipt.update(status='TIMED_OUT_OWN_GROUP_KILLED_AND_REAPED',exit_code=124,child_exit_code=p.returncode,checkpoint='retain child results.json/run_receipt.json exactly as written')
    receipt.update(stdout=out,stderr=err,elapsed_seconds=time.monotonic()-tic,finished_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),child_reaped=True,active_handles=[])
    save();print(json.dumps({k:receipt[k]for k in ['mode','status','exit_code','elapsed_seconds','owned_pid']}))
    if out:print(out,end='')
    if err:print(err,end='',file=sys.stderr)
    raise SystemExit(receipt['exit_code'])
except BaseException as e:
    if p is not None and p.poll()is None:
        os.killpg(p.pid,signal.SIGKILL);p.communicate()
        receipt.update(status='INTERRUPTED_OWN_GROUP_KILLED_AND_REAPED',error=repr(e),child_reaped=True,elapsed_seconds=time.monotonic()-tic,active_handles=[]);save()
    raise
