#!/usr/bin/env python3
"""Own-child170-second work timeout with10 seconds reserved for cleanup."""
import sys
sys.dont_write_bytecode=True
from pathlib import Path
import subprocess,os,signal,time,json,datetime
H=Path(__file__).resolve().parent
assert sys.argv[1:]in [['producer'],['checker']]
mode=sys.argv[1];path=H/(mode+'_process_receipt.json');assert not path.exists(),'no rerun'
cmd=[sys.executable,'-B',str(H/('check_fixed.py'if mode=='producer'else'verify_saved.py'))]
tic=time.monotonic();receipt=dict(mode=mode,status='STARTING',work_timeout_seconds=170,total_wall_limit_including_cleanup_seconds=180,cleanup_reserve_seconds=10,started_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),command=cmd,automatic_retries=0)
def save():path.write_text(json.dumps(receipt,indent=2)+'\n')
save();p=None
try:
    p=subprocess.Popen(cmd,cwd=H,stdout=subprocess.PIPE,stderr=subprocess.PIPE,text=True,start_new_session=True)
    receipt.update(status='RUNNING',owned_pid=p.pid,owned_process_group=p.pid);save()
    try:
        out,err=p.communicate(timeout=max(1,170-(time.monotonic()-tic)))
        receipt.update(status='COMPLETED'if p.returncode==0 else'FAILED',exit_code=p.returncode)
    except subprocess.TimeoutExpired:
        os.killpg(p.pid,signal.SIGKILL)
        out,err=p.communicate(timeout=max(1,178-(time.monotonic()-tic)))
        receipt.update(status='TIMED_OUT_OWN_GROUP_KILLED_AND_REAPED',exit_code=124,child_exit_code=p.returncode,checkpoint='retain existing results and run receipt')
    receipt.update(stdout=out,stderr=err,elapsed_seconds=time.monotonic()-tic,child_reaped=True,active_handles=[],finished_utc=datetime.datetime.now(datetime.timezone.utc).isoformat())
    receipt['within180_including_cleanup']=receipt['elapsed_seconds']<180
    save();print(json.dumps({k:receipt[k]for k in ['mode','status','exit_code','elapsed_seconds','owned_pid','within180_including_cleanup']}))
    if out:print(out,end='')
    if err:print(err,end='',file=sys.stderr)
    raise SystemExit(receipt['exit_code'])
except BaseException as ex:
    if p is not None and p.poll()is None:
        os.killpg(p.pid,signal.SIGKILL);p.communicate(timeout=max(1,179-(time.monotonic()-tic)))
        receipt.update(status='INTERRUPTED_OWN_GROUP_KILLED_AND_REAPED',error=repr(ex),child_reaped=True,elapsed_seconds=time.monotonic()-tic,active_handles=[]);save()
    raise
