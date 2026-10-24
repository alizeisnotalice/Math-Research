#!/usr/bin/env python3
"""One diagnostic process group, 180s hard including imports; no retries."""
from hashlib import sha256
import json
import os
from pathlib import Path
import signal
import subprocess
import sys
from time import monotonic

ROOT=Path(__file__).resolve().parent


def main():
    paths=[ROOT/n for n in ('diagnostic.run.json','diagnostic.log','result.json','calls.jsonl')]
    if any(p.exists() for p in paths):raise SystemExit('Existing diagnostic evidence; repetition forbidden.')
    source=ROOT/'diagnose.py';start=monotonic()
    record={'status':'started','hard_seconds':180,'source_sha256':sha256(source.read_bytes()).hexdigest(),
            'command_relative':['python','-B','diagnose.py']}
    paths[0].write_text(json.dumps(record,indent=2)+'\n')
    with paths[1].open('x') as stream:
        process=subprocess.Popen([sys.executable,'-B',str(source)],cwd=ROOT,
            stdout=stream,stderr=subprocess.STDOUT,start_new_session=True)
        record.update(pid=process.pid,owned_process_group=process.pid)
        paths[0].write_text(json.dumps(record,indent=2)+'\n')
        try:
            process.wait(timeout=max(.001,180-(monotonic()-start)))
            record['status']='returned'
        except subprocess.TimeoutExpired:
            os.killpg(process.pid,signal.SIGKILL);process.wait();record['status']='hard_timeout_killed'
        finally:
            try:os.killpg(process.pid,signal.SIGKILL);record['owned_group_cleanup']='remaining_owned_group_killed'
            except ProcessLookupError:record['owned_group_cleanup']='no_remaining_owned_group'
    record.update(exit_code=process.returncode,wall_seconds=monotonic()-start,
        source_sha256_after=sha256(source.read_bytes()).hexdigest(),log_sha256=sha256(paths[1].read_bytes()).hexdigest())
    record['source_unchanged']=record['source_sha256']==record['source_sha256_after']
    if paths[2].is_file():
        result=json.loads(paths[2].read_text());record.update(result_status=result['status'],counts=result['counts'])
    paths[0].write_text(json.dumps(record,indent=2)+'\n');print(json.dumps(record,indent=2),flush=True)


if __name__=='__main__':main()
