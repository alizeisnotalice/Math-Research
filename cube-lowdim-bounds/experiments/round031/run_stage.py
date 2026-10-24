"""One bounded run per stage. Never retries or controls unrelated processes."""
from pathlib import Path
import sys,subprocess,time,json,os,signal,hashlib
root=Path(__file__).resolve().parent
stage=sys.argv[1] if len(sys.argv)==2 else ''
assert stage in ('build','search')
budget=300 if stage=='build' else 180
source=root/('build_cells.py' if stage=='build' else 'search.py')
meta=root/(stage+'.run.json');log=root/(stage+'.log')
result=root/('build_result.json' if stage=='build' else 'result.json')
if any(p.exists() for p in (meta,log,result)):raise SystemExit('Existing evidence: do not repeat.')
start=time.monotonic();record={'stage':stage,'status':'started','hard_seconds':budget,
 'command':[sys.executable,'-B',str(source)],'source_sha256_before':hashlib.sha256(source.read_bytes()).hexdigest()}
meta.write_text(json.dumps(record,indent=2)+'\n')
with log.open('x') as stream:
 proc=subprocess.Popen(record['command'],cwd=root,stdout=stream,stderr=subprocess.STDOUT,start_new_session=True)
 record['owned_pgid']=proc.pid;meta.write_text(json.dumps(record,indent=2)+'\n')
 try:
  proc.wait(timeout=max(.001,budget-(time.monotonic()-start)));record['status']='returned'
 except subprocess.TimeoutExpired:
  os.killpg(proc.pid,signal.SIGKILL);proc.wait();record['status']='hard_timeout_killed'
try:os.killpg(proc.pid,0)
except ProcessLookupError:record['cleanup']='owned_group_absent'
else:os.killpg(proc.pid,signal.SIGKILL);record['cleanup']='remaining_owned_group_killed'
record.update(exit_code=proc.returncode,wall_seconds=time.monotonic()-start,
 source_sha256_after=hashlib.sha256(source.read_bytes()).hexdigest(),log_sha256=hashlib.sha256(log.read_bytes()).hexdigest())
meta.write_text(json.dumps(record,indent=2)+'\n');print(json.dumps(record,indent=2))
