"""One approved diagnostic only; external 60s hard process-group deadline."""
from pathlib import Path
from time import monotonic,time
from hashlib import sha256
import json,os,signal,subprocess,sys
root=Path(__file__).resolve().parent
meta=root/'quadrants.run.json';log=root/'quadrants.log';script=root/'quadrants_once.py'
if any(p.exists() for p in (meta,log,root/'result.json',root/'quadrants.json')):
 raise SystemExit('Existing evidence: do not repeat or overwrite.')
start=monotonic();record={'status':'started','budget_seconds':60,'started_unix':time(),
 'command_relative':['python','-B','quadrants_once.py'],'source_sha256_before':sha256(script.read_bytes()).hexdigest()}
meta.write_text(json.dumps(record,indent=2)+'\n')
with log.open('w') as stream:
 p=subprocess.Popen([sys.executable,'-B',str(script)],cwd=root,stdout=stream,stderr=subprocess.STDOUT,start_new_session=True)
 record['pid']=p.pid;meta.write_text(json.dumps(record,indent=2)+'\n')
 try:
  code=p.wait(timeout=max(.001,60-(monotonic()-start)))
  record.update(status='terminal',exit_code=code)
 except subprocess.TimeoutExpired:
  os.killpg(p.pid,signal.SIGKILL);p.wait()
  record.update(status='timed_out',exit_code=p.returncode)
try:os.killpg(p.pid,signal.SIGKILL)
except ProcessLookupError:pass
record.update(elapsed_seconds=monotonic()-start,log_sha256=sha256(log.read_bytes()).hexdigest(),
 source_sha256_after=sha256(script.read_bytes()).hexdigest())
meta.write_text(json.dumps(record,indent=2)+'\n');print(json.dumps(record))
