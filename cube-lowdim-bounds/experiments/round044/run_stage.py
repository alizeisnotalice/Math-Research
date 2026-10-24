"""One fixed stage only: build 300s, solve 180s process-group deadline."""
from pathlib import Path
from time import monotonic,time
from hashlib import sha256
import json,os,signal,subprocess,sys
root=Path(__file__).resolve().parent
stage=sys.argv[1];assert stage in ('build','solve')
seconds=300 if stage=='build' else 180
meta=root/(stage+'.run.json');log=root/(stage+'.log');script=root/('build_directions.py' if stage=='build' else 'solve.py')
if any(p.exists() for p in (meta,log,root/(stage+'_result.json'))):
 raise SystemExit('Existing evidence: do not repeat or overwrite.')
start=monotonic();record={'status':'started','budget_seconds':seconds,'started_unix':time(),
 'command_relative':['python','-B',script.name],'source_sha256_before':sha256(script.read_bytes()).hexdigest()}
meta.write_text(json.dumps(record,indent=2)+'\n')
with log.open('w') as stream:
 p=subprocess.Popen([sys.executable,'-B',str(script)],cwd=root,stdout=stream,stderr=subprocess.STDOUT,start_new_session=True)
 record['pid']=p.pid;meta.write_text(json.dumps(record,indent=2)+'\n')
 try:
  code=p.wait(timeout=max(.001,seconds-(monotonic()-start)))
  record.update(status='terminal',exit_code=code)
 except subprocess.TimeoutExpired:
  os.killpg(p.pid,signal.SIGKILL);p.wait()
  record.update(status='timed_out',exit_code=p.returncode)
try:os.killpg(p.pid,signal.SIGKILL)
except ProcessLookupError:pass
record.update(elapsed_seconds=monotonic()-start,log_sha256=sha256(log.read_bytes()).hexdigest(),
 source_sha256_after=sha256(script.read_bytes()).hexdigest())
meta.write_text(json.dumps(record,indent=2)+'\n');print(json.dumps(record))
