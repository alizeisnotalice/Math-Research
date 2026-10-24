"""One 120-second process-group bounded discovery run; refuse replacement."""
import sys,subprocess,time,os,signal,json,hashlib
from pathlib import Path
root=Path(__file__).resolve().parent
meta=root/'discovery.run.json';log=root/'discovery.log'
if meta.exists() or log.exists():raise SystemExit('Existing run: inspect; do not repeat.')
cmd=[sys.executable,'-B',str(root/'search.py')]
t=time.monotonic()
source_hash=hashlib.sha256((root/'search.py').read_bytes()).hexdigest()
record={'search_sha256_before':source_hash,'command_relative':['python','-B','search.py'],'budget_seconds':120,'started_unix':time.time(),'status':'running'}
meta.write_text(json.dumps(record,indent=2)+'\n')
with log.open('w') as out:
 p=subprocess.Popen(cmd,cwd=root,stdout=out,stderr=subprocess.STDOUT,start_new_session=True)
 try:
  code=p.wait(timeout=max(0,120-(time.monotonic()-t)));record.update(status='terminal',exit_code=code)
 except subprocess.TimeoutExpired:
  os.killpg(p.pid,signal.SIGKILL);p.wait();record.update(status='hard_timeout',exit_code=p.returncode)
record.update(elapsed_seconds=time.monotonic()-t,log_sha256=hashlib.sha256(log.read_bytes()).hexdigest(),search_sha256=hashlib.sha256((root/'search.py').read_bytes()).hexdigest())
record['source_unchanged']=record['search_sha256']==source_hash
meta.write_text(json.dumps(record,indent=2)+'\n');print(json.dumps(record))
