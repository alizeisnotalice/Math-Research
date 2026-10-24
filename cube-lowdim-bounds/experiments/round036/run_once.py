"""One fixed baseline full check; 180s wall budget for the owned process group."""
from pathlib import Path
from hashlib import sha256
import json,os,signal,subprocess,sys,time
root=Path(__file__).resolve().parent
meta=root/'first.run.json';log=root/'first.log';output=root/'first_report.json'
if any(p.exists() for p in (meta,log,output)):raise SystemExit('Existing evidence; no repetition.')
manifest=json.loads((root/'patch_manifest.json').read_text())
# Exact pinned hashes are added by root after the independent static review.
pins=json.loads((root/'approved_inputs.json').read_text())
for name,wanted in pins.items():
 assert sha256((root/name).read_bytes()).hexdigest()==wanted,name
spec=root.parent/'round01/scaled_certificate.json'
assert sha256(spec.read_bytes()).hexdigest()=='21d39978590d3659021c093b3100e4e4d5e8dafe4834d46aa55d6c86c6a5e172'
script=root/'fixed_checker/cert_run.py'
command=[sys.executable,'-I',str(script),str(spec),'--Q','8','--maxdepth','7','--Kdepth','6','--procs','6','--only','MKA','--out',str(output)]
start=time.monotonic()
record=dict(status='started',budget_seconds=180,started_unix=time.time(),command=command,source_hashes_before=pins,certificate_sha256=sha256(spec.read_bytes()).hexdigest())
meta.write_text(json.dumps(record,indent=2)+'\n')
with log.open('x') as stream:
 p=subprocess.Popen(command,cwd=root,stdout=stream,stderr=subprocess.STDOUT,start_new_session=True)
 record['pid']=p.pid;meta.write_text(json.dumps(record,indent=2)+'\n')
 try:
  code=p.wait(timeout=max(.001,180-(time.monotonic()-start)));record.update(status='terminal',exit_code=code)
 except subprocess.TimeoutExpired:
  os.killpg(p.pid,signal.SIGKILL);p.wait();record.update(status='timed_out',exit_code=p.returncode)
try:os.killpg(p.pid,signal.SIGKILL)
except ProcessLookupError:pass
record.update(elapsed_seconds=time.monotonic()-start,log_sha256=sha256(log.read_bytes()).hexdigest(),source_hashes_after={name:sha256((root/name).read_bytes()).hexdigest() for name in pins},certificate_sha256_after=sha256(spec.read_bytes()).hexdigest())
if output.exists():record['report_sha256']=sha256(output.read_bytes()).hexdigest()
meta.write_text(json.dumps(record,indent=2)+'\n');print(json.dumps(record))
