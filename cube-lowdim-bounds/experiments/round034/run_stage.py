"""Unique finite repair or global certificate stages, owned process groups only."""
from pathlib import Path
from hashlib import sha256
import json,os,signal,subprocess,sys,time
from fractions import Fraction as Q
root=Path(__file__).resolve().parent;stage=sys.argv[1]
assert stage in ('repair','first','second')
seconds=90 if stage=='repair' else 180
meta=root/(stage+'.run.json');log=root/(stage+'.log')
if meta.exists() or log.exists():raise SystemExit('Existing stage evidence; no retry.')
public=root.parent/'upstream-v4'
if stage=='repair':
 script=root/'repair_once.py';cwd=root;command=[sys.executable,'-B',str(script)]
 assert not (root/'result.json').exists()
else:
 result=json.loads((root/'result.json').read_text())
 assert result['candidate_exported'] is True
 candidate=root/'candidate.json';assert sha256(candidate.read_bytes()).hexdigest()==result['candidate_sha256']
 spec=json.loads(candidate.read_text())
 assert set(spec)=={'alpha','R','h','N','c','a'}
 cost_path=root/'exact_cost.json';assert sha256(cost_path.read_bytes()).hexdigest()==result['exact_cost_sha256']
 cost=json.loads(cost_path.read_text())
 structured_path=root/'structured_kernel.json';assert sha256(structured_path.read_bytes()).hexdigest()==result['structured_kernel_sha256']
 structured=json.loads(structured_path.read_text())
 assert all(spec[k]==structured[k] for k in spec)
 assert cost['certificate_sha256']==result['structured_kernel_sha256'] and cost['plain_candidate_sha256']==result['candidate_sha256']
 assert Q(cost['upper'])<Q(35909,10000) and Q(cost['upper'])<Q(cost['baseline_lower'])
 assert cost['improves_baseline'] is True
 pins={'checker/cert_run.py':'079906db354febdd57d537f0d747428bc90bd37fa88297449cd1d46c7446310d','checker/hlcert.py':'4bd7e823bdcb336f8e20c188d68d6d8c1e4aa22ec52353885059c24445317db5','checker/independent/hl_check.py':'0a5ee85e9b6c250494aaf88abce55bb28a4ba5fc46075c531eeefaa0cc2888f8','checker/independent/hl_core.py':'c3d864b1368ad8150fd3366bf77f687a9a3cdd49b924a20d576aa678aed42ebf'}
 for name,wanted in pins.items():assert sha256((public/name).read_bytes()).hexdigest()==wanted,name
 script=public/('checker/cert_run.py' if stage=='first' else 'checker/independent/hl_check.py')
 output=root/(stage+'_report.json');assert not output.exists();cwd=public
 args=['--Q','8','--maxdepth','7','--Kdepth','6','--procs','6','--only','MKA'] if stage=='first' else ['--claim','3.5909','--procs','6']
 command=[sys.executable,'-I',str(script),str(candidate),*args,'--out',str(output)]
start=time.monotonic();record={'status':'started','stage':stage,'budget_seconds':seconds,'started_unix':time.time(),'command':command,'source_sha256_before':sha256(script.read_bytes()).hexdigest()}
if stage!='repair':record.update(candidate_sha256=sha256(candidate.read_bytes()).hexdigest(),dependency_hashes=pins)
meta.write_text(json.dumps(record,indent=2)+'\n')
with log.open('x') as stream:
 p=subprocess.Popen(command,cwd=cwd,stdout=stream,stderr=subprocess.STDOUT,start_new_session=True)
 record['pid']=p.pid;meta.write_text(json.dumps(record,indent=2)+'\n')
 try:
  code=p.wait(timeout=max(.001,seconds-(time.monotonic()-start)));record.update(status='terminal',exit_code=code)
 except subprocess.TimeoutExpired:
  os.killpg(p.pid,signal.SIGKILL);p.wait();record.update(status='timed_out',exit_code=p.returncode)
try:os.killpg(p.pid,signal.SIGKILL)
except ProcessLookupError:pass
record.update(elapsed_seconds=time.monotonic()-start,source_sha256_after=sha256(script.read_bytes()).hexdigest(),log_sha256=sha256(log.read_bytes()).hexdigest())
if stage!='repair':
 record['candidate_sha256_after']=sha256(candidate.read_bytes()).hexdigest()
 assert record['candidate_sha256_after']==record['candidate_sha256']
 if output.exists():record['report_sha256']=sha256(output.read_bytes()).hexdigest()
meta.write_text(json.dumps(record,indent=2)+'\n');print(json.dumps(record))
