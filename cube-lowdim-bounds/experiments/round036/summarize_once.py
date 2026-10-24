"""Read-only terminal/provenance summary. No mathematical evaluation."""
from pathlib import Path
from hashlib import sha256
import json
r=Path(__file__).resolve().parent
out=r/'result.json'
if out.exists():raise SystemExit('Existing summary; no overwrite.')
run=json.loads((r/'first.run.json').read_text());assert run['status'] in ('terminal','timed_out')
pins=json.loads((r/'approved_inputs.json').read_text())
assert pins==run['source_hashes_before']==run['source_hashes_after']
for n,h in pins.items():assert sha256((r/n).read_bytes()).hexdigest()==h
assert run['certificate_sha256']==run['certificate_sha256_after']=='21d39978590d3659021c093b3100e4e4d5e8dafe4834d46aa55d6c86c6a5e172'
log=(r/'first.log').read_bytes();assert sha256(log).hexdigest()==run['log_sha256']
j=dict(status='not_closed',scope='One corrected first checker on unchanged B1; no new bound, optimization or exact cost run.',elapsed_seconds=run['elapsed_seconds'],exit_code=run['exit_code'],run_status=run['status'],input_hashes_consistent=True,old_source_unchanged=True)
p=r/'first_report.json'
if p.exists():
 assert sha256(p.read_bytes()).hexdigest()==run['report_sha256']
 report=json.loads(p.read_text());j['report']=report
 j['status']='passed_corrected_first' if run['status']=='terminal' and run['exit_code']==0 and report.get('K_fail')==[] and report.get('A_fail')==[] and report.get('A_cells',0)>0 and 'PASS' in log.decode().splitlines() else 'not_closed'
original=r.parent/'upstream-v4/checker/hlcert.py'
if not original.exists():original=r/'original_checker/hlcert.py'
assert sha256(original.read_bytes()).hexdigest()=='4bd7e823bdcb336f8e20c188d68d6d8c1e4aa22ec52353885059c24445317db5'
j['run_sha256']=sha256((r/'first.run.json').read_bytes()).hexdigest();j['report_sha256']=run.get('report_sha256');j['log_sha256']=run['log_sha256']
out.write_text(json.dumps(j,indent=2)+'\n');print(json.dumps(j))
