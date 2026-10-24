"""Exact decimal-ball sign/overlap audit only; no kernel/Arb evaluation."""
from pathlib import Path
from fractions import Fraction as Q
from hashlib import sha256
import json,re
root=Path(__file__).resolve().parent
out=root/'saved_verification.json'
if out.exists():raise SystemExit('Existing audit: do not overwrite.')
def ball(s):
 m=re.fullmatch(r'\[([^\[\]]+) \+/- ([^\[\]]+)\]',s)
 assert m,s
 c,r=map(Q,m.groups());assert r>=0
 return c-r,c+r
j=json.loads((root/'result.json').read_text());run=json.loads((root/'diagnose.run.json').read_text())
assert j['status']=='four_fixed_point_calls_closed'
assert j['counts']==dict(started=4,completed=4,failed=0,interrupted=0,not_started=0)
assert run['status']=='terminal' and run['exit_code']==0 and run['elapsed_seconds']<=60
assert run['source_sha256_before']==run['source_sha256_after']==j['source_sha256']==sha256((root/'diagnose_once.py').read_bytes()).hexdigest()
assert j['trace_sha256']==sha256((root/'calls.jsonl').read_bytes()).hexdigest()
assert j['points_sha256']==sha256((root/'points.json').read_bytes()).hexdigest()
trace=[json.loads(s) for s in (root/'calls.jsonl').read_text().splitlines()]
assert len(trace)==8
rows=[];intervals=[]
for i,row in enumerate(j['rows']):
 assert row['call']==i+1 and row['status']=='completed'
 assert trace[2*i]['status']=='started' and trace[2*i+1]==row
 lo,hi=ball(row['ball']);assert hi<0
 threshold=-Q(1,10**6 if row['quantity']=='K' else 10**4)
 assert hi<threshold
 rows.append(dict(call=i+1,lower=str(lo),upper=str(hi),strict_upper_bound=str(threshold)))
 intervals.append((lo,hi))
for a,b in ((0,1),(2,3)):
 assert max(intervals[a][0],intervals[b][0])<=min(intervals[a][1],intervals[b][1])
result=dict(status='passed',scope='Only exact Fraction parsing of saved decimal balls, trace and provenance. No new point, integral, or global evaluation.',rows=rows,overlap_pairs=[[1,2],[3,4]],result_sha256=sha256((root/'result.json').read_bytes()).hexdigest(),run_sha256=sha256((root/'diagnose.run.json').read_bytes()).hexdigest())
out.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(dict(status='passed',saved_balls=4,kernel_calls=0)))
