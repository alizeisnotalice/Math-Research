"""One <=10s saved physical-box knot classification; no kernel evaluation."""
from collections import Counter
from fractions import Fraction as Q
from hashlib import sha256
import json
from pathlib import Path
from time import monotonic
R=Path(__file__).resolve().parent
out=R/'first_failure_geometry_audit.json'
if out.exists():raise SystemExit('No repeat or overwrite.')
t=monotonic()
p=R/'first_report.json';raw=p.read_bytes()
assert sha256(raw).hexdigest()=='cde8d4c4547a706d420df92a0c776643aaa86d82a71a7144cf9a1d584b506b77'
fails=json.loads(raw)['A_fail'];counts=Counter();ambiguous=[]
for index,box in enumerate(fails):
 if monotonic()-t>=10:raise TimeoutError('10s saved geometry budget; no extension.')
 a,b,c,d=(32*Q(s) for s in box[:4])
 assert a<b and c<d
 floor=lambda x:x.numerator//x.denominator
 ceil=lambda x:-((-x.numerator)//x.denominator)
 lo=(floor(a),floor(c));hi=(ceil(b)-1,ceil(d)-1)
 if lo!=hi:ambiguous.append(index)
 counts[lo]+=1
answer=dict(status='passed_saved_physical_geometry_classification',budget_seconds=10,
 failures_classified=len(fails),parent_counts=[dict(parent=list(k),count=v) for k,v in sorted(counts.items())],
 boxes_crossing_knot_interior=ambiguous,all_in_parent_11_5=(not ambiguous and set(counts)=={(11,5)}),
 scope='Only exact Fraction classification of saved physical box endpoints; no new points or generator values.',
 source_report_sha256=sha256(raw).hexdigest(),script_sha256=sha256(Path(__file__).read_bytes()).hexdigest(),
 elapsed_seconds=monotonic()-t)
out.write_text(json.dumps(answer,indent=2)+'\n');print(json.dumps(answer))
