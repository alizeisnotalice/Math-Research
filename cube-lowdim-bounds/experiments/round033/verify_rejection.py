"""Standard-library exact check of one archived rejection witness, no new model."""
from pathlib import Path
from fractions import Fraction as Q
from hashlib import sha256
import json,time
r=Path(__file__).resolve().parent;out=r/'rejection_verification.json'
if out.exists():raise SystemExit('Existing one-row diagnosis; do not repeat.')
start=time.monotonic();raw=(r/'rejection_witness.json').read_bytes();d=json.loads(raw)
row=d['row'];q=list(map(Q,d['q']));assert len(q)==26
assert q[0]>0 and q[1]>=0 and q[14:]==list(map(abs,q[2:14]))
ij=row['bernstein_index'];pos=4*ij[0]+ij[1]
value=q[0]*Q(row['c_column'])+q[1]*Q(row['s_column'])
for j,model in enumerate(d['models']):
 assert model['direction']==d['direction_names'][j] and model['leaf_node_id']==row['leaf_node_id']
 assert model['parent']==d['leaf']['parent'] and model['local']==d['leaf']['local']
 b=model['bernstein'][pos];assert b['index']==ij
 lo,hi=Q(b['lower']),Q(b['upper']);e=Q(model['remainder_upper']);assert lo<=hi and e>=0
 mid=(lo+hi)/2;penalty=(hi-lo)/2+e
 assert mid==Q(row['t_midpoints'][j]) and penalty==Q(row['t_penalties'][j])
 value+=q[2+j]*mid-q[14+j]*penalty
half=Q(row['baseline_slack'])/2;residual=value-half
assert half>0 and residual==Q(d['expected_half_margin_residual'])<0 and value>0
record={'status':'one_half_margin_rejection_verified_but_original_row_positive','half_margin':str(half),'original_safe_row':str(value),'strengthened_residual':str(residual),'scope':'Only this row. Other rows and global positivity are not recalculated here.','witness_sha256':sha256(raw).hexdigest(),'source_sha256':sha256(Path(__file__).read_bytes()).hexdigest(),'elapsed_seconds':time.monotonic()-start}
assert record['elapsed_seconds']<10
out.write_text(json.dumps(record,indent=2)+'\n');print(json.dumps(record,indent=2))
