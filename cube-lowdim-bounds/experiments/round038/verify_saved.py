"""One exact audit of saved Bernstein combinations and call ledger; no Arb."""
from pathlib import Path
from fractions import Fraction as Q
from hashlib import sha256
import json,time
r=Path(__file__).resolve().parent;out=r/'saved_verification.json'
if out.exists():raise SystemExit('Existing verification; no repetition.')
t=time.monotonic()
def read(p):return json.loads(p.read_text())
j=read(r/'result.json');m=read(r/'new_model.json');run=read(r/'compare.run.json')
assert j['status']=='completed' and j['accepted'] is True and j['in_flight']==[]
assert run['status']=='terminal' and run['exit_code']==0 and run['elapsed_seconds']<=60
assert run['source_sha256_before']==run['source_sha256_after']==j['source_sha256']==sha256((r/'compare_once.py').read_bytes()).hexdigest()
assert sha256((r/'new_model.json').read_bytes()).hexdigest()==j['new_model_sha256']
assert sha256((r/'calls.jsonl').read_bytes()).hexdigest()==j['calls_sha256']
old_path=r.parent/'round18/models.local/node_00615.json';old=read(old_path)
assert sha256(old_path.read_bytes()).hexdigest()==m['source_AS_model_sha256']
b=r.parent/'round01/scaled_certificate.json';assert sha256(b.read_bytes()).hexdigest()==m['baseline_sha256']
c1=Q(read(b)['c']);assert c1==Q(j['baseline_c'])
assert m['physical']==old['physical'] and m['local']==old['local']
e=Q(m['Psi_remainder_upper']);es=Q(old['AS']['remainder_upper'])
assert 0<=Q(m['Psi_remainder_lower'])<=e and es==Q(m['AS_remainder_upper'])
slacks=[];diffs=[]
for k,row in enumerate(m['rows']):
 assert len(m['rows'])==len(m['Psi_bernstein'])==16
 i,jj=divmod(k,4);B=m['Psi_bernstein'][k];S=old['AS']['bernstein'][k]
 assert B['index']==S['index']==row['bernstein_index']==[i,jj]
 assert Q(B['lower'])<=Q(B['upper'])
 a=Q(B['lower'])-e;bb=Q(S['lower'])-es;s=c1*a+bb;old_s=Q(old['rows'][k]['baseline_slack'])
 assert (Q(row['c_column']),Q(row['s_column']),Q(row['baseline_slack']),Q(row['old_const_slack']),Q(row['difference']))==(a,bb,s,old_s,s-old_s)
 assert s>Q(7,100)
 slacks.append(s);diffs.append(s-old_s)
assert min(slacks)==Q(j['minimum_slack'])
assert j['minimum_index']==list(divmod(slacks.index(min(slacks)),4))
trace=[json.loads(x) for x in (r/'calls.jsonl').read_text().splitlines()]
assert len(j['calls'])==25 and len(trace)==50
for kind,n in [('Psi_model',1),('Hint',12),('acb_integral',12)]:
 assert j['counts'][kind]==dict(started=n,completed=n,failed=0,interrupted=0)
 assert sum(x['kind']==kind for x in j['calls'])==n
for call in j['calls']:
 events=[x for x in trace if x['id']==call['id']]
 assert len(events)==2 and events[0]['status']=='started' and events[1]==call and call['status']=='completed'
assert all(j[k]==0 for k in ['new_AS_calls','new_const_calls','point_calls','subdivisions','T_models','LP_calls','full_checker_calls'])
result=dict(status='passed',scope='Exact 16-row arithmetic and execution evidence only; no independent integration or model enclosure.',rows=16,all_strictly_above='7/100',minimum_slack=str(min(slacks)),old_minimum=str(min(Q(x['baseline_slack']) for x in old['rows'])),minimum_gain=str(min(slacks)-min(Q(x['baseline_slack']) for x in old['rows'])),rows_improved=sum(x>0 for x in diffs),rows_worsened=sum(x<0 for x in diffs),trace_events=50,elapsed_seconds=time.monotonic()-t,result_sha256=sha256((r/'result.json').read_bytes()).hexdigest(),new_model_sha256=j['new_model_sha256'])
out.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result))
