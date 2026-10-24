"""One <=10s standard-library audit of four saved AK balls; zero new AK calls."""
from fractions import Fraction as Q
from hashlib import sha256
import json
from pathlib import Path
from time import monotonic
ROOT=Path(__file__).resolve().parent
OUT=ROOT/'saved_verification.json'
SOURCE_SHA='53195f55dc40dad9da6749e2bea7ad55d748cd9d3597fdea53d38f36948dad3a'
WRAPPER_SHA='5bb9140fc42aac834c134209f4a20e5d8e13602fdf82ce471d1ae643e11fb6a6'
if OUT.exists():raise SystemExit('No repeated audit or overwrite.')
start=monotonic()
def guard():
 if monotonic()-start>=10:raise TimeoutError('10s saved-value budget; no extension.')
def read(p):guard();return json.loads(p.read_text())
def digest(p):guard();return sha256(p.read_bytes()).hexdigest()
j,run,points=(read(ROOT/n) for n in ('result.json','diagnostic.run.json','points.json'))
assert run['status']=='terminal' and run['exit_code']==0 and 0<=run['elapsed_seconds']<60
assert digest(ROOT/'run_once.py')==WRAPPER_SHA
assert run['source_sha256_before']==run['source_sha256_after']==j['source_sha256_before']==j['source_sha256_after']==digest(ROOT/'diagnose_once.py')==SOURCE_SHA
assert run['log_sha256']==digest(ROOT/'diagnostic.log')
assert j['source_unchanged'] and j['inputs_sources_unchanged'] and j['input_source_sha256_before']==j['input_source_sha256_after']
for relative,wanted in j['input_source_sha256_before'].items():assert digest(ROOT/relative)==wanted
assert j['status']=='four_fixed_AK_calls_closed' and j['counts']==dict(started=4,completed=4,failed=0,interrupted=0,not_started=0)
assert j['maximum_calls']==4 and j['precision_bits']==192 and j['budget_seconds']==60
assert all(j[k]==0 for k in ('new_K_calls','new_box_models','LP_calls','full_checker_calls'))
assert digest(ROOT/'points.json')==j['points_sha256'] and points['points']==j['points']
assert j['trace_sha256']==digest(ROOT/'calls.jsonl')
events=[json.loads(line) for line in (ROOT/'calls.jsonl').read_text().splitlines()]
assert len(events)==8 and len(j['rows'])==4
intervals={};facts=[]
for index,row in enumerate(j['rows']):
 guard();a,b=events[2*index:2*index+2]
 assert a['call']==b['call']==row['call']==index+1 and a['status']=='started' and b==row and row['status']=='completed'
 assert all(row[k]==v for k,v in a.items() if k!='status')
 assert row['point_id']==('P1' if index<2 else 'P2') and row['quantity']=='AK' and row['threshold']=='0'
 assert row['method']==('corrected_first_exact_B_adapter' if index%2==0 else 'independent_core_point_methods')
 u,v=Q(row['u']),Q(row['v']);assert 0<v<u and u+v<Q(33,20)
 assert (u,v)==((Q(23497,65536),Q(10819,65536)) if index<2 else (Q(47105,131072),Q(21905,131072)))
 lo,hi=Q(row['lower']),Q(row['upper']);assert lo<=hi
 assert all(x.denominator&(x.denominator-1)==0 for x in (lo,hi))
 if index<2:
  assert hi<-Q(1,2500) and row['classification']=='strict_negative_AK' and row['strict_negative_neighborhood_justified'] is True
 else:
  assert lo>Q(1,250000) and row['classification']=='strict_positive_AK' and row['strict_negative_neighborhood_justified'] is False
 for value in row['parts'].values():assert Q(value['lower'])<=Q(value['upper'])
 intervals.setdefault(row['point_id'],[]).append((lo,hi));facts.append(dict(point=row['point_id'],method=row['method'],lower=str(lo),upper=str(hi)))
for cp in j['cross_implementation_comparisons']:
 guard();a,b=intervals[cp['point_id']]
 assert max(a[0],b[0])<=min(a[1],b[1])
 assert cp['intervals_overlap'] and cp['difference_contains_zero'] and Q(cp['difference']['lower'])<=0<=Q(cp['difference']['upper'])
assert len(j['cross_implementation_comparisons'])==2
guard()
answer=dict(status='passed_saved_four_ball_audit',budget_seconds=10,scalar_balls_checked=4,trace_events_checked=8,
 P1_upper_less_than='-1/2500',P2_lower_greater_than='1/250000',both_implementation_pairs_overlap=True,
 new_AK_calls=0,new_integrals=0,new_LP=0,scope='Saved Fraction endpoints, trace/counts and input/source hashes only; no kernel recomputation.',
 facts=facts,result_sha256=digest(ROOT/'result.json'),run_sha256=digest(ROOT/'diagnostic.run.json'),trace_sha256=j['trace_sha256'],
 verifier_sha256=digest(Path(__file__)),elapsed_seconds=monotonic()-start)
OUT.write_text(json.dumps(answer,indent=2)+'\n');print(json.dumps(answer))
