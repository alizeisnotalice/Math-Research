"""One exact saved-model/partial-tree audit; no Arb or kernel evaluation."""
from pathlib import Path
from fractions import Fraction as Q
from hashlib import sha256
import json,time
r=Path(__file__).resolve().parent;out=r/'saved_verification.json'
if out.exists():raise SystemExit('Existing audit; do not repeat.')
start=time.monotonic()
def guard():
 if time.monotonic()-start>90:raise TimeoutError('90s exact saved-audit budget')
def read(p):return json.loads(p.read_text())
d=read(r/'leaf_cover.json');j=read(r/'result.json');run=read(r/'export.run.json')
assert sha256((r/'leaf_cover.json').read_bytes()).hexdigest()==j['cover_sha256']
assert run['source_sha256_before']==run['source_sha256_after']==j['source_sha256']
assert sha256((r/'export_cover.py').read_bytes()).hexdigest()==j['source_sha256']
b=r.parent/'round01/scaled_certificate.json';assert sha256(b.read_bytes()).hexdigest()==d['baseline_sha256']=='21d39978590d3659021c093b3100e4e4d5e8dafe4834d46aa55d6c86c6a5e172'
c1=Q(read(b)['c']);assert d['parents']==[[52,1]] and len(d['roots'])==324
nodes=d['nodes'];assert all(n['id']==i for i,n in enumerate(nodes))
models=0;rows=0;worst=None;accepted=[]
for n in nodes:
 guard();assert n['parent']==[52,1]
 x0,x1,y0,y1=map(Q,n['local']);assert 0<=x0<x1<=1 and 0<=y0<y1<=1
 physical=[(52+x0)/32,(52+x1)/32,(1+y0)/32,(1+y1)/32]
 assert physical==list(map(Q,n['physical'])) and physical[0]+physical[2]>Q(33,20)
 if 'children' in n:
  xm,ym=(x0+x1)/2,(y0+y1)/2
  boxes=[(x0,xm,y0,ym),(xm,x1,y0,ym),(x0,xm,ym,y1),(xm,x1,ym,y1)]
  assert len(n['children'])==4
  for child,box in zip(n['children'],boxes):
   z=nodes[child];assert tuple(map(Q,z['local']))==box and z['parent_id']==n['id'] and z['depth']==n['depth']+1
 if 'model_file' not in n:continue
 p=r/n['model_file'];assert sha256(p.read_bytes()).hexdigest()==n['model_sha256'];m=read(p)
 assert m['mode']=='const' and m['node_id']==n['id'] and m['local']==n['local']
 assert len(m['rows'])==16;E=Q(m['AS']['remainder_upper']);assert E>=0
 assert Q(m['const']['lower'])<=Q(m['const']['upper'])
 slacks=[]
 for k,row in enumerate(m['rows']):
  guard();i,jj=divmod(k,4);assert row['bernstein_index']==[i,jj]
  bern=m['AS']['bernstein'][k];assert bern['index']==[i,jj]
  lo,hi=Q(bern['lower']),Q(bern['upper']);assert lo<=hi
  a=Q(m['const']['lower']);bb=lo-E;slack=c1*a+bb
  assert (Q(row['c_column']),Q(row['s_column']),Q(row['baseline_slack']))==(a,bb,slack)
  slacks.append(slack);rows+=1
  if worst is None or slack<Q(worst['slack']):worst=dict(node=n['id'],index=[i,jj],slack=str(slack),before_spline_remainder=str(c1*a+lo),spline_remainder=str(E),c_column=str(a),AS_bern_lower=str(lo))
 assert m['accepted']==all(s>0 for s in slacks)
 assert (n['status']=='accepted')==m['accepted']
 if m['accepted']:accepted.append(n['id'])
 models+=1
assert models==j['counts']['boxes_completed'] and rows==16*models
assert sorted(accepted)==sorted(leaf['node_id'] for leaf in d['leaves'])
for leaf in d['leaves']:
 n=nodes[leaf['node_id']];m=read(r/n['model_file'])
 assert leaf['rows']==m['rows'] and leaf['local']==n['local']
rects=[tuple(map(Q,nodes[i]['local'])) for i in accepted]
area=sum(((b-a)*(dd-c) for a,b,c,dd in rects),Q(0));slabs=0
if d['complete']:
 assert not d['pending'] and not d['unresolved'] and d['in_flight'] is None
 xs=sorted({Q(0),Q(1)}|{v for box in rects for v in box[:2]})
 for a,b in zip(xs,xs[1:]):
  guard();cursor=Q(0)
  for c,dd in sorted((c,dd) for x0,x1,c,dd in rects if x0<=a and b<=x1):assert c==cursor;cursor=dd
  assert cursor==1;slabs+=1
 assert area==1
result=dict(status='passed_saved_arithmetic_audit',scope='Exact rows/tree from saved interval enclosures; not independent transcendental evaluation.',models_checked=models,rows_checked=rows,accepted_leaves=len(accepted),local_accepted_area=str(area),complete=d['complete'],slabs_checked=slabs,worst_saved_row=worst,elapsed_seconds=time.monotonic()-start,cover_sha256=j['cover_sha256'],result_sha256=sha256((r/'result.json').read_bytes()).hexdigest())
out.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result))
