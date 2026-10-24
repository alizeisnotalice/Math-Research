"""Independent rational slab coverage and row verification; no Arb/model calls."""
from pathlib import Path
from fractions import Fraction as Q
from hashlib import sha256
import json,time,sys
ROOT=Path(__file__).resolve().parent
OUT=ROOT/'cover_verification.json'
if OUT.exists():raise SystemExit('Existing verification; no automatic repetition.')
start=time.monotonic()
r={'status':'started','scope':'Exact geometry and rational-row arithmetic only; no independent transcendental model enclosure.','leaves_checked':0,'rows_checked':0,'slabs_checked':0}
def save():OUT.write_text(json.dumps(r,indent=2)+'\n')
def guard():
 if time.monotonic()-start>=90:raise TimeoutError('90-second verification budget.')
save()
try:
 raw=(ROOT/'leaf_cover.json').read_bytes();d=json.loads(raw)
 assert d['complete'] is True
 assert not d['pending'] and not d['unresolved'] and d['in_flight'] is None
 assert (Q(d['alpha']),Q(d['R']),Q(d['h']))==(Q(47,40),Q(33,20),Q(1,32))
 braw=(ROOT.parent/'round01/scaled_certificate.json').read_bytes()
 bsha=sha256(braw).hexdigest()
 assert bsha=='21d39978590d3659021c093b3100e4e4d5e8dafe4834d46aa55d6c86c6a5e172'
 assert d['baseline_sha256']==bsha
 c1=Q(json.loads(braw)['c'])
 leaves=d['leaves'];assert 1<=len(leaves)<=1500
 parents={(7,5):[],(8,6):[],(10,6):[]}
 r.update(cover_sha256=sha256(raw).hexdigest(),baseline_sha256=bsha,parent_results=[])
 for leaf in leaves:
  guard();key=tuple(leaf['parent']);assert key in parents
  x0,x1,y0,y1=map(Q,leaf['local'])
  assert 0<=x0<x1<=1 and 0<=y0<y1<=1
  h=Q(1,32);p0,q0=key
  assert list(map(Q,leaf['physical']))==[h*(p0+x0),h*(p0+x1),h*(q0+y0),h*(q0+y1)]
  model_path=ROOT/leaf['model_file'];assert model_path.is_file()
  mraw=model_path.read_bytes();assert sha256(mraw).hexdigest()==leaf['model_sha256']
  model=json.loads(mraw)
  assert model['parent']==leaf['parent'] and model['local']==leaf['local'] and model['mode']==leaf['mode']
  assert model['rows']==leaf['rows']
  assert leaf['mode'] in ('const','Taylor')
  rows=leaf['rows'];assert len(rows)==16
  indices=set()
  for row in rows:
   guard();ij=tuple(row['bernstein_index']);assert ij not in indices and ij in {(i,j) for i in range(4) for j in range(4)}
   indices.add(ij)
   pos=ij[0]*4+ij[1]
   ss=model['AS'];assert Q(ss['remainder_upper'])>=0
   assert Q(ss['bernstein'][pos]['lower'])<=Q(ss['bernstein'][pos]['upper'])
   assert tuple(ss['bernstein'][pos]['index'])==ij
   assert Q(row['s_column'])==Q(ss['bernstein'][pos]['lower'])-Q(ss['remainder_upper'])
   if leaf['mode']=='const': assert Q(row['c_column'])==Q(model['const']['lower'])
   else:
    pp=model['Taylor'];assert Q(pp['remainder_upper'])>=0
    assert Q(pp['bernstein'][pos]['lower'])<=Q(pp['bernstein'][pos]['upper'])
    assert tuple(pp['bernstein'][pos]['index'])==ij
    assert Q(row['c_column'])==Q(pp['bernstein'][pos]['lower'])-Q(pp['remainder_upper'])
   val=c1*Q(row['c_column'])+Q(row['s_column'])
   assert val==Q(row['baseline_slack']) and val>0
   r['rows_checked']+=1
  parents[key].append((x0,x1,y0,y1));r['leaves_checked']+=1
 for key,rects in parents.items():
  guard();assert rects
  xs=sorted({Q(0),Q(1)}|{a for box in rects for a in box[:2]})
  area=sum((b-a)*(d-c) for a,b,c,d in rects)
  assert area==1
  for left,right in zip(xs,xs[1:]):
   guard();assert left<right
   intervals=sorted((c,d) for a,b,c,d in rects if a<=left and right<=b)
   cursor=Q(0)
   for bottom,top in intervals:
    assert bottom==cursor,'Gap or overlapping interiors in rational slab.'
    cursor=top
   assert cursor==1,'Missing vertical coverage.'
   r['slabs_checked']+=1
  r['parent_results'].append({'parent':list(key),'leaves':len(rects),'local_area':str(area),'slabs':len(xs)-1,
   'coverage':'Every rational x slab has exact disjoint y partition [0,1].'})
  save()
 assert r['rows_checked']<=24000
 r.update(status='exact_slab_cover_and_all_positive_rows_verified',physical_area=str(Q(3,1024)),
  mathematical_scope='Complete local cover of three parent knot cells and strict baseline Bernstein rows; not full-domain certificate or new C2 bound.')
except Exception as e:r.update(status='failed_or_incomplete',reason=f'{type(e).__name__}: {e}')
r.update(wall_seconds=time.monotonic()-start,source_sha256=sha256(Path(__file__).read_bytes()).hexdigest(),python=sys.version)
save();print(json.dumps({k:v for k,v in r.items() if k not in ['parent_results','python']},indent=2))
