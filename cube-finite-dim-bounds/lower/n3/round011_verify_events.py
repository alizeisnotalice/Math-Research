"""Independent rational-cube emptiness certificate and periodic event closure."""
from pathlib import Path
from fractions import Fraction as F
from itertools import product
import hashlib,json,time
import round011_geometry as geo
ROOT=Path(__file__).resolve().parent

P=tuple(map(F,geo.SEED['periods']));PH=tuple(tuple(map(F,q)) for q in geo.SEED['phases']);WT=tuple(map(F,geo.SEED['weights']))
def floor(x):return x.numerator//x.denominator
def ceil(x):return -((-x.numerator)//x.denominator)

def data(key):
 lo=[key[3*b:3*b+3] for b in range(3)];hi=[key[9+3*b:12+3*b] for b in range(3)];nums=[]
 for b in range(3):
  n=1
  for j in range(3):n*=max(0,hi[b][j]-lo[b][j]+1)
  nums.append(n)
 active=[b for b in range(3) if nums[b]];mass=sum((nums[b]*WT[b] for b in active),F(0));assert mass>0
 LC=[max(P[j]*(hi[b][j]+PH[b][j]) for b in active) for j in range(3)];HC=[min(P[j]*(lo[b][j]+PH[b][j]) for b in active) for j in range(3)]
 return lo,hi,nums,mass,LC,HC

def failure(mass,LC,HC):
 out=[]
 for j in range(3):
  span=LC[j]-HC[j];assert span>=0
  if mass<=span**3:out.append((span**3-mass,j,0))
  upper=LC[j]-P[j]
  if upper>=0 and mass<=8*upper**3:out.append((8*upper**3-mass,j,1))
  lower=-HC[j]
  if lower>=0 and mass<=8*lower**3:out.append((8*lower**3-mass,j,2))
 return max(out) if out else None

def shifted(key,q):
 lo,hi,num,mass,LC,HC=data(key);a=[];b=[]
 for orbit in range(3):
  if num[orbit]:a.extend(lo[orbit][j]+q[j] for j in range(3));b.extend(hi[orbit][j]+q[j] for j in range(3))
  else:a.extend([0,0,0]);b.extend([-1,-1,-1])
 return tuple(a+b)

def main():
 start=time.monotonic();bankfile=ROOT/'round011_centers.json';bank=json.loads(bankfile.read_text());capture=json.loads((ROOT/'round011_capture_result.json').read_text());assert capture['center_bank_sha256']==hashlib.sha256(bankfile.read_bytes()).hexdigest();old=set(map(tuple,geo.SEED['keys']));events=set();event_count=0
 for row in bank['centers']:
  x=tuple(map(F,row['point']));radii={abs(P[j]*(k+PH[b][j])-x[j]) for b in range(3) for j in range(3) for k in range(-3,4) if abs(P[j]*(k+PH[b][j])-x[j])<3};event_count+=len(radii)
  for r in radii:
   lo=[];hi=[];mass=F(0)
   for b in range(3):
    l=[ceil((x[j]-r)/P[j]-PH[b][j]) for j in range(3)];h=[floor((x[j]+r)/P[j]-PH[b][j]) for j in range(3)];n=1
    for a,z in zip(l,h):n*=max(0,z-a+1)
    if not n:l=[0,0,0];h=[-1,-1,-1]
    lo.extend(l);hi.extend(h);mass+=n*WT[b]
   if mass>0:events.add(tuple(lo+hi))
 raw=sorted(events-old);assert event_count==capture['total_radius_events'] and len(events)==capture['distinct_event_keys'] and len(raw)==capture['new_event_keys_before_geometry'];fp=ROOT/'round011_new_event_keys.json';fp.write_text(json.dumps([list(k) for k in raw],separators=(',',':'))+'\n');witness=[];strict=zero=global_empty=global_positive=translated_tests=0;gaps=[];positive_translates=set();old_translates=set();translated_empty=0;global_degenerate=0;translated_equal=0
 for key in raw:
  lo,hi,num,mass,LC,HC=data(key);assert mass<216
  proof=failure(mass,LC,HC);assert proof is not None;gap,j,kind=proof;witness.append([j,kind]);gaps.append(gap)
  if gap>0:strict+=1
  else:zero+=1
  intrinsic=[(LC[j]-HC[j])**3-mass for j in range(3)]
  if max(intrinsic)>=0:
   global_empty+=1
   if max(intrinsic)==0:global_degenerate+=1
   continue
  global_positive+=1;ranges=[range(floor(-(HC[j]+3)/P[j])+1,ceil((P[j]-LC[j]+3)/P[j])) for j in range(3)];hits=0
  for q in product(*ranges):
   translated_tests+=1;LL=[LC[j]+q[j]*P[j] for j in range(3)];HH=[HC[j]+q[j]*P[j] for j in range(3)];proof2=failure(mass,LL,HH)
   if proof2 is not None:
    translated_empty+=1
    if proof2[0]==0:translated_equal+=1
    continue
   hits+=1;k2=shifted(key,q);assert all(-3<=z<=3 for z in k2);b2=geo.receiver_box(k2);assert b2 is not None
   if k2 in old:old_translates.add(k2)
   else:positive_translates.add(k2)
  assert hits>0
 new_boxes=[];seen=set()
 for k in sorted(positive_translates):
  b=geo.receiver_box(k)
  if b not in seen:new_boxes.append(b);seen.add(b)
 old_boxes,oldids=geo.old_geometry();gate=geo.exact_coverage_test(old_boxes,new_boxes);addedfile=ROOT/'round011_added_keys.json';addedfile.write_text(json.dumps([list(k) for k in sorted(positive_translates)],separators=(',',':'))+'\n')
 result={'status':'no new finite coverage after complete periodic closure; no optimizer' if not gate['positive_new_cell'] else 'positive new finite coverage; candidate gate open','seed_input_sha256':bank['seed_input_sha256'],'center_bank_sha256':capture['center_bank_sha256'],'raw_new_event_keys_sha256':hashlib.sha256(fp.read_bytes()).hexdigest(),'raw_new_event_key_count':len(raw),'rational_empty_witness_codes':{'0':'M <= (LC-HC)^3: intrinsic coordinate zero width','1':'M <=8*(LC-P)^3,LC>=P: above fundamental cell','2':'M <=8*(-HC)^3,HC<=0: below fundamental cell'},'witness_axis_type_in_raw_key_order':witness,'strict_raw_empty_certificates':strict,'zero_gap_raw_certificates':zero,'minimum_chosen_cubic_gap':str(min(gaps)),'globally_zero_volume_blocks':global_empty,'globally_degenerate_zero_width_blocks':global_degenerate,'globally_positive_but_raw_clipped_empty_blocks':global_positive,'periodic_shift_tests':translated_tests,'periodic_shift_proved_empty':translated_empty,'translated_zero_volume_equalities':translated_equal,'distinct_positive_translations_already_old':len(old_translates),'distinct_positive_translations_new':len(positive_translates),'new_translated_keys_sha256':hashlib.sha256(addedfile.read_bytes()).hexdigest(),'exact_coverage_gate':gate,'all_old639_retained':True,'finite_unions_exactly_equal':not gate['positive_new_cell'],'optimization_evaluations':0,'elapsed_seconds':time.monotonic()-start};(ROOT/'round011_verification.json').write_text(json.dumps(result,indent=2)+'\n');print({k:v for k,v in result.items() if k!='witness_axis_type_in_raw_key_order'})
if __name__=='__main__':main()
