"""All common-radius events at the frozen receiver-face center bank; exact gate."""
from pathlib import Path
from fractions import Fraction as F
import hashlib,json,time
import round011_geometry as geo
ROOT=Path(__file__).resolve().parent

def floor(q):return q.numerator//q.denominator

def ceil(q):return -((-q.numerator)//q.denominator)

def event_key(point,r):
 los=[];his=[];nums=[]
 for b in range(3):
  lo=[ceil((point[j]-r)/geo.P[j]-geo.PH[b][j]) for j in range(3)];hi=[floor((point[j]+r)/geo.P[j]-geo.PH[b][j]) for j in range(3)];n=1
  for l,h in zip(lo,hi):n*=max(0,h-l+1)
  if n:assert all(-3<=l<=h<=3 for l,h in zip(lo,hi))
  else:lo=[0,0,0];hi=[-1,-1,-1]
  los.extend(lo);his.extend(hi);nums.append(n)
 return tuple(los+his),nums,sum((nums[b]*geo.WTS[b] for b in range(3)),F(0))

def main():
 start=time.monotonic();bankfile=ROOT/'round011_centers.json';bank=json.loads(bankfile.read_text());assert bank['center_count']<=128;assert bank['seed_input_sha256']==hashlib.sha256((ROOT/'round008_certificate_input.json').read_bytes()).hexdigest();assert bank['generator_sha256']==hashlib.sha256((ROOT/'round011_centers.py').read_bytes()).hexdigest();assert bank['geometry_helper_sha256']==hashlib.sha256((ROOT/'round011_geometry.py').read_bytes()).hexdigest();tail=geo.W
 for p in geo.P:tail*=1/p+F(1,6)
 assert tail<1 and min(geo.P)>1
 oldkeys=set(map(tuple,geo.SEED['keys']));eventkeys=set();percenter=[];all_events=positive=passing=tied=0
 for center in bank['centers']:
  point=tuple(map(F,center['point']));assert all(0<x<geo.P[j] for j,x in enumerate(point));radii=sorted({abs(geo.P[j]*(F(k)+geo.PH[b][j])-point[j]) for b in range(3) for j in range(3) for k in range(-3,4) if abs(geo.P[j]*(F(k)+geo.PH[b][j])-point[j])<3});assert len(radii)<=63;new_before=len(eventkeys);npos=npass=ntie=0;best=None;zero_mass_events=0
  for r in radii:
   key,nums,mass=event_key(point,r)
   if mass==0:zero_mass_events+=1;continue
   eventkeys.add(key);npos+=1;delta=mass-8*r**3
   if delta>0:npass+=1
   elif delta==0:ntie+=1
   ratio=None if r==0 else mass/(8*r**3)
   if best is None or r==0 or (best['r']!=0 and ratio>best['ratio']):best={'r':r,'mass':mass,'counts':nums,'key':key,'ratio':ratio}
  all_events+=len(radii);positive+=npos;passing+=npass;tied+=ntie
  percenter.append({'center_id':center['center_id'],'candidate_id':center['candidate_id'],'events':len(radii),'positive_mass_events':npos,'zero_mass_events':zero_mass_events,'strict_passing_events':npass,'tie_events':ntie,'new_distinct_event_keys':len(eventkeys)-new_before,'maximum_density_event':{'radius':str(best['r']),'mass':str(best['mass']),'orbit_counts':best['counts'],'key':list(best['key']),'density':'+infinity' if best['ratio'] is None else str(best['ratio']),'density_display':None if best['ratio'] is None else float(best['ratio'])} if best else None})
 allnew=sorted(eventkeys-oldkeys);nonempty=[];empty=0;bygeometry={};oldboxes,oldids=geo.old_geometry();oldgeoset=set(oldboxes);identical=0
 for key in allnew:
  box=geo.receiver_box(key)
  if box is None:empty+=1;continue
  nonempty.append(key)
  if box in oldgeoset:identical+=1
  else:bygeometry.setdefault(box,key)
 newfile=ROOT/'round011_new_keys.json';newfile.write_text(json.dumps([list(k) for k in nonempty],separators=(',',':'))+'\n');gate=geo.exact_coverage_test(oldboxes,list(bygeometry));gate['all_old639_retained']=True;gate['all_added_nonempty_boxes_contained_in_old_union']=not gate['positive_new_cell'];gate['finite_unions_exactly_equal_if_no_positive_cell']=not gate['positive_new_cell']
 output={'status':'positive new coverage; candidate certification/optional refinement permitted' if gate['positive_new_cell'] else 'no new finite coverage; stop with no optimizer','seed_input_sha256':bank['seed_input_sha256'],'center_bank_sha256':hashlib.sha256(bankfile.read_bytes()).hexdigest(),'new_keys_sha256':hashlib.sha256(newfile.read_bytes()).hexdigest(),'centers':len(bank['centers']),'total_radius_events':all_events,'positive_mass_events':positive,'strict_passing_events':passing,'tie_events':tied,'distinct_event_keys':len(eventkeys),'new_event_keys_before_geometry':len(allnew),'new_proved_empty_boxes':empty,'new_nonempty_keys':len(nonempty),'new_boxes_identical_to_old_geometry':identical,'distinct_new_box_geometries':len(bygeometry),'full_exact_uniform_tail_bound':str(tail),'tail_radius':3,'event_index_range':[-3,3],'all_scales_complete_at_selected_centers':True,'exact_coverage_gate':gate,'optimization_evaluations':0,'per_center':percenter,'elapsed_seconds':time.monotonic()-start};(ROOT/'round011_capture_result.json').write_text(json.dumps(output,indent=2)+'\n');print({k:v for k,v in output.items() if k!='per_center'})
if __name__=='__main__':main()
