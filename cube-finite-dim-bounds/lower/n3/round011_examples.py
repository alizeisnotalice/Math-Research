"""Stable fixed examples for main's serialized MATLAB interface cross-check."""
from pathlib import Path
from fractions import Fraction as F
import hashlib,json
import round011_geometry as geo
from round011_capture import event_key
from round011_verify_events import data
ROOT=Path(__file__).resolve().parent

def main():
 bank=json.loads((ROOT/'round011_centers.json').read_text());raw=json.loads((ROOT/'round011_new_event_keys.json').read_text());proof=json.loads((ROOT/'round011_verification.json').read_text());known=set(map(tuple,geo.SEED['keys']))|set(map(tuple,raw));probes=[]
 for centerid in [0,64]:
  row=bank['centers'][centerid];point=tuple(map(F,row['point']));entries={abs(geo.P[j]*(k+geo.PH[b][j])-point[j]) for b in range(3) for j in range(3) for k in range(-3,4)}
  for r in [F(1,4),F(1,3),F(1,2),F(2,3),F(3,4),F(1),F(5,4),F(3,2)]:
   key,counts,mass=event_key(point,r);margin=min(abs(r-e) for e in entries)
   if mass>0 and margin>F(1,10**6):
    assert key in known;probes.append({'center_id':centerid,'type':row['type'],'point':row['point'],'radius':str(r),'minimum_event_separation':str(margin),'orbit_counts':counts,'mass':str(mass),'key':list(key)});break
  else:raise AssertionError('No stable predeclared diagnostic radius')
 empty=[]
 for kind in range(3):
  for i,(axis,k) in enumerate(proof['witness_axis_type_in_raw_key_order']):
   if k!=kind:continue
   key=raw[i];lo,hi,counts,mass,LC,HC=data(key);rhs=(LC[axis]-HC[axis])**3 if kind==0 else 8*(LC[axis]-geo.P[axis])**3 if kind==1 else 8*(-HC[axis])**3;assert rhs>mass;empty.append({'raw_key_index':i,'key':key,'axis':axis,'kind':kind,'orbit_counts':counts,'mass':str(mass),'cubic_rhs':str(rhs),'strict_cubic_gap':str(rhs-mass)});break
 out={'center_bank_sha256':hashlib.sha256((ROOT/'round011_centers.json').read_bytes()).hexdigest(),'seed_input_sha256':bank['seed_input_sha256'],'stable_center_probes':probes,'strict_empty_examples':empty,'tail_radius':3,'exact_tail_bound':bank['exact_uniform_density_tail_bound'],'note':'diagnostics reuse fixed centers and previously recorded block states; no new search/capture center'};(ROOT/'round011_examples.json').write_text(json.dumps(out,indent=2)+'\n');print(out)
if __name__=='__main__':main()
