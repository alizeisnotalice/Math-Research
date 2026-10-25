"""New candidate verification, complete012b target and all619 seed keys retained."""
from pathlib import Path
from fractions import Fraction as F
import hashlib,json,time
from round006_certify import verify
ROOT=Path(__file__).resolve().parent

def main():
 start=time.monotonic();fp=ROOT/'round007_certificate_input.json';inp=json.loads(fp.read_text());target=F(json.loads((ROOT/'round007_target.json').read_text())['ratio_interval'][1]);out=verify(inp,target);seed=json.loads((ROOT/'round006_certificate_input.json').read_text());old=set(map(tuple,seed['keys']));current=set(map(tuple,inp['keys']));assert old<=current and len(old)==619
 out.update({'status':'exact three-orbit finite-union candidate versus full012b','old619_preserved':True,'old_key_count':619,'new_key_count':len(current-old),'input_sha256':hashlib.sha256(fp.read_bytes()).hexdigest(),'seed_input_sha256':hashlib.sha256((ROOT/'round006_certificate_input.json').read_bytes()).hexdigest(),'elapsed_seconds':time.monotonic()-start});lo,hi=map(F,out['ratio_interval']);out['lower_above_full_012b_upper']=lo>target;out['rigorous_gain_above_full_012b_upper']=str(lo-target)
 if lo>target:
  simple=F((lo.numerator*10**16)//lo.denominator,10**16);out['simple_rational_lower']=str(simple);out['simple_above_full_target']=simple>target
 raw=out.pop('radius_records');table=[];lookup={};indices=[]
 for r in raw:
  k=(r['mass'],tuple(r['orbit_counts']))
  if k not in lookup:lookup[k]=len(table);table.append(r)
  else:assert table[lookup[k]]==r
  indices.append(lookup[k])
 out['distinct_radius_records']=table;out['key_radius_record_indices']=indices
 (ROOT/'round007_certificate_result.json').write_text(json.dumps(out,indent=2)+'\n');print({k:v for k,v in out.items() if k not in ['distinct_radius_records','key_radius_record_indices']})
if __name__=='__main__':main()
