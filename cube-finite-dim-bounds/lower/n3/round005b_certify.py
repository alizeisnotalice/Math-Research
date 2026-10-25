"""Read-only Fraction certificate for round005 versus full T3 exact upper."""
from pathlib import Path
from fractions import Fraction as F
import hashlib,json,time
from certify import verify
ROOT=Path(__file__).resolve().parent

def main():
 start=time.monotonic();fp=ROOT/'round005b_certificate_input.json';inp=json.loads(fp.read_text());out=verify(inp);old=json.loads((ROOT/'round005_certificate_input.json').read_text());oldkeys={tuple(k) for k in old['witness_keys']};keys={tuple(k) for k in inp['witness_keys']};assert oldkeys<=keys
 target=F(json.loads((ROOT/'round005b_results.json').read_text())['target_checkpoint_upper']);lo,hi=map(F,out['ratio_interval']);out.pop('strictly_below_original_R3',None);out.pop('rigorous_gap_lower',None)
 out.update({'status':'complete exact finite-union enclosure; compare full checkpoint upper','target_exact_upper':str(target),'lower_above_checkpoint_exact_upper':lo>target,'rigorous_improvement_lower':str(lo-target),'original_170_preserved':oldkeys<=keys,'original_key_count':len(oldkeys),'new_key_count':len(keys-oldkeys),'input_sha256':hashlib.sha256(fp.read_bytes()).hexdigest(),'seed_input_sha256':hashlib.sha256((ROOT/'round005_certificate_input.json').read_bytes()).hexdigest(),'elapsed_seconds':time.monotonic()-start})
 checks=0;b=out['b']
 for key in inp['witness_keys']:
  for orbit in range(2):
   low=key[3*orbit:3*orbit+3];high=key[6+3*orbit:9+3*orbit]
   if all(l<=h for l,h in zip(low,high)):
    for N in [0,1,7]:assert all(-N-b<=-N+low[j]<=N+high[j]<=N+b for j in range(3));checks+=1
 out['translation_endpoint_checks']={'N':[0,1,7],'checks':checks,'all_passed':True};assert lo>target
 simple=F('1.7939758676554228');assert target<simple<lo
 out['simple_rational_lower']=str(simple);out['simple_lower_below_certified_lower']=simple<lo
 previous=F(json.loads((ROOT/'round005_results.json').read_text())['target_exact_upper']);out['gain_above_round004_T3_upper']=str(lo-previous);assert lo>previous
 (ROOT/'round005b_certificate_result.json').write_text(json.dumps(out,indent=2)+'\n');print({k:v for k,v in out.items() if k!='radius_enclosures'})
if __name__=='__main__':main()
