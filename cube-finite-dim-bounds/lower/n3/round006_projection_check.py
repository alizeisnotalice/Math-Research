"""Counterfactual finite-witness attribution at fixed original two-orbit input."""
from pathlib import Path
from fractions import Fraction as F
import hashlib,json
from certify import verify
ROOT=Path(__file__).resolve().parent

def main():
 new=json.loads((ROOT/'round006_structural_input.json').read_text());old=json.loads((ROOT/'round005c_certificate_input.json').read_text());keys=set(map(tuple,old['witness_keys']))
 for k in new['keys']:
  z=k[:6]+k[9:15]
  if any(all(z[6+3*b+j]>=z[3*b+j] for j in range(3)) for b in range(2)):keys.add(tuple(z))
 inp={'periods':old['periods'],'phases':old['phases'],'weights':old['weights'],'witness_keys':[list(k) for k in sorted(keys)]};fp=ROOT/'round006_projection_input.json';fp.write_text(json.dumps(inp,indent=2)+'\n');out=verify(inp);three=F(json.loads((ROOT/'round006_structural_certificate.json').read_text())['ratio_interval'][0]);projected_upper=F(out['ratio_interval'][1]);assert three>projected_upper
 result={k:v for k,v in out.items() if k not in ['radius_enclosures','rigorous_gap_lower','strictly_below_original_R3']};result.update({'status':'finite projected two-orbit counterfactual at unchanged original input','input_sha256':hashlib.sha256(fp.read_bytes()).hexdigest(),'projected_key_count':len(keys),'structural_three_orbit_lower':str(three),'three_orbit_lower_above_projected_two_orbit_upper':three>projected_upper,'rigorous_finite_attribution_gain_lower':str(three-projected_upper),'interpretation':'strict comparison only between these finite witness unions; not full superlevel sets'});(ROOT/'round006_projection_result.json').write_text(json.dumps(result,indent=2)+'\n');print(result)
if __name__=='__main__':main()
