"""Explicit finite-block strict L1 improvement using conservative shrink losses."""
from pathlib import Path
from fractions import Fraction as F
import hashlib,json
from round005c_grid_check import build,G
ROOT=Path(__file__).resolve().parent

def main():
 fp=ROOT/'round005c_certificate_input.json';inp=json.loads(fp.read_text());grid=json.loads((ROOT/'round005c_grid_result.json').read_text());assert hashlib.sha256(fp.read_bytes()).hexdigest()==grid['input_sha256'];inner,outer,mass=build(inp)
 # Every box lies in periods <2, hence each face area <=4. Inward face
 # displacement tau costs at most 6*4*tau per box, even with overlaps.
 assert all(F(0)<F(p)<2 for p in inp['periods']);N=10**7;b=2;tau=F(1,10**12);eta=tau/4
 side=min(box[2*j+1]-box[2*j] for box in inner for j in range(3));assert 2*tau<F(side,G);assert min(map(F,inp['weights']))>tau**3
 loss=24*len(inner)*tau/mass;factor=F(2*N+1,2*N+2*b+1)**3;lower=(F(grid['ratio_interval'][0])-loss)*factor;target=F(grid['target_exact_upper']);assert lower>target
 out={'status':'explicit finite strict L1 input improves full audited checkpoint upper','input_sha256':grid['input_sha256'],'N':N,'b':b,'tau':str(tau),'bump_support_infinity_radius':str(eta),'minimum_inner_grid_side':str(F(side,G)),'box_count':len(inner),'normalized_union_shrink_loss_upper':str(loss),'finite_index_factor':str(factor),'finite_input_mass':str(mass*(2*N+2*b+1)**3),'finite_weak_ratio_lower':str(lower),'finite_lower_decimal_display':float(lower),'target_exact_upper':str(target),'finite_ratio_above_full_checkpoint_upper':lower>target,'finite_improvement_lower':str(lower-target)}
 (ROOT/'round005c_checkpoint_transfer.json').write_text(json.dumps(out,indent=2)+'\n');print({k:v for k,v in out.items() if k not in ['finite_weak_ratio_lower','target_exact_upper','finite_improvement_lower']})
if __name__=='__main__':main()
