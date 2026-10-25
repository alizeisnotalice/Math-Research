"""Explicit finite strict L1 lower ratio for three labeled orbits."""
from pathlib import Path
from fractions import Fraction as F
import argparse,hashlib,json
from round006_grid_check import build,G
ROOT=Path(__file__).resolve().parent

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--input',default='round006_certificate_input.json');ap.add_argument('--grid',default='round006_grid_result.json');ap.add_argument('--output',default='round006_finite_result.json');args=ap.parse_args();fp=ROOT/args.input;inp=json.loads(fp.read_text());grid=json.loads((ROOT/args.grid).read_text());assert hashlib.sha256(fp.read_bytes()).hexdigest()==grid['input_sha256'];inner,outer,W=build(inp);assert W==sum(map(F,inp['weights'])) and len(inp['weights'])==3
 lows=[];highs=[]
 for key in inp['keys']:
  for orbit in range(3):
   lo=key[3*orbit:3*orbit+3];hi=key[9+3*orbit:12+3*orbit]
   if all(l<=h for l,h in zip(lo,hi)):lows.append(lo);highs.append(hi)
 imin=[min(l[j] for l in lows) for j in range(3)];imax=[max(h[j] for h in highs) for j in range(3)];b=max(abs(x) for x in imin+imax);assert all(-b<=l[j]<=h[j]<=b for l,h in zip(lows,highs) for j in range(3));N=10**9;tau=F(1,10**13);eta=tau/4
 assert all(0<F(p)<2 for p in inp['periods']);side=min(box[2*j+1]-box[2*j] for box in inner for j in range(3));assert 2*tau<F(side,G) and min(map(F,inp['weights']))>tau**3
 loss=24*len(inner)*tau/W;factor=F(2*N+1,2*N+2*b+1)**3;lower=(F(grid['ratio_interval'][0])-loss)*factor;target=F(grid['target_exact_upper']);assert lower>target
 out={'status':'explicit finite strict L1 triple-orbit improvement PASS','input_sha256':grid['input_sha256'],'N':N,'b':b,'index_min':imin,'index_max':imax,'tau':str(tau),'bump_support_infinity_radius':str(eta),'minimum_inner_side':str(F(side,G)),'nonempty_inner_boxes':len(inner),'mass_per_index_triple':str(W),'finite_input_mass':str(W*(2*N+2*b+1)**3),'normalized_union_loss_upper':str(loss),'finite_index_factor':str(factor),'finite_weak_ratio_lower':str(lower),'decimal_display_only':float(lower),'target_exact_upper':str(target),'finite_ratio_above_full_target_upper':lower>target,'finite_gain_lower':str(lower-target),'label_convention':'for every integer k in [-N-b,N+b]^3 retain all3 orbit representatives with their exact masses; receive q in [-N,N]^3'};(ROOT/args.output).write_text(json.dumps(out,indent=2)+'\n');print({k:v for k,v in out.items() if k not in ['finite_weak_ratio_lower','target_exact_upper','finite_gain_lower']})
if __name__=='__main__':main()
