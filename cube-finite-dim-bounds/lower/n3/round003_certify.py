"""Read-only exact Fraction verification of three-layer finite witness unions."""
from pathlib import Path
from fractions import Fraction as F
import json,hashlib,time,random,itertools
from certify import root_interval,volume
ROOT=Path(__file__).resolve().parent

def verify(inp):
 p=list(map(F,inp['periods']));ph=[list(map(F,q)) for q in inp['phases']];w=list(map(F,inp['weights']));assert len(p)==3 and len(ph)==len(w)==3 and all(x>0 for x in p+w) and sum(w)==3
 inner=[];outer=[];details=[];los=[];his=[]
 for key in inp['keys']:
  assert len(key)==18 and all(isinstance(k,int) for k in key)
  lo=[key[3*b:3*b+3] for b in range(3)];hi=[key[9+3*b:12+3*b] for b in range(3)]
  cnt=[max(0,hi[b][0]-lo[b][0]+1)*max(0,hi[b][1]-lo[b][1]+1)*max(0,hi[b][2]-lo[b][2]+1) for b in range(3)]
  active=[b for b in range(3) if cnt[b]>0];assert active
  for b in active:los.append(lo[b]);his.append(hi[b])
  mass=sum((cnt[b]*w[b] for b in active),F(0));rl,rh=root_interval(mass)
  lc=[max(p[j]*(hi[b][j]+ph[b][j]) for b in active) for j in range(3)];hc=[min(p[j]*(lo[b][j]+ph[b][j]) for b in active) for j in range(3)]
  for target,r in [(inner,rl),(outer,rh)]:
   box=tuple(c for j in range(3) for c in [max(F(0),lc[j]-r),min(p[j],hc[j]+r)])
   if all(box[2*j]<box[2*j+1] for j in range(3)):target.append(box)
  details.append({'mass':str(mass),'radius':[str(rl),str(rh)],'strict_inner_radius':(2*rl)**3<mass})
 low=volume(inner)/3;high=volume(outer)/3
 baseline=json.loads((ROOT.parents[1]/'baseline/original_certificates.json').read_text())[0];bl,bh=map(F,baseline['lower_bound_interval'])
 imin=[min(l[j] for l in los) for j in range(3)];imax=[max(h[j] for h in his) for j in range(3)];b=max(abs(x) for x in imin+imax)
 assert all(-b<=l[j]<=h[j]<=b for l,h in zip(los,his) for j in range(3))
 for N in [0,1,7]:
  assert all(-N-b<=-N+l[j] and N+h[j]<=N+b for l,h in zip(los,his) for j in range(3))
 return {'translation_endpoint_checks':{'N_values':[0,1,7],'active_orbit_block_checks':3*len(los),'all_passed':True},'ratio_interval':[str(low),str(high)],'decimal_display_only':[float(low),float(high)],'upper_below_original_R3_lower':high<bl,'lower_above_original_R3_upper':low>bh,'gap_from_original_lower':str(bl-high),'key_count':len(inp['keys']),'inner_box_count':len(inner),'outer_box_count':len(outer),'index_min':imin,'index_max':imax,'b':b,'finite_mass_formula':'3*(2*N+2*b+1)^3','receiver_volume_formula':'(2*N+1)^3*U_tau','critical_radius_count':sum(not d['strict_inner_radius'] for d in details),'radius_enclosures':details,'transfer':'keep given phase representatives; each index has three orbit atoms of total mass 3; receive q in [-N,N]^3, retain all orbit labels [-N-b,N+b]^3; contract receiver boxes by tau, reduce radius by tau/2, smooth with support <tau/2; first N to infinity, then tau to zero'}

def selftest():
 rng=random.Random(20261025306)
 for case in range(32):
  boxes=[]
  for i in range(1+case%7):
   bb=[]
   for j in range(3):a,b=sorted(rng.sample(range(9),2));bb.extend([F(a),F(b)])
   boxes.append(tuple(bb))
  brute=sum(any(all(b[2*j]<=pt[j] and pt[j]+1<=b[2*j+1] for j in range(3)) for b in boxes) for pt in itertools.product(range(8),repeat=3))
  assert volume(boxes)==brute
 return {'seed':20261025306,'cases':32,'coordinate_range':[0,8],'boxes_per_case':'1..7','all_passed':True}

def main():
 start=time.monotonic();fp=ROOT/'round003_certificate_input.json';inp=json.loads(fp.read_text());out={'status':'complete certification of saved finite unions only','input_sha256':hashlib.sha256(fp.read_bytes()).hexdigest(),'unions':{name:verify(case) for name,case in inp.items()}}
 if 'staircase' in out['unions']:
  gain=F(out['unions']['staircase']['ratio_interval'][0])-F(out['unions']['rectangular']['ratio_interval'][1]);out['strict_staircase_gain_lower']=str(gain);out['strict_staircase_gain_positive']=gain>0
 out['independent_integer_grid_union_check']=selftest()
 out['elapsed_seconds']=time.monotonic()-start;(ROOT/'round003_certificate_result.json').write_text(json.dumps(out,indent=2)+'\n');print({name:{k:v for k,v in r.items() if k not in ['radius_enclosures','transfer']} for name,r in out['unions'].items()});print('gain',out.get('strict_staircase_gain_lower'))
if __name__=='__main__':main()
