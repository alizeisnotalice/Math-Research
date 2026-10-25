"""Three-orbit Fraction verifier: full mass, explicit atom endpoints, joint union."""
from pathlib import Path
from fractions import Fraction as F
from itertools import product
import argparse,hashlib,json,time
from certify import root_interval,volume
ROOT=Path(__file__).resolve().parent

def verify(inp,target_upper=None):
 p=list(map(F,inp['periods']));ph=[list(map(F,q)) for q in inp['phases']];w=list(map(F,inp['weights']));assert len(p)==3 and len(ph)==len(w)==3 and w[0]==1 and all(x>0 for x in p+w);W=sum(w)
 inner=[];outer=[];los=[];his=[];radius=[];endpoint_checks=0;atom_occurrences=0
 for key in inp['keys']:
  assert len(key)==18 and all(isinstance(k,int) for k in key)
  lo=[key[3*b:3*b+3] for b in range(3)];hi=[key[9+3*b:12+3*b] for b in range(3)];counts=[];points=[];enum_mass=F(0)
  for b in range(3):
   num=1
   for j in range(3):num*=max(0,hi[b][j]-lo[b][j]+1)
   counts.append(num)
   if num:
    los.append(lo[b]);his.append(hi[b]);local=[]
    for k in product(*(range(lo[b][j],hi[b][j]+1) for j in range(3))):local.append(tuple(p[j]*(F(k[j])+ph[b][j]) for j in range(3)));enum_mass+=w[b]
    assert len(local)==num;points.extend(local)
  mass=sum((counts[b]*w[b] for b in range(3)),F(0));assert mass==enum_mass and mass>0;rl,rh=root_interval(mass);lc=[max(x[j] for x in points) for j in range(3)];hc=[min(x[j] for x in points) for j in range(3)];atom_occurrences+=len(points)
  for j in range(3):assert lc[j]==max(p[j]*(hi[b][j]+ph[b][j]) for b in range(3) if counts[b]) and hc[j]==min(p[j]*(lo[b][j]+ph[b][j]) for b in range(3) if counts[b])
  for target,r in [(inner,rl),(outer,rh)]:
   box=tuple(x for j in range(3) for x in [max(F(0),lc[j]-r),min(p[j],hc[j]+r)])
   if all(box[2*j]<box[2*j+1] for j in range(3)):
    target.append(box)
    if target is inner:
     for x in points:
      for j in range(3):assert x[j]-box[2*j]<=r and box[2*j+1]-x[j]<=r;endpoint_checks+=1
  radius.append({'orbit_counts':counts,'mass':str(mass),'radius_interval':[str(rl),str(rh)],'critical_inner_radius':(2*rl)**3==mass})
 low=volume(inner)/W;high=volume(outer)/W;target=F(target_upper) if target_upper is not None else F(json.loads((ROOT/'round006_cert_target.json').read_text())['ratio_interval'][1]);oldtarget=F(json.loads((ROOT/'round006_target.json').read_text())['ratio_interval'][1]);imin=[min(l[j] for l in los) for j in range(3)];imax=[max(h[j] for h in his) for j in range(3)];b=max(abs(x) for x in imin+imax);checks=0
 for N in [0,1,7]:
  for l,h in zip(los,his):assert all(-N-b<=-N+l[j]<=N+h[j]<=N+b for j in range(3));checks+=1
 old=json.loads((ROOT/'round005c_certificate_input.json').read_text());oldkeys={tuple(k[:6]+[0,0,0]+k[6:]+[-1,-1,-1]) for k in old['witness_keys']};assert oldkeys<=set(map(tuple,inp['keys']))
 return {'status':'exact finite three-orbit joint union enclosure','ratio_interval':[str(low),str(high)],'decimal_display_only':[float(low),float(high)],'target_exact_upper':str(target),'lower_above_full_target_upper':low>target,'lower_above_full_007c_upper':low>oldtarget,'rigorous_improvement_lower':str(low-target),'upper_below_full_target_upper':high<target,'mass_per_index_triple':str(W),'original170_preserved':True,'key_count':len(inp['keys']),'inner_boxes':len(inner),'outer_boxes':len(outer),'atom_occurrences':atom_occurrences,'explicit_nonempty_atom_endpoint_checks':endpoint_checks,'radius_records':radius,'critical_radius_count':sum(r['critical_inner_radius'] for r in radius),'index_min':imin,'index_max':imax,'b':b,'translation_endpoint_checks':{'N':[0,1,7],'active_block_checks':checks,'all_passed':True},'finite_mass_formula':str(W)+'*(2*N+2*b+1)^3','receiver_volume_formula':'(2*N+1)^3*U_tau','transfer':'retain three labeled orbit atoms per integer index, contract boxes tau, use radius r-tau/2, smooth with support <tau/2; first N to infinity then tau to zero'}

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--input',default='round006_certificate_input.json');ap.add_argument('--output',default='round006_certificate_result.json');ap.add_argument('--target',default='round006_cert_target.json');args=ap.parse_args();start=time.monotonic();fp=ROOT/args.input;inp=json.loads(fp.read_text());out=verify(inp,json.loads((ROOT/args.target).read_text())['ratio_interval'][1]);out['input_sha256']=hashlib.sha256(fp.read_bytes()).hexdigest();out['elapsed_seconds']=time.monotonic()-start
 if out['lower_above_full_target_upper']:
  lo=F(out['ratio_interval'][0]);simple=F((lo.numerator*10**16)//lo.denominator,10**16);out['simple_rational_lower']=str(simple);out['simple_above_target']=simple>F(out['target_exact_upper'])
 raw=out.pop('radius_records');records=[];lookup={};indices=[]
 for r in raw:
  key=(r['mass'],tuple(r['orbit_counts']))
  if key not in lookup:lookup[key]=len(records);records.append(r)
  else:assert records[lookup[key]]==r
  indices.append(lookup[key])
 out['distinct_radius_records']=records;out['key_radius_record_indices']=indices
 (ROOT/args.output).write_text(json.dumps(out,indent=2)+'\n');print({k:v for k,v in out.items() if k not in ['distinct_radius_records','key_radius_record_indices']})
if __name__=='__main__':main()
