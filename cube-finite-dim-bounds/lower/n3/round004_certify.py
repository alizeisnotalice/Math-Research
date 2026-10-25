"""Exact saved free-satellite witness certificate, no optimization on verification."""
from pathlib import Path
from fractions import Fraction as F
import argparse,hashlib,json,time
from certify import verify
ROOT=Path(__file__).resolve().parent

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--build-input',action='store_true');args=ap.parse_args()
 if args.build_input:
  r=json.loads((ROOT/'round004_results.json').read_text());c=r['final']
  inp={'periods':[format(x,'.12f') for x in c['periods']],'phases':[['0','0','0'],[format(x,'.12f') for x in c['phases'][1]]],'weights':['1',format(c['weights'][1],'.12f')],'witness_keys':json.loads((ROOT/'round004_keys.json').read_text())}
  (ROOT/'round004_certificate_input.json').write_text(json.dumps(inp,indent=2)+'\n')
 fp=ROOT/'round004_certificate_input.json';inp=json.loads(fp.read_text());start=time.monotonic();out=verify(inp)
 cert=json.loads((ROOT.parents[1]/'baseline/original_certificates.json').read_text())[0];bl,bh=map(F,cert['lower_bound_interval']);low,high=map(F,out['ratio_interval'])
 out.update({'status':'complete exact enclosure of saved satellite finite-witness union','lower_above_original_R3_upper':low>bh,'rigorous_improvement_lower':str(low-bh),'mass_per_index_pair':str(sum(map(F,inp['weights']))),'input_sha256':hashlib.sha256(fp.read_bytes()).hexdigest(),'elapsed_seconds':time.monotonic()-start})
 # Endpoints of every translated active block lie in retained [-N-b,N+b]^3.
 checks=0;b=out['b']
 for key in inp['witness_keys']:
  for orbit in range(2):
   lo=key[3*orbit:3*orbit+3];hi=key[6+3*orbit:9+3*orbit]
   if all(l<=h for l,h in zip(lo,hi)):
    for N in [0,1,7]:assert all(-N-b<=-N+lo[j]<=N+hi[j]<=N+b for j in range(3));checks+=1
 out['translation_endpoint_checks']={'N':[0,1,7],'checks':checks,'all_passed':True}
 simple=F('1.7939726044040153');assert bh<simple<low
 out['simple_rational_lower']=str(simple);out['simple_lower_above_original_upper']=simple>bh;out['simple_lower_below_certified_lower']=simple<low
 out['compatibility_field_notes']={'rigorous_gap_lower':'legacy field original R3 lower minus new upper; negative is expected for an improvement','strictly_below_original_R3':'legacy comparison flag; false is expected for an improvement','lower_above_original_R3_upper':'authoritative strict improvement comparison','rigorous_improvement_lower':'authoritative new lower minus original upper'}
 (ROOT/'round004_certificate_result.json').write_text(json.dumps(out,indent=2)+'\n')
 print({k:v for k,v in out.items() if k!='radius_enclosures'})
if __name__=='__main__':main()
