#!/usr/bin/env python3
"""Pure postprocessing of completed time nodes/gradients with F's sparse shift."""
from fractions import Fraction as Q
from pathlib import Path
import hashlib,json
from flint import arb
import run_moments as m

HERE=Path(__file__).resolve().parent
PROOF=HERE.parents[2]/'audit/round007/sparse_shift_bound.md'
PROOF_SHA='c386ed307d9df8379c4dcca9983bf1571405b112b093c4bb9c9c0cd1dc04fede'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def phi(u):return(12*(u/12).asinh()).exp()
def sparse_u2(r,l,vals,masses,details=False):
 parts=[m.A(m.TAUS[k])**-2*phi(m.A(r)/m.A(m.TAUS[k]))*masses[k][1]for k in range(l,7)]
 tail=m.A(m.TAUS[7])**-2*phi(m.A(r)/m.A(m.TAUS[7]))*vals[7].upper()
 value=sum(parts,arb(0))+tail
 return(value,parts,tail)if details else value
def sparse_band_lower(r,l,vals,masses):
 return sum((masses[k][0]/phi(m.A(r)/m.A(m.TAUS[k]))for k in range(l,7)),arb(0))+vals[7].lower()/phi(m.A(r)/m.A(m.TAUS[7]))
def main():
 assert sha(PROOF)==PROOF_SHA
 original=json.loads((HERE/'results.json').read_text());assert original['status']=='COMPLETE_LOCAL_ONLY'
 out=json.loads(json.dumps(original));m.U2=sparse_u2;m.band_lower=sparse_band_lower
 m.finalize(out,out['center_numerators'],out['denominator'])
 vals,masses=m.bands(out['time_nodes'])
 for box,gr in zip(out['taylor_boxes'],out['gradients']):
  lower=m.taylor_lower(Q(box['radius']),gr['time_index'],out['center_numerators'],out['denominator'],
   vals,masses,[arb(v)for v in gr['gradient_x2_to_x6']])[0]
  lo,hi=0,1000000000
  while hi-lo>1:
   mid=(lo+hi)//2
   if lower>arb(mid)/10**12:lo=mid
   else:hi=mid
  box['G_lower_exact_rational']=str(Q(lo,10**12))
  box['G_lower_endpoint_arb384']=lower.lower().str(100,more=True)
  assert Q(lo,10**12)>Q(1,60000)
 out['postprocessing_only']=True;out['new_kernel_nodes_in_postprocessing']=0
 out['space_shift_factor']='Phi(u)=exp(12*asinh(u/12))'
 out['sparse_shift_proof']=str(PROOF);out['sparse_shift_proof_sha256']=sha(PROOF)
 out['original_results_sha256']=sha(HERE/'results.json');out['postprocessor_sha256']=sha(Path(__file__))
 out['original_exponential_taylor_radii']=[{'epsilon':r['epsilon'],'radius':r['radius']}for r in original['taylor_boxes']]
 out['best_sparse_over_best_exponential']=str(Q(out['best_taylor_box']['radius'])/Q(original['best_taylor_box']['radius']))
 (HERE/'sparse_shift_results.json').write_text(json.dumps(out,indent=2)+'\n')
 print('original multiplicative r',out['best_original_radius'])
 print('old exponential Taylor r',original['best_taylor_box']['radius'])
 print('new sparse Taylor r',out['best_taylor_box']['radius'])
 print('epsilon',out['best_taylor_box']['epsilon'],'gain over original',out['best_taylor_box']['radius_over_best_original'])
if __name__=='__main__':main()
