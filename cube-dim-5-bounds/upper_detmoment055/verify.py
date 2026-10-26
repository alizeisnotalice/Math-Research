#!/usr/bin/env python3
"""Independent Fraction-only reconstruction of saved quadratic/root/gate bounds."""
from fractions import Fraction as F
from pathlib import Path
import hashlib,json,os,time
HERE=Path(__file__).resolve().parent
MAIN=Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube-dim-5/cube-dim-5-bounds')
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 start=time.monotonic();r=json.loads((HERE/'results.json').read_text());reg=json.loads((HERE/'registration.json').read_text())
 assert r['source_sha256']==sha(HERE/'certify.py') and r['registration_sha256']==sha(HERE/'registration.json')
 for path,record in reg['sources'].items():assert sha(Path(path))==record['sha256']
 for directory,names in [('astra_upper015',['graph_moments_beta_half.hpp']),('astra_upper016',['hermite_moments.json'])]:
  manifest=json.loads((MAIN/directory/'manifest.json').read_text())
  for name in names:
   row=next(x for x in manifest['white_list'] if x['path']==name)
   assert row['sha256']==sha(MAIN/directory/name)
 raw=json.loads((MAIN/'astra_upper016/hermite_moments.json').read_text())
 m=[F(1)]+[F(raw['moments'][str(k)][0]) for k in range(1,4)];assert m==list(map(F,r['unweighted_moments_0_to_3']))
 variance=m[2]-m[1]**2;A=(m[3]-m[1]*m[2])/variance;B=(m[1]*m[3]-m[2]**2)/variance
 assert [F(1),-A,B]==list(map(F,r['p2_coefficients']))
 assert m[2]==A*m[1]-B and m[3]==A*m[2]-B*m[1]
 disc=A*A-4*B;sl,sh=map(F,r['exact_sqrt_discriminant_interval']);assert sl*sl<=disc<sh*sh
 nodes=[tuple(map(F,x)) for x in r['nodes_rational_intervals']]
 assert nodes==[((A-sh)/2,(A-sl)/2),((A+sl)/2,(A+sh)/2)]
 (l1,h1),(l2,h2)=nodes;p=lambda x:x*x-A*x+B
 assert 0<l1<h1<m[1]<l2<h2<1 and p(l1)>0>p(h1) and p(l2)<0<p(h2)
 expected_weights=[((l2-m[1])/(h2-l1),(h2-m[1])/(l2-h1)),((m[1]-h1)/(h2-l1),(m[1]-l1)/(l2-h1))]
 weights=[tuple(map(F,x)) for x in r['weights_rational_intervals']];assert weights==expected_weights
 fs=[tuple(map(F,x)) for x in r['node_function_tenth_root_intervals']]
 for (xlo,xhi),(flo,fhi),(wlo,whi) in zip(nodes,fs,weights):
  assert 0<wlo<whi<1 and flo**10<=1/xhi and fhi**10>1/xlo
 jlo=sum(w[0]*f[0] for w,f in zip(weights,fs));jhi=sum(w[1]*f[1] for w,f in zip(weights,fs))
 assert [jlo,jhi]==list(map(F,r['Gauss_value_rational_interval']))
 budget=json.loads((MAIN/'round010/orbit_budget.json').read_text());pl,ph=map(F,budget['pi_machin80_interval'])
 gate_lo,gate_hi=map(F,r['J_gate_rational_interval']);g=F(1,10000);coef=F(60**4*32**9,3125)
 assert gate_lo**5<=g**5*coef*pl**7 and gate_hi**5>g**5*coef*ph**7
 assert jhi<gate_lo
 assert F(r['J_lower_minus_gate_upper_exact'])==jlo-gate_hi<0
 assert r['new_graph_states']==0 and r['graph_backend_calls']==0 and not r['cache_expanded']
 out={'status':'PASS_independent_fraction_checker_local_no_exclusion','pid':os.getpid(),'source_sha256':sha(Path(__file__)),'results_sha256':sha(HERE/'results.json'),
  'locked_beta_half_source_hashes_match':True,'unweighted_existing_moment_selection_matches':True,'quadratic_orthogonality_exact':True,'positive_nodes_and_weights_verified':True,
  'integer_tenth_root_comparisons_verified':True,'Gauss_upper_strictly_below_gate_lower':True,'conclusion':'No separation from this two-node lower; do not infer actual I/C5 or all-orbit exclusion',
  'Gauss_upper_minus_gate_lower_exact':str(jhi-gate_lo),'wall_seconds':time.monotonic()-start,'new_graph_states':0}
 (HERE/'verification.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({k:out[k] for k in ('status','pid','wall_seconds')}))
if __name__=='__main__':main()
