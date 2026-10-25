#!/usr/bin/env python3
from pathlib import Path
from fractions import Fraction as Q
import hashlib,json,sys
HERE=Path(__file__).resolve().parent;ROOT=HERE.parent
sys.path.insert(0,str(ROOT/'lower'))
from certify_rectangular import field
F,u,v,bounds,_,_=field()
def decode(d):
 a=[Q(0)]*25
 for key,value in d.items():i,j=map(int,key.split(','));a[5*i+j]=Q(value)
 return F(a=a)
cell=json.loads((HERE/'event_cell.json').read_text());res=json.loads((HERE/'cell_derivative_results.json').read_text());jet=json.loads((HERE/'jet_check_results.json').read_text());C=decode(json.loads((ROOT/'lower_round003/a2_equals_a3_results.json').read_text())['coefficients'])
L=decode(cell['lower_coefficients']);U=decode(cell['upper_coefficients'])
assert L==(-v+u*v**2)/2 and U==(1+u+u**2-u**3*v**3)/2
assert L<(1+u)/8 and (1+u)/8<U and F(0)<L and U<F(Q(1,2))
ds=[]
for i,node in enumerate(cell['nodes']):
 t=decode(node['position_coefficients']);assert t==L+F(Q(node['theta']))*(U-L)
 half=[]
 for sign in ('plus','minus'):
  data=json.loads((HERE/f'kernel_{i}_{sign}.json').read_text());val=decode(data['kernel_coefficients'])
  assert sum((decode(x['contribution']) for x in data['nonzero_faces']),F(0))==val
  assert sum((decode(x) for x in data['axis_contributions']),F(0))==val
  half.append(val)
 d=16*sum(half,F(0))-C;ds.append(d)
 assert d==decode(res['node_derivatives'][i]['coefficients'])
assert ds[0]+ds[2]==2*ds[1]
assert decode(jet['derivative_coefficients'])==ds[1]
for key,d in [('left_limit',2*ds[0]-ds[1]),('right_limit',2*ds[2]-ds[1])]:
 assert decode(res[key]['coefficients'])==d
 lo,hi=bounds(d.a);assert hi<0
 assert list(map(str,(lo,hi)))==res[key]['interval']
out={'status':'all_saved_coefficient_identities_and_strict_signs_verified','source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'checks':['exact interval endpoints and3 registered positions','oriented face and coordinate sums','reflection factor16 and denominator derivative-C','exact affine midpoint identity','independent mixed-box midpoint jet equality','both limiting endpoint isolating upper bounds strictly negative']}
(HERE/'verification.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
