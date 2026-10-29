"""Independent small exact-Fraction checks of normalization and count-polynomial expansion."""
from fractions import Fraction as F
from itertools import product
from pathlib import Path
import json, math, hashlib
import probe as p
H=Path(__file__).resolve().parent
checks=[]
# Source constants and all cells, as rationals, not decimal binary approximation.
eps=F(2,25);delta=F(1,25)
for phases in p.CFG['phases']:
 for j,phase in enumerate(phases):
  center=(j+F(str(phase)))/2
  assert F(j,2)<center-delta/2<center+delta/2<F(j+1,2)
  assert eps/2+(1-eps)/(2*delta)*delta==F(1,2)
checks.append('All 8 internal boxes and exact mass 1/2 per cell')
# Actual source overlap formula independently against rational period intersections.
import numpy as np
xs=[F(i,16)for i in range(1,16)]
ff,_,_=p.features(np.array([[float(x)]for x in xs]),p.CFG)
max_kernel_error=0.
for b,phases in enumerate(p.CFG['phases']):
 for qi,x in enumerate(xs):
  for li,ll in enumerate(p.CFG['lengths']):
   L=F(str(ll))
   def ov(left,right):
    return sum(max(F(0),min(x+L/2,right+m)-max(x-L/2,left+m))for m in range(-4,5))
   for j,phase in enumerate(phases):
    center=(j+F(str(phase)))/2
    exact=(eps*ov(F(j,2),F(j+1,2))+(1-eps)/(2*delta)*ov(center-delta/2,center+delta/2))/L
    err=abs(float(exact)-ff[b,qi,li,j]);max_kernel_error=max(max_kernel_error,err)
    assert err<2e-14
checks.append(f'2160 actual n1 source kernels versus independent rational overlap: max double error {max_kernel_error}')

for n in p.CFG['dimensions']:
 f=json.loads((H/f'FROZEN_n{n}.json').read_text());s=f['selector'];nums=[sum(math.comb(n,k)for k in range(n+1)if s[k]==b)for b in range(4)]
 assert nums==f['prior_exact_numerators'];assert sum(nums)==2**n
 g=json.loads((H/f'BALANCED_FROZEN_n{n}.json').read_text());weights=[[F(str(v))for v in arr]for arr in g['selector_weights']]
 assert all(sum(weights[b][k]for b in range(2))==1 for k in range(n+1))
 assert all(sum(F(math.comb(n,k),2**n)*weights[b][k]for k in range(n+1))==F(1,2)for b in range(2))
 checks.append(f'n={n}: exact primary prior and balanced prior mass')
# Arbitrary exact positive z, direct full labels against polynomial recurrence and selector.
examples=0
for n in range(1,8):
 for rep in range(5):
  zs=[[(F(1+(i+rep+2*b)%7,19),F(2+(2*i+rep+b)%5,17))for i in range(n)]for b in range(2)]
  co=[]
  for b in range(2):
   cur=[F(1)]
   for z0,z1 in zs[b]:
    nxt=[F(0)]*(len(cur)+1)
    for k,v in enumerate(cur):nxt[k]+=v*z0;nxt[k+1]+=v*z1
    cur=nxt
   co.append(cur)
  weights=[[F((k+rep)%3,2)for k in range(n+1)]];weights.append([1-w for w in weights[0]])
  dp=sum(co[b][k]*weights[b][k]for b in range(2)for k in range(n+1));direct=F(0)
  for labels in product(range(2),repeat=n):
   k=sum(labels)
   for b in range(2):direct+=weights[b][k]*math.prod(zs[b][i][labels[i]]for i in range(n))
  assert direct==dp;examples+=1
checks.append(f'{examples} exact positive coefficient-DP versus full-label expansions, n<=7')
for base in ['', 'BALANCED_']:
 for n in p.CFG['dimensions']:
  r=json.loads((H/f'{base}RESULT_n{n}.json').read_text());assert all(t['near_guard']==0 and t['k_guard']==t['k_raw']for t in r['rows'])
checks.append('All six frozen runs: no selected-lambda guard crossings; raw=guard count')
record=dict(passed=True,checks=checks,not_checked=['Continuous query cell classification','Outward interval integration and DP','PRNG ideal iid assumption','Global or asymptotic growth optimization'],sha256={f.name:hashlib.sha256(f.read_bytes()).hexdigest()for f in H.glob('*.py')})
p.dump(H/'VERIFICATION.json',record);print(json.dumps(record,indent=2))
