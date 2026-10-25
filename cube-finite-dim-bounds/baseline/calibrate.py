"""Exact geometric calibrations; numerical values are diagnostic only."""
from fractions import Fraction as F
from itertools import product
import math,json
from pathlib import Path
out=[]
for n in (3,4):
 vertices=list(product((-1,1),repeat=n))
 eta=F(1,10); r=1+2*eta; eps=F(1,100); delta=F(1,10000)
 for v in vertices:
  assert max(abs(t) for t in v)==1
  assert [w for w in vertices if max(abs(a-b) for a,b in zip(v,w))<=1]==[v]
 # entire origin and own vertex blobs included; every other blob excluded
 assert 1+2*eta<=r and 2*eta<=r and 2-2*eta>r
 smooth=F(2**n)/(1+eps+delta*(2*r)**n)
 u=2**(1/n);a=(1+u)/2
 al=1.5*a**(n-1)-.5*((3*a-2*1.5**(1/n))**(n-1)-((3-u)/2)**(n-1))
 out.append(dict(dimension=n,atomic_ratio=str(F(2**n)/(1+eps)),smooth_ratio=str(smooth),aldaz_float=al))
# Frobenius-sign average v^T B v = tr B checked over rational non-diagonal SPD example.
for n in (3,4):
 B=[[F(i==j)+F(1,10) for j in range(n)] for i in range(n)]
 vals=[sum(v[i]*B[i][j]*v[j] for i in range(n) for j in range(n)) for v in product((-1,1),repeat=n)]
 assert sum(vals)/2**n==sum(B[i][i] for i in range(n))
 assert max(vals)>=sum(B[i][i] for i in range(n))
result={'calibrations':out,'radial_cost_float':[27*math.pi*math.sqrt(3)/16,2*math.pi**2], 'exact_assertions':'PASS','note':'finite examples supplement, not replace, analytic proof'}
print(json.dumps(result,indent=2))
if __name__=='__main__':
 import sys
 if '--save' in sys.argv: (Path(__file__).resolve().parents[1]/'results/calibration_run.json').write_text(json.dumps(result,indent=2)+'\n')
