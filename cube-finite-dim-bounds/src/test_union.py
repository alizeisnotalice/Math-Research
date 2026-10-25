"""Deterministic dyadic geometry checks against independent coordinate-cell enumeration.
This checks the numerical search engine, not an entire candidate certificate.
"""
from itertools import product
from fractions import Fraction as F
from pathlib import Path
import random,json,time
from lattice_search import volume

def grid_volume(boxes):
 boxes=[b for b in boxes if all(b[2*j]<b[2*j+1] for j in range(len(b)//2))]
 if not boxes:return F(0)
 n=len(boxes[0])//2
 cuts=[sorted({b[2*j+k] for b in boxes for k in [0,1]}) for j in range(n)]
 total=F(0)
 for ids in product(*[range(len(c)-1) for c in cuts]):
  if any(all(b[2*j]<=cuts[j][i] and cuts[j][i+1]<=b[2*j+1] for j,i in enumerate(ids)) for b in boxes):
   v=F(1)
   for j,i in enumerate(ids):v*=cuts[j][i+1]-cuts[j][i]
   total+=v
 return total

def main():
 start=time.monotonic();rng=random.Random(20261025);count=0
 for n in [2,3,4]:
  for _ in range(20):
   boxes=[]
   for k in range(5):
    b=[]
    for j in range(n):
     x,y=sorted([rng.randrange(-4,5),rng.randrange(-4,5)])
     b.extend([F(x,8),F(y,8)])
    boxes.append(tuple(b))
   boxes+=boxes[:1]  # duplicates; input also permits degeneracy and crossing the origin
   exact=grid_volume(boxes);actual=volume([tuple(float(x) for x in b) for b in boxes])
   assert actual==float(exact),(n,boxes,actual,exact)
   count+=1
 result={'tests':count,'dimensions':[2,3,4],'status':'PASS','seconds':time.monotonic()-start,'arithmetic':'Fraction grid vs float recursive sweep; dyadic operands yield exact binary arithmetic at these sizes'}
 print(json.dumps(result,indent=2))
 (Path(__file__).resolve().parents[1]/'results/round002_union_tests.json').write_text(json.dumps(result,indent=2)+'\n')
if __name__=='__main__':main()
