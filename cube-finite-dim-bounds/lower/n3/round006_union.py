"""Low-memory full floating box-union integration by axis masks and z intervals."""
from functools import lru_cache
import numpy as np

def volume(boxes):
 boxes=np.asarray(boxes,dtype=float)
 if len(boxes)==0:return 0.
 def axis(j):
  events={}
  for i,b in enumerate(boxes):
   for x in b[j]:events[x]=events.get(x,0)^(1<<i)
  cuts=sorted(events);mask=0;weights={}
  for u,v in zip(cuts,cuts[1:]):
   mask^=events[u]
   if mask:weights[mask]=weights.get(mask,0.)+(v-u)
  return list(weights.items())
 xs=axis(0);ys=axis(1)
 @lru_cache(None)
 def zlength(mask):
  if not mask:return 0.
  spans=[]
  while mask:
   bit=mask&-mask;i=bit.bit_length()-1;spans.append(tuple(boxes[i,2]));mask^=bit
  spans.sort();a,b=spans[0];total=0.
  for c,d in spans[1:]:
   if c>b:total+=b-a;a,b=c,d
   else:b=max(b,d)
  return float(total+b-a)
 return float(sum(dx*sum(dy*zlength(x&y) for y,dy in ys if x&y) for x,dx in xs))

def selftest():
 from pathlib import Path
 import json
 from witness import union as grid_union
 rng=np.random.default_rng(20261025605);differences=[]
 for i in range(32):
  b=np.sort(rng.random((1+i%9,3,2)),axis=2);difference=abs(volume(b)-grid_union(b)[0]);assert difference<1e-12;differences.append(difference)
 out={'seed':20261025605,'cases':32,'all_passed':True,'largest_absolute_difference':max(differences),'method':'axis masks/cached z interval union vs original3D grid union'}
 (Path(__file__).resolve().parent/'round006_union_check.json').write_text(json.dumps(out,indent=2)+'\n');print(out)
if __name__=='__main__':selftest()
