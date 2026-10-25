#!/usr/bin/env python3
"""Bounded nonlocal 2+3 plane screen, with exact-reconstructible event provenance."""
from pathlib import Path
from fractions import Fraction as Q
from itertools import product,combinations
from functools import lru_cache
import numpy as np,json,time,hashlib
HERE=Path(__file__).resolve().parent
LS=np.array(list(product(range(1,5),repeat=5)),dtype=int)
N=np.prod(LS,axis=1);RAD=N.astype(float)**0.2/2
started=time.monotonic();evaluations=0

def check():
 if time.monotonic()-started>240:raise TimeoutError('240-second float cap')
def union_ratio(b,c):
 global evaluations
 evaluations+=1;check()
 if evaluations>80008:raise RuntimeError('union evaluation cap')
 a=np.array([b,b,c,c,c]);odd=LS%2==1
 low=np.where(odd,0,np.maximum(0,LS*a/2-RAD[:,None]))
 high=np.where(odd,np.minimum(a/2,RAD[:,None]-(LS-1)*a/2),a/2)
 valid=np.all(high>low,axis=1);raw=np.stack((low[valid],high[valid]),axis=2).reshape(-1,10)
 raw=np.unique(raw,axis=0)
 dom=np.all(raw[:,None,0::2]>=raw[None,:,0::2],axis=2)&np.all(raw[:,None,1::2]<=raw[None,:,1::2],axis=2)
 np.fill_diagonal(dom,False);raw=raw[~np.any(dom,axis=1)]
 nodes=0
 @lru_cache(None)
 def sweep(boxes):
  nonlocal nodes
  nodes+=1
  if nodes>20000:raise RuntimeError('float union node cap')
  if not boxes:return 0.
  if len(boxes[0])==2:
   ans=0.;left=right=None
   for lo,hi in sorted(boxes):
    if left is None:left,right=lo,hi
    elif right<lo:ans+=right-left;left,right=lo,hi
    elif right<hi:right=hi
   return ans+right-left
  cuts=sorted({q for box in boxes for q in box[:2]});ans=0.
  for lo,hi in zip(cuts,cuts[1:]):
   sub=tuple(sorted({box[2:] for box in boxes if box[0]<=lo and hi<=box[1]}))
   if sub:ans+=(hi-lo)*sweep(sub)
  return ans
 return 32*sweep(tuple(map(tuple,raw)))

def events():
 # intercept sign*root(N)/2, slope slope_num/2; root(0)=0.
 lines={(0,0,0),(0,0,1)}
 for ls,n in zip(LS,N):
  ell=int(ls[0]);n=int(n)
  lines.add((n,1,-(ell-1)));lines.add((n,-1,ell))
 lines=sorted(lines);seen={}
 def val(line):n,s,q=line;return s*n**0.2/2 if n else 0.
 for x,y in combinations(lines,2):
  check()
  if x[2]==y[2]:continue
  a=2*(val(y)-val(x))/(x[2]-y[2])
  if .8<=a<=2.25:seen.setdefault(round(a,12),{'display_only':a,'line_left':x,'line_right':y})
 if len(seen)>20000:raise RuntimeError('event cap')
 out=sorted(seen.values(),key=lambda x:x['display_only'])
 # Include exact rational domain endpoints for bounded coordinate scans.
 out=[{'display_only':.8,'rational':'4/5'}]+out+[{'display_only':2.25,'rational':'9/4'}]
 return out

def main():
 reg=json.loads((HERE/'registration.json').read_text());evs=events();rows=[]
 for i,(b,c) in enumerate(reg['seeds']):
  row={'type':'seed','index':i,'b':{'rational':b},'c':{'rational':c},'ratio_display_only':union_ratio(float(Q(b)),float(Q(c)))}
  rows.append(row);print(json.dumps(row),flush=True)
 current=max(rows,key=lambda x:x['ratio_display_only']).copy()
 for step,axis in enumerate(('b','c','b','c')):
  best=current.copy();start_eval=evaluations
  other='c' if axis=='b' else 'b';fixed=current[other]
  def value(q):return float(Q(q['rational'])) if 'rational'in q else q['display_only']
  for event in evs:
   b=value(event) if axis=='b' else value(fixed);c=value(event) if axis=='c' else value(fixed)
   ratio=union_ratio(b,c)
   if ratio>best['ratio_display_only']+1e-12:best={**current,axis:event,'ratio_display_only':ratio}
  best.update(type='coordinate_event_step',index=step,axis=axis,event_evaluations=evaluations-start_eval)
  rows.append(best.copy());current=best
  print(json.dumps({k:best[k] for k in ('type','index','axis','b','c','ratio_display_only','event_evaluations')}),flush=True)
  (HERE/'float_results.json').write_text(json.dumps({'status':'partial','rows':rows,'event_count':len(evs),'union_evaluations':evaluations},indent=2)+'\n')
 out={'status':'completed_bounded_numerical_screen','rows':rows,'event_count':len(evs),'union_evaluations':evaluations,'wall_seconds':time.monotonic()-started,
 'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'registration_sha256':hashlib.sha256((HERE/'registration.json').read_bytes()).hexdigest(),'scope':'float diagnostics only; event values reconstructed exactly before certification'}
 (HERE/'float_results.json').write_text(json.dumps(out,indent=2)+'\n')
 print(json.dumps({k:out[k] for k in ('status','event_count','union_evaluations','wall_seconds')}),flush=True)
if __name__=='__main__':main()
