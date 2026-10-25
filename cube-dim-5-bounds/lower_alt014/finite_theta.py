#!/usr/bin/env python3
"""Exact finite two-axis parity certificate. Only lower_alt013 is written."""
from fractions import Fraction as Q
from functools import lru_cache
from itertools import product
from math import prod
from pathlib import Path
import time,json,sys,hashlib,platform
HERE=Path(__file__).resolve().parent;ROOT=HERE.parent
D=1<<60

def radius(m,outer=False):
 t=m.numerator*D**5;d=32*m.denominator;lo,hi=0,1
 while d*hi**5<=t:hi*=2
 while hi-lo>1:
  z=(hi+lo)//2
  if d*z**5<=t:lo=z
  else:hi=z
 if d*lo**5==t:lo-=1
 assert lo>0 and d*lo**5<t
 return Q(lo+1 if outer else lo,D)

def run(spec,outer=False):
 start=time.monotonic();nodes=0
 sys.path.insert(0,str(ROOT/'lower'))
 from certify_rectangular import field
 F,u,v,bounds,_,_=field()
 ae=(1+u**3/2,)+((1+u)/2,)*4
 isolations=tuple(bounds(x.a) for x in ae)
 ai=tuple((Q(lo.numerator*D//lo.denominator,D),Q(-((-hi.numerator*D)//hi.denominator),D)) for lo,hi in isolations)
 a=tuple(x[0] for x in ai)
 def point(k,j,want_upper):
  low,high=ai[k]
  return j*(high if (want_upper == (j>=0)) else low)
 axes=tuple(spec['axes']);theta=Q(spec['theta']);dom=tuple(ai[k][1 if outer else 0]/(1 if k in axes else 2) for k in range(5))
 raw=[]; rs={}; enum=0; checks=0
 def cap():
  if time.monotonic()-start>120:raise TimeoutError('120 second cap')
 for ls in product(range(1,5),repeat=5):
  cap();base=prod(ls)
  # Every possible start on active axes is covered; r <= 2 for mass <= 1024.
  # Domain is [0,a], a>=1, so a positive receiver intersection requires -2<j<3.
  for js in product(range(-3,4),repeat=2):
   enum+=1
   alt=((-1 if (sum(js)%2) else 1)*prod(ls[k] for k in range(5) if k not in axes)) if all(ls[k]%2 for k in axes) else 0
   mass=Q(base)+theta*alt
   if mass<=0:continue
   if mass not in rs:rs[mass]=radius(mass,outer)
   r=rs[mass];box=[];starts=[]
   for k,(aa,l) in enumerate(zip(a,ls)):
    j=js[axes.index(k)] if k in axes else (-(l-1)//2 if l%2 else 1-l//2)
    domainbound=ai[k][1 if outer else 0]/(1 if k in axes else 2)
    lo=max(Q(0),point(k,j+l-1,not outer)-r)
    hi=min(domainbound,point(k,j,outer)+r)
    if hi<=lo:break
    box.extend((lo,hi));starts.append(j)
   else:
    # Distinct site mass checked by explicit parity counts, independent of formula.
    nplus=nminus=0
    for x,y in product(range(starts[axes[0]],starts[axes[0]]+ls[axes[0]]),range(starts[axes[1]],starts[axes[1]]+ls[axes[1]])):
     if (x+y)%2:nminus+=1
     else:nplus+=1
    direct=((1+theta)*nplus+(1-theta)*nminus)*prod(ls[k] for k in range(5) if k not in axes)
    assert direct==mass
    assert (32*r**5>=mass) if outer else (32*r**5<mass)
    # One common radius contains every selected block atom for every receiver.
    assert all(point(k,starts[k]+ls[k]-1,not outer)-r<=box[2*k] and box[2*k+1]<=point(k,starts[k],outer)+r for k in range(5))
    checks+=1;raw.append((tuple(box),{'starts':starts,'lengths':ls,'mass':str(mass),'radius':str(r),'distinct_positive_atoms':(nplus+(nminus if theta<1 else 0))*prod(ls[k] for k in range(5) if k not in axes)}))
 unique={b:meta for b,meta in raw};ret=[]
 for b in unique:
  cap()
  if not any(o!=b and all(o[2*k]<=b[2*k] and b[2*k+1]<=o[2*k+1] for k in range(5)) for o in unique):ret.append(b)
 @lru_cache(None)
 def sweep(boxes):
  nonlocal nodes
  nodes+=1;cap()
  if nodes>200000:raise RuntimeError('node cap')
  if not boxes:return Q(0)
  if len(boxes[0])==2:
   total=Q(0);lo=hi=None
   for x,y in sorted(boxes):
    if lo is None:lo,hi=x,y
    elif x>hi:total+=hi-lo;lo,hi=x,y
    else:hi=max(hi,y)
   return total+hi-lo
  cuts=sorted({x for b in boxes for x in b[:2]});total=Q(0)
  for x,y in zip(cuts,cuts[1:]):
   ss=tuple(sorted({b[2:] for b in boxes if b[0]<=x and y<=b[1]}))
   total+=(y-x)*sweep(ss)
  return total
 vol=sweep(tuple(sorted(ret)))
 seg=[]
 for k in range(5):
  cuts=sorted({b[2*k+d] for b in ret for d in (0,1)});ss=[]
  for x,y in zip(cuts,cuts[1:]):
   mask=sum(1<<i for i,b in enumerate(ret) if b[2*k]<=x and y<=b[2*k+1])
   if mask:ss.append((y-x,mask))
  seg.append(ss)
 @lru_cache(None)
 def masks(k,mask):
  cap()
  if k==5:return Q(1)
  return sum((width*masks(k+1,mask&mm) for width,mm in seg[k] if mask&mm),Q(0))
 vol2=masks(0,(1<<len(ret))-1);assert vol==vol2
 ratio=8*vol
 C=json.loads((ROOT/'lower_round003/a2_equals_a3_results.json').read_text());clo,chi=map(Q,C['ratio_interval'])
 out={'radius_mode':'outer' if outer else 'inner','status':'STRICT_BELOW_C_finite_family' if ratio<clo else 'STRICT_ABOVE_C' if ratio>chi else 'UNDECIDED_C','specification':spec,'exact_spacing_intervals':[[str(q) for q in pair] for pair in ai],'domain_upper':list(map(str,dom)),'period_mass':'4','normalization':'8*exact_union_volume','finite_lengths':'1..4 each coordinate','active_start_range':[-3,3],'enumerated_length_start_pairs':enum,'raw_feasible_boxes':len(raw),'unique_boxes':len(unique),'retained_boxes':len(ret),'direct_distinct_mass_checks':checks,'common_radius_containment_checks':checks,'radius_by_mass':{str(k):str(v) for k,v in rs.items()},'union_engines_agree':True,'volume':str(vol),'ratio':str(ratio),'ratio_display_only':float(ratio),'gap_to_C_lower':str(clo-ratio),'gap_display_only':float(clo-ratio),'C_source':'lower_round003/a2_equals_a3_results.json','C_source_sha256':hashlib.sha256((ROOT/'lower_round003/a2_equals_a3_results.json').read_bytes()).hexdigest(),'recursive_nodes':nodes,'segments':[len(s) for s in seg],'seconds':time.monotonic()-start,'python':sys.executable,'python_version':platform.python_version(),'witnesses':[dict(unique[b],box=list(map(str,b))) for b in ret],'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
 return out
if __name__=='__main__':
 reg=json.loads((HERE/'registration.json').read_text());spec=next(x for x in reg['inputs'] if x['name']==sys.argv[1])
 outer=len(sys.argv)>2 and sys.argv[2]=='outer'
 try:out=run(spec,outer)
 except (RuntimeError,TimeoutError) as e:out={'status':'COMPUTE_LIMIT','specification':spec,'reason':str(e)}
 (HERE/(spec['name']+('_outer' if outer else '_inner')+'_theta_results.json')).write_text(json.dumps(out,indent=2)+'\n')
 print(json.dumps({k:v for k,v in out.items() if k not in ['witnesses','radius_by_mass','gap_to_C_lower','ratio','volume']},indent=2))
