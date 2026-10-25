#!/usr/bin/env python3
"""Exact eight-candidate algebraic-event batch, preregistered in registration.json."""
from fractions import Fraction as Q
from functools import lru_cache
from itertools import product, combinations
from math import prod
from pathlib import Path
import hashlib,json,sys,time
HERE=Path(__file__).resolve().parent; ROOT=HERE.parent
sys.path.insert(0,str(ROOT/'lower'))
from certify_rectangular import field

def sha(path): return hashlib.sha256(path.read_bytes()).hexdigest()
def coefficients(x): return {f'{i//5},{i%5}':str(c) for i,c in enumerate(x.a) if c}
def decode(F,data):
 a=[Q(0)]*25
 for k,c in data.items():
  i,j=map(int,k.split(','));a[5*i+j]=Q(c)
 return F(a=a)
def root(n,u,v):
 e2=e3=0
 while n%2==0:n//=2;e2+=1
 while n%3==0:n//=3;e3+=1
 assert n==1
 return u**e2*v**e3

def event(F,u,v,bounds,base,direction):
 start=time.monotonic(); pairs=0;winner=F(Q(1,64)); selected=None;counts=[]
 for j,(aa,dj) in enumerate(zip(base,direction)):
  lines={(F(0),Q(0)),(aa/2,Q(dj,2))}
  for ls in product(range(1,5),repeat=5):
   ell=ls[j];r=root(prod(ls),u,v)/2
   lines.add((r-(ell-1)*aa/2,-Q((ell-1)*dj,2)))
   lines.add((ell*aa/2-r,Q(ell*dj,2)))
  counts.append(len(lines))
  for (x,s),(y,z) in combinations(sorted(lines,key=lambda q:(q[0].a,q[1])),2):
   pairs+=1
   if pairs>1000000 or time.monotonic()-start>120:raise TimeoutError('event selection cap')
   if s==z:continue
   cross=(y-x)/(s-z)
   if F(0)<cross and cross<=winner:
    winner=cross;selected={'axis':j,'left_intercept':coefficients(x),'left_slope':str(s),'right_intercept':coefficients(y),'right_slope':str(z)}
 lo,hi=bounds(winner.a)
 return winner,{'event_coefficients':coefficients(winner),'event_interval':[str(lo),str(hi)],'display_only':float(lo),'event_pair':selected,'raw_endpoint_counts':counts,'pair_tests':pairs,'wall_seconds':time.monotonic()-start}

def calculate(F,u,v,bounds,C,base,direction,t):
 started=time.monotonic();a=tuple(aa+d*t for aa,d in zip(base,direction));raw=[]; nodes=0
 def check():
  if time.monotonic()-started>120:raise TimeoutError('candidate 120-second cap')
 for ls in product(range(1,5),repeat=5):
  r=root(prod(ls),u,v)/2;box=[]
  for aa,ell in zip(a,ls):
   if ell%2:lo,hi=F(0),min(aa/2,r-(ell-1)*aa/2)
   else:lo,hi=max(F(0),ell*aa/2-r),aa/2
   if hi<=lo:break
   box.extend((lo,hi))
  else:raw.append((ls,tuple(box)))
 unique={b:ls for ls,b in raw};kept=[]
 for b,ls in unique.items():
  check()
  if not any(o!=b and all(o[2*j]<=b[2*j] and b[2*j+1]<=o[2*j+1] for j in range(5)) for o in unique):kept.append((ls,b))
 retained=tuple(sorted(b for ls,b in kept))
 @lru_cache(None)
 def rec(boxes):
  nonlocal nodes
  nodes+=1;check()
  if nodes>200000:raise RuntimeError('200000-node cap')
  if not boxes:return F(0)
  if len(boxes[0])==2:
   out=F(0);left=right=None
   for lo,hi in sorted(boxes):
    if left is None:left,right=lo,hi
    elif right<lo:out+=right-left;left,right=lo,hi
    elif right<hi:right=hi
   return out+right-left
  cuts=sorted({x for b in boxes for x in b[:2]});out=F(0)
  for lo,hi in zip(cuts,cuts[1:]):out+=(hi-lo)*rec(tuple(sorted({b[2:] for b in boxes if b[0]<=lo and hi<=b[1]})))
  return out
 volume=rec(retained);segments=[]
 for j in range(5):
  cuts=sorted({b[2*j+k] for b in retained for k in (0,1)});axis=[]
  for lo,hi in zip(cuts,cuts[1:]):
   mask=sum(1<<i for i,b in enumerate(retained) if b[2*j]<=lo and hi<=b[2*j+1])
   if mask:axis.append((hi-lo,mask))
  segments.append(axis)
 @lru_cache(None)
 def grid(j,mask):
  check()
  if j==5:return F(1)
  out=F(0)
  for width,nxt in segments[j]:
   both=nxt if mask<0 else mask&nxt
   if both:out+=width*grid(j+1,both)
  return out
 assert volume==grid(0,-1)
 # Verify actual unit atoms and a single radius for every retained block.
 witnesses=[]
 for ls,b in kept:
  N=prod(ls);r=root(N,u,v)/2
  assert 32*r**5==F(N)
  block=[list(range(-(ell//2),ell//2+1)) if ell%2 else list(range(1-ell//2,1+ell//2)) for ell in ls]
  atoms=list(product(*block));assert len(atoms)==N
  for j in range(5):
   lo,hi=b[2*j:2*j+2];mid=(lo+hi)/2
   assert F(0)<=lo and hi<=a[j]/2 and lo<hi
   for index in block[j]:
    atom=index*a[j]
    assert atom-r<=lo and hi<=atom+r
    assert atom-r<mid and mid<atom+r
  witnesses.append({'lengths':ls,'mass':N,'common_radius_coefficients':coefficients(r),'atom_index_ranges':[[q[0],q[-1]] for q in block],
   'atom_index_list_sha256':hashlib.sha256(json.dumps(atoms,separators=(',',':')).encode()).hexdigest(),'box_coefficients':[coefficients(x) for x in b]})
 ratio=32*volume;diff=ratio-C;rl,rh=bounds(ratio.a);dl,dh=bounds(diff.a)
 status='PASS_exact_finite_model_improvement' if dl>0 else ('NO_IMPROVEMENT_exact' if dh<=0 else 'UNRESOLVED_comparison')
 return {'status':status,'spacing_coefficients':[coefficients(x) for x in a],'ratio_coefficients':coefficients(ratio),'ratio_interval':[str(rl),str(rh)],'display_only':float(rl),
  'difference_from_exact_C_coefficients':coefficients(diff),'difference_from_exact_C_interval':[str(dl),str(dh)],'difference_display_only':float(dl),'strictly_above_exact_C':dl>0,
  'union_engines_agree':True,'raw_boxes':len(raw),'retained_boxes':len(kept),'coordinate_segment_counts':[len(q) for q in segments],'recursive_nodes':nodes,'wall_seconds':time.monotonic()-started,
  'real_atom_same_radius_containment_verified':True,'witnesses':witnesses}

def main():
 reg=json.loads((HERE/'registration.json').read_text());F,u,v,bounds,roots,den=field()
 Csource=ROOT/'lower_round003/a2_equals_a3_results.json';C=decode(F,json.loads(Csource.read_text())['coefficients'])
 base=(1+u**3/2,)+((1+u)/2,)*4
 result={'registration_sha256':sha(HERE/'registration.json'),'source_sha256':sha(Path(__file__)),'field_source_sha256':sha(ROOT/'lower/certify_rectangular.py'),'C_source_sha256':sha(Csource),'C_coefficients':coefficients(C),'root_lower_numerators':roots,'root_denominator':str(den),'candidates':[]}
 for base_d in reg['unoriented_directions']:
  for sign in (1,-1):
   d=tuple(sign*x for x in base_d);row={'direction':d}
   try:
    t,ev=event(F,u,v,bounds,base,d);row['event_selection']=ev;row.update(calculate(F,u,v,bounds,C,base,d,t))
   except (RuntimeError,TimeoutError) as e:row.update(status='COMPUTE_LIMIT',reason=str(e))
   result['candidates'].append(row)
   (HERE/'results.json').write_text(json.dumps(result,indent=2)+'\n')
   print(json.dumps({k:row[k] for k in ('direction','status','display_only','difference_display_only','raw_boxes','retained_boxes','wall_seconds') if k in row}),flush=True)
 result['status']='completed_8_preregistered_candidates' if len(result['candidates'])==8 and all(x['status']!='COMPUTE_LIMIT' for x in result['candidates']) else 'incomplete'
 (HERE/'results.json').write_text(json.dumps(result,indent=2)+'\n')
if __name__=='__main__':main()
