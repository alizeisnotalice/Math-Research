"""Exact critical receiver endpoints and low-memory geometry for fixed022 seed."""
from pathlib import Path
from fractions import Fraction as F
from functools import lru_cache,cmp_to_key
from sympy import integer_nthroot
import json
ROOT=Path(__file__).resolve().parent
SEED=json.loads((ROOT/'round008_certificate_input.json').read_text())
P=tuple(map(F,SEED['periods']));PH=tuple(tuple(map(F,q)) for q in SEED['phases']);WTS=tuple(map(F,SEED['weights']));W=sum(WTS);D=10**36;Q=4*D

@lru_cache(None)
def rational_radius(mass):
 a,ea=integer_nthroot(mass.numerator,3);b,eb=integer_nthroot(mass.denominator,3)
 return F(int(a),2*int(b)) if ea and eb else None

@lru_cache(None)
def radius_interval(mass,digits=36):
 den=10**digits;q=int(integer_nthroot(mass.numerator*den**3//mass.denominator,3)[0]);assert q**3*mass.denominator<=mass.numerator*den**3<(q+1)**3*mass.denominator
 return F(q,2*den),F(q+1,2*den)

def expr(base,sign=0,mass=F(0)):
 if sign==0:return ('q',F(base))
 r=rational_radius(mass)
 return ('q',base+sign*r) if r is not None else ('r',F(base),int(sign),F(mass))

@lru_cache(None)
def bounds(e,digits=36):
 if e[0]=='q':return e[1],e[1]
 lo,hi=radius_interval(e[3],digits)
 return (e[1]+lo,e[1]+hi) if e[2]>0 else (e[1]-hi,e[1]-lo)

@lru_cache(None)
def compare(a,b):
 if a==b:return 0
 digits=36
 while True:
  al,ah=bounds(a,digits);bl,bh=bounds(b,digits)
  if ah<bl:return -1
  if bh<al:return 1
  # Distinct canonical pure-cubic affine endpoints cannot be equal: if
  # two signed cube roots differ by nonzero rational, a cubic expansion
  # gives a quadratic equation, forcing both roots rational (handled above).
  digits*=2

def midpoint_integer(e):
 lo,hi=bounds(e,36);q=(lo+hi)*Q/2;assert q.denominator==1
 return int(q)

def block_counts(key):
 lo=[key[3*b:3*b+3] for b in range(3)];hi=[key[9+3*b:12+3*b] for b in range(3)];count=[]
 for b in range(3):
  n=1
  for j in range(3):n*=max(0,hi[b][j]-lo[b][j]+1)
  count.append(n)
 return lo,hi,count

def receiver_box(key):
 lo,hi,n=block_counts(key);active=[b for b in range(3) if n[b]];mass=sum((n[b]*WTS[b] for b in active),F(0))
 if not active:return None
 box=[]
 for j in range(3):
  high=max(P[j]*(hi[b][j]+PH[b][j]) for b in active);low=min(P[j]*(lo[b][j]+PH[b][j]) for b in active);a=expr(high,-1,mass);b=expr(low,1,mass);zero=expr(0);edge=expr(P[j])
  if compare(a,zero)<0:a=zero
  if compare(b,edge)>0:b=edge
  if compare(a,b)>=0:return None
  box.append((a,b))
 return tuple(box)

def old_geometry():
 boxes=[];ids=[]
 for i,key in enumerate(SEED['keys']):
  b=receiver_box(key)
  if b is not None:boxes.append(b);ids.append(i)
 return boxes,ids

def exact_coverage_test(old_boxes,new_boxes):
 """Decide new finite union subset old, using ordered endpoints and z masks."""
 if not new_boxes:return {'positive_new_cell':False,'new_nonempty_boxes':0,'xy_checks':0}
 boxes=old_boxes+new_boxes;cuts=[sorted({e for b in boxes for e in b[j]},key=cmp_to_key(compare)) for j in range(3)];look=[{e:i for i,e in enumerate(c)} for c in cuts]
 def axis(j):
  events={}
  for i,b in enumerate(boxes):
   for e in b[j]:
    at=look[j][e];events[at]=events.get(at,0)^(1<<i)
  current=0;out=[]
  for q in range(len(cuts[j])-1):current^=events.get(q,0);out.append(current)
  return out
 masks=[axis(j) for j in range(3)];zspan=[]
 for b in boxes:
  lo,hi=look[2][b[2][0]],look[2][b[2][1]];zspan.append(((1<<(hi-lo))-1)<<lo)
 @lru_cache(None)
 def zu(mask):
  answer=0
  while mask:
   bit=mask&-mask;i=bit.bit_length()-1;answer|=zspan[i];mask^=bit
  return answer
 oldmask=(1<<len(old_boxes))-1;newmask=((1<<len(boxes))-1)^oldmask;checks=0
 for i,x in enumerate(masks[0]):
  for j,y in enumerate(masks[1]):
   active=x&y
   if not active&newmask:continue
   checks+=1;uncovered=zu(active&newmask)&~zu(active&oldmask)
   if uncovered:
    bit=uncovered&-uncovered;k=bit.bit_length()-1
    return {'positive_new_cell':True,'new_nonempty_boxes':len(new_boxes),'xy_checks':checks,'cell_indices':[i,j,k],'coordinate_cut_counts':[len(c) for c in cuts]}
 return {'positive_new_cell':False,'new_nonempty_boxes':len(new_boxes),'xy_checks':checks,'coordinate_cut_counts':[len(c) for c in cuts],'old_nonempty_boxes':len(old_boxes)}
