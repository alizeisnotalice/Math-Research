#!/usr/bin/env python3
"""Event-driven discovery screen for the a2-only piecewise-affine union."""
from functools import lru_cache
from itertools import product
from math import prod
from pathlib import Path
import hashlib,json,time
TOL=1e-12

def constants():
 u=2**(1/6);v=3**(1/6);return (1+u**4/2,(2+u**5*v)/4,*((1+u)/2 for _ in range(4)))

def union_ratio(t):
 a=list(constants());a[1]*=t;raw=[]
 for ell in product(range(1,5),repeat=6):
  r=prod(ell)**(1/6)/2;b=[]
  for aj,k in zip(a,ell):
   lo,hi=(0.,min(aj/2,r-(k-1)*aj/2)) if k%2 else (max(0.,k*aj/2-r),aj/2)
   if hi<=lo+TOL:break
   b.extend((lo,hi))
  else:raw.append((ell,tuple(b)))
 U={b:e for e,b in raw};keep=[]
 for b,e in U.items():
  if any(c!=b and all(c[2*j]<=b[2*j]+TOL and b[2*j+1]<=c[2*j+1]+TOL for j in range(6)) for c in U):continue
  keep.append((e,b))
 @lru_cache(None)
 def rec(B):
  if not B:return 0.
  if len(B[0])==2:
   S=sorted(B);x,y=S[0];z=0.
   for a,b in S[1:]:
    if a>y+TOL:z,x,y=z+y-x,a,b
    elif b>y:y=b
   return z+y-x
  z=0.;cuts=sorted({x for b in B for x in b[:2]})
  for x,y in zip(cuts,cuts[1:]):
   C=tuple(sorted({b[2:] for b in B if b[0]<=x+TOL and y<=b[1]+TOL}))
   z+=(y-x)*rec(C)
  return z
 return {'t':t,'raw_boxes':len(raw),'kept_boxes':len(keep),'ratio_float':64*rec(tuple(sorted(b for _,b in keep)))}

def events(t0):
 a=constants();cand=[]
 for ell in product(range(1,5),repeat=6):
  r=prod(ell)**(1/6)/2
  # Only retain types whose other five coordinate intersections can exist.
  viable=True
  for j in range(6):
   if j==1:continue
   aj=a[j];k=ell[j];lo,hi=(0.,min(aj/2,r-(k-1)*aj/2)) if k%2 else (max(0.,k*aj/2-r),aj/2)
   if hi<=lo+TOL:viable=False;break
  if not viable:continue
  q=ell[1];s=a[1]/2
  cand.extend([(0.,0.),(s,0.),(-(q-1)*s,r),(q*s,-r)])
 vals=[]
 for i,(p,x) in enumerate(cand):
  for q,y in cand[i+1:]:
   if abs(p-q)<=TOL:continue
   z=(y-x)/(p-q)
   if .5<z<1.5:vals.append(z)
 vals.sort();unique=[]
 for z in vals:
  if not unique or abs(z-unique[-1])>1e-10:unique.append(z)
 below=[z for z in unique if z<t0-1e-10][-8:];above=[z for z in unique if z>t0+1e-10][:8]
 return below+above,{'all_event_count':len(unique),'below_count':len([z for z in unique if z<t0]),'above_count':len([z for z in unique if z>t0])}

def main():
 start=time.monotonic();t0=1.0016;pts,counts=events(t0);out=[union_ratio(t) for t in pts];ref=2.06563255516386
 for x in out:x['margin_over_k16_display_upper']=x['ratio_float']-ref
 winner=max(out,key=lambda x:x['ratio_float']) if out else None
 result={'status':'completed','scope':'event-driven numerical discovery only','t0':t0,'event_inventory':counts,'evaluated_event_count':len(out),'events':out,'winner':winner,'promotion_required':bool(winner and winner['margin_over_k16_display_upper']>1e-7),'elapsed_seconds':time.monotonic()-start,'code_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
 Path(__file__).with_name('second_spacing_breakpoint_results.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
if __name__=='__main__':main()
