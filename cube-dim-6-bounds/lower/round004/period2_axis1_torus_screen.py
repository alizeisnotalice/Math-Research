#!/usr/bin/env python3
"""Correct torus-cell screen for period-two masses in coordinate 1."""
from functools import lru_cache
from itertools import product
from math import prod,floor,ceil
from pathlib import Path
import hashlib,json,time
TOL=1e-12; CAP=5000
def spacing():
 u=2**(1/6);v=3**(1/6);return (1+u**4/2,(u**3+u**4*v**3)/4,*((1+u)/2 for _ in range(4)))
def torus_pieces(lo,hi,P):
 out=[]
 for q in range(floor(-hi/P)-1,ceil((P-lo)/P)+2):
  a,b=lo+q*P,hi+q*P
  x,y=max(0.,a),min(P,b)
  if y>x+TOL:out.append((x,y,q))
 return out
def boxes(delta):
 a=spacing();per=(2*a[0],*a[1:]);raw=[]
 for ell in product(range(1,5),repeat=6):
  rest=prod(ell[j] for j in range(1,6));k=ell[0]
  for parity in (0,1):
   selected=k if k%2==0 else k+(delta if parity==0 else -delta)
   r=(rest*selected)**(1/6)/2;arcs=[]
   for j,(aj,P) in enumerate(zip(a,per)):
    start=parity if j==0 else 0
    arcs.append(torus_pieces((start+ell[j]-1)*aj-r,start*aj+r,P))
   for pieces in product(*arcs):
    b=tuple(z for part in pieces for z in part[:2]);sh=[x[2] for x in pieces]
    # A period shift of the center corresponds to the same period-shifted
    # source block; axis 1 shifts by two and keeps its mass parity.
    for j,aj in enumerate(a):
     start=parity if j==0 else 0; jump=per[j]*sh[j]
     for c in b[2*j:2*j+2]:
      for z in (start*aj+jump,(start+ell[j]-1)*aj+jump):
       if abs(c-z)>r+3e-12: raise AssertionError((ell,parity,delta,j,c,z,r,sh,pieces))
    raw.append((ell,parity,rest*selected,b))
    if len(raw)>CAP:raise RuntimeError('raw-box cap exceeded')
 U={b:(e,p,m) for e,p,m,b in raw};keep=[]
 for b,item in U.items():
  if any(c!=b and all(c[2*j]<=b[2*j]+TOL and b[2*j+1]<=c[2*j+1]+TOL for j in range(6)) for c in U):continue
  keep.append((*item,b))
 return raw,keep
def union(B):
 @lru_cache(None)
 def rec(X):
  if not X:return 0.
  if len(X[0])==2:
   S=sorted(X);a,b=S[0];z=0.
   for c,d in S[1:]:
    if c>b+TOL:z,a,b=z+b-a,c,d
    elif d>b:b=d
   return z+b-a
  z=0.;cuts=sorted({x for b in X for x in b[:2]})
  for a,b in zip(cuts,cuts[1:]):
   Y=tuple(sorted({x[2:] for x in X if x[0]<=a+TOL and b<=x[1]+TOL}));z+=(b-a)*rec(Y)
  return z
 return rec(tuple(sorted(B)))
def evaluate(delta):
 raw,keep=boxes(delta);return {'delta':str(delta),'raw_boxes':len(raw),'kept_boxes':len(keep),'ratio_float':union(tuple(b for *_,b in keep))/2}
def main():
 t=time.monotonic();cal=evaluate(0.);pilot=evaluate(1/20);ref=2.06692382601089
 out={'status':'completed','scope':'torus geometry/discovery only; no C6 claim','axis_one_based':1,'calibration_delta_zero':cal,'pilot_delta_1_20':pilot,'pilot_margin_over_breakpoint_display_upper':pilot['ratio_float']-ref,'promotion_required':pilot['ratio_float']>ref+1e-7,'elapsed_seconds':time.monotonic()-t,'code_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
 Path(__file__).with_name('period2_axis1_torus_results.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
if __name__=='__main__':main()
