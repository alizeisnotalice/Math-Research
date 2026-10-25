#!/usr/bin/env python3
"""Exact lower-bound certificates for centered cube maximal operators.
Pure Python standard library. No numerical library or floating point is used
in certificate generation, witness verification, union volume, or comparisons.
"""
from itertools import product
from functools import lru_cache
from fractions import Fraction
import json
from pathlib import Path
D=10**12
PARAMETERS={
 3:['1.630000','1.072357','1.129960'],
 4:['1.094603','1.053341','1.094603','1.707107'],
}
THRESHOLDS={3:Fraction(17939,10000),4:Fraction(19017,10000)}

def floor_root(n,k):
 lo=0;hi=1
 while hi**k<=n:hi*=2
 while hi-lo>1:
  m=(hi+lo)//2
  if m**k<=n:lo=m
  else:hi=m
 assert lo**k<=n<(lo+1)**k
 return lo

def build(d):
 A=[int(Fraction(s)*D) for s in PARAMETERS[d]]
 records=[]
 for ls in product(range(1,5),repeat=d):
  N=1
  for l in ls:N*=l
  t=floor_root(N*D**d,d)
  B=[]
  for a,l in zip(A,ls):
   # Quarter-cell coordinates are encoded with denominator 2D.
   if l%2:lo=0;hi=min(a,t-(l-1)*a)
   else:lo=max(0,l*a-t);hi=a
   if lo>=hi:break
   B.extend((lo,hi))
  else:
   # Verify every corner of this box lies in the witness center region.
   # l odd: rows -(l-1)/2,...,(l-1)/2.
   # l even: rows 1-l/2,...,l/2.
   starts=[-(l//2) if l%2 else 1-l//2 for l in ls]
   for j,(a,l,s) in enumerate(zip(A,ls,starts)):
    atom_lo=2*s*a;atom_hi=2*(s+l-1)*a
    assert max(abs(B[2*j]-atom_lo),abs(B[2*j]-atom_hi),
               abs(B[2*j+1]-atom_lo),abs(B[2*j+1]-atom_hi))<=t
   assert t**d<=N*D**d
   records.append({'counts':ls,'mass':N,'side_numerator':t,'box':tuple(B)})
 # Removing contained boxes is exact and only shortens the proof certificate.
 records.sort(key=lambda r:r['box'])
 keep=[]
 for i,r in enumerate(records):
  b=r['box'];discard=False
  for j,s in enumerate(records):
   if i==j:continue
   c=s['box']
   if all(c[2*k]<=b[2*k] and b[2*k+1]<=c[2*k+1] for k in range(d)):
    if b!=c or j<i:discard=True;break
  if not discard:keep.append(r)
 return A,keep

@lru_cache(maxsize=None)
def union_sweep(B):
 if not B:return 0
 if len(B[0])==2:
  ss=sorted(B);a,b=ss[0];v=0
  for c,e in ss[1:]:
   if c>b:v+=b-a;a,b=c,e
   else:b=max(b,e)
  return v+b-a
 cuts=sorted({x for b in B for x in b[:2]});v=0
 for a,b in zip(cuts,cuts[1:]):
  sub=tuple(sorted({c[2:] for c in B if c[0]<=a and b<=c[1]}))
  v+=(b-a)*union_sweep(sub)
 return v

def union_grid(B):
 """Independent coordinate-cell enumeration with bit-mask coverage."""
 d=len(B[0])//2
 cuts=[sorted({b[2*j+k] for b in B for k in (0,1)}) for j in range(d)]
 segments=[]
 for j,cs in enumerate(cuts):
  ss=[]
  for a,b in zip(cs,cs[1:]):
   mask=sum(1<<i for i,c in enumerate(B) if c[2*j]<=a and b<=c[2*j+1])
   if mask:ss.append((b-a,mask))
  segments.append(ss)
 def rec(j,mask):
  if j==d:return 1
  return sum(length*rec(j+1,mask&m) for length,m in segments[j] if mask&m)
 return rec(0,(1<<len(B))-1)

def run():
 results=[]
 for d in (3,4):
  A,recs=build(d);B=tuple(r['box'] for r in recs)
  V=union_sweep(tuple(sorted(B)))
  independent=union_grid(B)
  assert independent==V
  ratio=Fraction(V,D**d)
  assert ratio>THRESHOLDS[d]
  print(f'd={d}, {len(B)} witnesses, exact ratio={ratio}')
  print('decimal lower enclosure:',ratio.numerator*10**12//ratio.denominator,'/ 1000000000000')
  print('threshold comparison PASS; independent union verification PASS')
  results.append({'dimension':d,'spacing':PARAMETERS[d],'denominator':D,
    'coordinate_denominator':2*D,'witnesses':recs,'integer_union_volume':V,
    'lower_bound_numerator':ratio.numerator,'lower_bound_denominator':ratio.denominator,
    'certified_threshold':str(THRESHOLDS[d])})
 Path(__file__).with_name('certificates.json').write_text(json.dumps(results,indent=2))
if __name__=='__main__':run()
