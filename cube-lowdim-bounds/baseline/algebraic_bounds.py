#!/usr/bin/env python3
"""Exact radical construction and finite box-union certificate.
All arithmetic and comparisons use fractions and rational isolating intervals.
No floating-point computation is used in the proof.
"""
from fractions import Fraction as Q
from functools import lru_cache
from itertools import product
import json
from pathlib import Path
from verify_bounds import floor_root


def calculate(d):
 D=10**45
 roots=[floor_root(n*D**d,d) for n in (2,3)]
 ul,vl=[Q(n,D) for n in roots];uh,vh=[Q(n+1,D) for n in roots]
 size=d*d
 @lru_cache(None)
 def add(a,b):return tuple(x+y for x,y in zip(a,b))
 @lru_cache(None)
 def neg(a):return tuple(-x for x in a)
 @lru_cache(None)
 def mul(a,b):
  z=[Q(0)]*size
  for i,x in enumerate(a):
   if not x:continue
   ei,ej=divmod(i,d)
   for j,y in enumerate(b):
    if not y:continue
    fi,fj=divmod(j,d);ii=ei+fi;jj=ej+fj
    z[(ii%d)*d+jj%d]+=x*y*2**(ii//d)*3**(jj//d)
  return tuple(z)
 @lru_cache(None)
 def bounds(a):
  low=Q(0);high=Q(0)
  for k,c in enumerate(a):
   if not c:continue
   i,j=divmod(k,d);p=ul**i*vl**j;q=uh**i*vh**j
   if c>0:low+=c*p;high+=c*q
   else:low+=c*q;high+=c*p
  return low,high
 @lru_cache(None)
 def positive(a):
  if not any(a):return False
  low,high=bounds(a)
  if low>0:return True
  if high<0:return False
  raise ArithmeticError('Isolating intervals need refinement')
 class F:
  __slots__=('a',)
  def __init__(self,n=0,a=None):
   if a is not None:self.a=tuple(a)
   else:self.a=(Q(n),)+(Q(0),)*(size-1)
  def __hash__(self):return hash(self.a)
  def __eq__(self,b):
   if not isinstance(b,F):b=F(b)
   return self.a==b.a
  def __add__(self,b):
   if not isinstance(b,F):b=F(b)
   return F(a=add(self.a,b.a))
  __radd__=__add__
  def __neg__(self):return F(a=neg(self.a))
  def __sub__(self,b):return self+-b if isinstance(b,F) else self+F(-b)
  def __rsub__(self,b):return F(b)+-self
  def __mul__(self,b):
   if not isinstance(b,F):b=F(b)
   return F(a=mul(self.a,b.a))
  __rmul__=__mul__
  def __truediv__(self,n):return F(a=tuple(x/Q(n) for x in self.a))
  def __pow__(self,n):
   z=F(1)
   for _ in range(n):z=z*self
   return z
  def __lt__(self,b):
   if not isinstance(b,F):b=F(b)
   return positive((b-self).a)
  def __le__(self,b):return self==b or self<b
  def __repr__(self):return str(self.a)
 u0=[Q(0)]*size;u0[d]=Q(1);u=F(a=u0)
 v0=[Q(0)]*size;v0[1]=Q(1);v=F(a=v0)
 assert u**d==F(2) and v**d==F(3)
 a=[1+u**(d-2)/2,(2+u**(d-1)*v)/4]+[(1+u)/2]*(d-2)
 def root(n):
  i=j=0
  while n%2==0:i+=1;n//=2
  while n%3==0:j+=1;n//=3
  assert n==1
  return u**i*v**j
 records=[]
 for ls in product(range(1,5),repeat=d):
  N=1
  for l in ls:N*=l
  side=root(N);assert side**d==F(N);r=side/2;B=[]
  for aa,l in zip(a,ls):
   if l%2:lo=F(0);hi=min(aa/2,r-(l-1)*aa/2)
   else:lo=max(F(0),l*aa/2-r);hi=aa/2
   if hi<=lo:break
   B.extend((lo,hi))
  else:
   starts=[-(l//2) if l%2 else 1-l//2 for l in ls]
   for j,(aa,l,s) in enumerate(zip(a,ls,starts)):
    for x in B[2*j:2*j+2]:
     for z in (s*aa,(s+l-1)*aa):
      assert -r<=x-z<=r
   records.append((ls,tuple(B)))
 uniq={b:ls for ls,b in records};recs=[]
 for b,ls in uniq.items():
  if any(c!=b and all(c[2*j]<=b[2*j] and b[2*j+1]<=c[2*j+1] for j in range(d)) for c in uniq):continue
  recs.append((ls,b))
 @lru_cache(None)
 def volume(B):
  if not B:return F(0)
  if len(B[0])==2:
   ss=sorted(B);a,b=ss[0];z=F(0)
   for c,e in ss[1:]:
    if b<c:z+=b-a;a,b=c,e
    else:b=max(b,e)
   return z+b-a
  cuts=sorted({x for b in B for x in b[:2]});z=F(0)
  for a,b in zip(cuts,cuts[1:]):
   sub=tuple(sorted({c[2:] for c in B if c[0]<=a and b<=c[1]}))
   z+=(b-a)*volume(sub)
  return z
 R=2**d*volume(tuple(sorted(b for ls,b in recs)))
 lo,hi=bounds(R.a);scale=10**14
 if d==3:assert lo>Q('1.7939332277947')
 if d==4:assert lo>Q('1.9017426111009')
 print('dimension',d,'witnesses',len(recs),flush=True)
 print('lower interval',lo.numerator*scale//lo.denominator,'/',scale,flush=True)
 print('upper interval',(hi.numerator*scale+hi.denominator-1)//hi.denominator,'/',scale,flush=True)
 coefs={f'{i},{j}':str(R.a[i*d+j]) for i in range(d) for j in range(d) if R.a[i*d+j]}
 print('polynomial coefficients',coefs,flush=True)
 return {'dimension':d,'parameter_formula':['1+u^(d-2)/2','(2+u^(d-1)*v)/4']+['(1+u)/2']*(d-2),
         'u_definition':f'2^(1/{d})','v_definition':f'3^(1/{d})','coefficients':coefs,
         'witness_counts':[ls for ls,b in recs],'lower_bound_interval':[str(lo),str(hi)],
         'root_denominator':D,'root_lower_numerators':roots}

if __name__=='__main__':
 result=[calculate(d) for d in (2,3,4)]
 Path(__file__).with_name('algebraic_certificates.json').write_text(json.dumps(result,indent=2))
