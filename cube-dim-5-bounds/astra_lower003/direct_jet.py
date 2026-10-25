#!/usr/bin/env python3
"""Direct mixed-point receiver jets at a previously selected free-h event."""
from pathlib import Path
from functools import lru_cache
from itertools import product
from fractions import Fraction as Q
from math import prod,comb
import json,hashlib,time,sys
HERE=Path(__file__).resolve().parent;ROOT=HERE.parent
sys.path.insert(0,str(ROOT/'lower'))
from certify_rectangular import field
F,u,v,bounds,_,_=field();ZERO=F(0);ONE=F(1)
def decode(d):
 a=[Q(0)]*25
 for k,c in d.items():i,j=map(int,k.split(','));a[5*i+j]=Q(c)
 return F(a=a)
def encode(x):return {f'{i//5},{i%5}':str(c) for i,c in enumerate(x.a) if c}
class J:
 __slots__=('x','d')
 def __init__(self,x=0,d=0):self.x=x if isinstance(x,F) else F(x);self.d=d if isinstance(d,F) else F(d)
 def __hash__(self):return hash((self.x,self.d))
 def __eq__(self,b):return self.x==b.x and self.d==b.d
 def __lt__(self,b):return self.x<b.x if self.x!=b.x else self.d<b.d
 def __le__(self,b):return self==b or self<b
 def __add__(self,b):return J(self.x+b.x,self.d+b.d)
 def __sub__(self,b):return J(self.x-b.x,self.d-b.d)
 def __mul__(self,b):return J(self.x*b.x,self.d*b.x+self.x*b.d)
JZERO=J();JONE=J(1)
A=(1+u**3/2,)+((1+u)/2,)*4
CELL=json.loads((HERE/'selected_cell.json').read_text());Z=tuple(decode(x) for x in CELL['shift_coefficients'])
RESULT=json.loads((HERE/'free_threshold_results.json').read_text());CAND=next(x for x in RESULT['candidates'] if x['event_index']==1)
HNUM=decode(CAND['h_numerator']);HDEN=decode(CAND['h_denominator'])
@lru_cache(None)
def radius(n,m):
 t=n;i=j=0
 while t%2==0:t//=2;i+=1
 while t%3==0:t//=3;j+=1
 assert t==1
 r=u**i*v**j/2;d=r*(HDEN*m-HNUM*n)/(5*n)
 assert 5*n*d==r*(HDEN*m-HNUM*n)
 return J(r,d)
def ref(index):
 f,p,q=list(product(range(2),range(3),range(3)))[index]
 return (-1 if f else 1,)+(-1,)*p+(1,)*(2-p)+(-1,)*q+(1,)*(2-q),comb(2,p)*comb(2,q)
def run(index):
 begin=time.monotonic();nodes=0;signs,multiplicity=ref(index);z=tuple(s*x for s,x in zip(signs,Z));raw={};patterns=0
 def cap():
  if time.monotonic()-begin>120 or nodes>200000:raise TimeoutError('registered direct jet cap')
 @lru_cache(None)
 def satellite_patterns(axis,lo,hi,r):
  cuts={lo,hi}
  for k in range(-4,5):
   point=k*A[axis]+z[axis]
   for e in (point-r,point+r):
    if lo<e and e<hi:cuts.add(e)
  cuts=sorted(cuts);out=set()
  for left,right in zip(cuts,cuts[1:]):
   x=(left+right)/2
   out.add(tuple(k for k in range(-4,5) if x-r<k*A[axis]+z[axis] and k*A[axis]+z[axis]<x+r))
  return tuple(sorted(out))
 for lengths in product(range(1,5),repeat=5):
  n=prod(lengths);r0=radius(n,0).x;primary=[];base=[]
  for axis,l in enumerate(lengths):
   start=-((l-1)//2);inds=tuple(range(start,start+l));primary.append(inds)
   lo=max(ZERO,max(inds)*A[axis]-r0);hi=min(A[axis]/2,min(inds)*A[axis]+r0)
   if hi<=lo:break
   base.extend((lo,hi))
  else:
   options=[satellite_patterns(axis,base[2*axis],base[2*axis+1],r0) for axis in range(5)]
   configurations={None}
   for sat in product(*options):configurations.add(sat if all(sat) else None)
   for satellite in configurations:
    patterns+=1;m=prod(map(len,satellite)) if satellite else 0;r=radius(n,m);b=[]
    for axis in range(5):
     coords=[k*A[axis] for k in primary[axis]]
     if satellite:coords.extend(k*A[axis]+z[axis] for k in satellite[axis])
     lo=max(JZERO,J(max(coords))-r);hi=min(J(A[axis]/2),J(min(coords))+r)
     if hi<=lo:break
     b.extend((lo,hi))
    else:raw[tuple(b)]=(tuple(primary),satellite,n,m)
   cap()
 kept=[]
 for b in raw:
  cap()
  if not any(o!=b and all(o[2*j]<=b[2*j] and b[2*j+1]<=o[2*j+1] for j in range(5)) for o in raw):kept.append(b)
 # Distinct orbit geometry: first-coordinate shifts are nonintegral periods.
 assert all(z[0]!=k*A[0] for k in range(-8,9))
 guards=0;atom_count=0
 for b in kept:
  primary,satellite,n,m=raw[b];r=radius(n,m)
  tags=[(0,idx) for idx in product(*primary)]
  if satellite:tags.extend((1,idx) for idx in product(*satellite))
  assert len(tags)==n+m==len(set(tags));atom_count+=len(tags)
  for orbit,coords in ((0,primary),(1,satellite)):
   if not coords:continue
   for axis in range(5):
    for k in coords[axis]:
     point=J(k*A[axis]+(z[axis] if orbit else ZERO))
     assert point-r<=b[2*axis] and b[2*axis+1]<=point+r;guards+=2
 @lru_cache(None)
 def volume(boxes):
  nonlocal nodes
  nodes+=1;cap()
  if not boxes:return JZERO
  cuts=sorted({x for b in boxes for x in b[:2]});total=JZERO
  for lo,hi in zip(cuts,cuts[1:]):
   sub=tuple(sorted({b[2:] for b in boxes if b[0]<=lo and hi<=b[1]}))
   if sub:total=total+(hi-lo)*(JONE if not sub[0] else volume(sub))
  return total
 V=volume(tuple(sorted(kept)));C=decode(json.loads((ROOT/'lower_round003/a2_equals_a3_results.json').read_text())['coefficients']);assert 32*V.x==C
 out={'status':'exact_direct_receiver_jet','reflection_index':index,'signs':signs,'multiplicity':multiplicity,'constant_halfcell_equals_C_over32':True,'mixed_volume_derivative':encode(V.d),'patterns':patterns,'raw_boxes':len(raw),'retained_boxes':len(kept),'tagged_atom_count':atom_count,'coordinate_endpoint_jet_guards':guards,'nodes':nodes,'seconds':time.monotonic()-begin,'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
 (HERE/f'direct_jet_{index:02}.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({k:v for k,v in out.items() if k!='mixed_volume_derivative'},indent=2))
def finish():
 value=ZERO;mult=0;guards=atoms=0
 for i in range(18):
  d=json.loads((HERE/f'direct_jet_{i:02}.json').read_text());mult+=d['multiplicity'];value+=d['multiplicity']*decode(d['mixed_volume_derivative']);guards+=d['coordinate_endpoint_jet_guards'];atoms+=d['tagged_atom_count']
 assert mult==32
 C=decode(json.loads((ROOT/'lower_round003/a2_equals_a3_results.json').read_text())['coefficients']);numerator=value+(HNUM-HDEN)*C+HDEN
 assert numerator==decode(CAND['D_numerator'])
 out={'status':'EXACT_FREE_THRESHOLD_CONTACT_FORMULA_MATCHES_DIRECT_MIXED_RECEIVER_JETS','event_index':1,'h_display_only':CAND['h_display_only'],'D_display_only':CAND['D_display_only'],'D_numerator':encode(numerator),'D_denominator':encode(HDEN),'normalization':'tau parameter: epsilon=h_den*tau,lambda=1+h_num*tau; add pure halo h_den and(h_num-h_den)*C','exact_coefficient_match':True,'reflection_multiplicity_sum':mult,'tagged_atoms_checked':atoms,'coordinate_endpoint_jet_guards':guards,'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
 (HERE/'direct_jet_verification.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({k:v for k,v in out.items() if k not in ('D_numerator','D_denominator')},indent=2))
if __name__=='__main__':
 if sys.argv[1]=='finish':finish()
 else:run(int(sys.argv[1]))
