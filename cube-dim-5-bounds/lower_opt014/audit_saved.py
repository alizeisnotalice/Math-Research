#!/usr/bin/env python3
"""Read-only reconstruction of the saved event, atom and volume certificates."""
from fractions import Fraction as Q
from functools import lru_cache
from itertools import product
from math import prod
from pathlib import Path
import json,hashlib,sys,time
HERE=Path(__file__).resolve().parent;ROOT=HERE.parent
sys.path.insert(0,str(ROOT/'lower'))
from certify_rectangular import field

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 start=time.monotonic();F,u,v,bounds,roots,den=field()
 def dec(data):
  a=[Q(0)]*25
  for k,c in data.items():i,j=map(int,k.split(','));a[5*i+j]=Q(c)
  return F(a=a)
 result=json.loads((HERE/'exact_results.json').read_text());reg=json.loads((HERE/'registration.json').read_text())
 assert result['source_sha256']==sha(HERE/'certify_selected.py')
 assert result['registration_sha256']==sha(HERE/'registration.json')
 cp=ROOT/'lower_round003/a2_equals_a3_results.json'
 assert result['C_source_sha256']==sha(cp)
 C=dec(json.loads(cp.read_text())['coefficients']);assert C==dec(result['C_coefficients'])
 allspacings=set();rows=[]
 for i,row in enumerate(result['candidates']):
  assert time.monotonic()-start<120
  a=tuple(dec(q) for q in row['spacing_coefficients'])
  assert a not in allspacings;allspacings.add(a)
  def root(n):
   e2=e3=0
   while n%2==0:n//=2;e2+=1
   while n%3==0:n//=3;e3+=1
   assert n==1
   return u**e2*v**e3
  def coordinate(q):
   if 'rational' in q:return F(Q(q['rational']))
   x,sx,qx=q['line_left'];y,sy,qy=q['line_right']
   xx=F(0) if x==0 else sx*root(x)/2
   yy=F(0) if y==0 else sy*root(y)/2
   t=(yy-xx)/Q(qx-qy,2)
   assert xx+Q(qx,2)*t==yy+Q(qy,2)*t
   return t
  selected=row['selected_float_row'];b,c=coordinate(selected['b']),coordinate(selected['c'])
  assert a==(b,b,c,c,c)
  boxes=[];atoms_checked=0
  for w in row['witnesses']:
   ls=w['lengths'];N=prod(ls);assert N==w['mass']
   r=dec(w['common_radius_coefficients']);assert 32*r**5==F(N)
   b=tuple(dec(q) for q in w['box_coefficients']);boxes.append(b)
   indices=[list(range(lo,hi+1)) for lo,hi in w['atom_index_ranges']]
   atoms=list(product(*indices));assert len(atoms)==N
   assert hashlib.sha256(json.dumps(atoms,separators=(',',':')).encode()).hexdigest()==w['atom_index_list_sha256']
   for k in atoms:
    atoms_checked+=1
    for j in range(5):
     x=k[j]*a[j];lo,hi=b[2*j:2*j+2];mid=(lo+hi)/2
     assert x-r<=lo and hi<=x+r
     assert x-r<mid and mid<x+r
  # Sweep the five axes in reverse order to reconstruct saved volume.
  reversed_boxes=tuple(sorted(tuple(z for j in (4,3,2,1,0) for z in b[2*j:2*j+2]) for b in boxes))
  @lru_cache(None)
  def sweep(bs):
   assert time.monotonic()-start<120
   if not bs:return F(0)
   cuts=sorted({x for b in bs for x in b[:2]});answer=F(0)
   for lo,hi in zip(cuts,cuts[1:]):
    selected=tuple(sorted({b[2:] for b in bs if b[0]<=lo and hi<=b[1]}))
    if selected:answer+=(hi-lo)*(F(1) if len(bs[0])==2 else sweep(selected))
   return answer
  ratio=32*sweep(reversed_boxes);assert ratio==dec(row['ratio_coefficients'])
  diff=ratio-C;assert diff==dec(row['difference_from_exact_C_coefficients'])
  dl,dh=bounds(diff.a);assert dh<0
  assert [str(dl),str(dh)]==row['difference_from_exact_C_interval']
  rows.append({'candidate':i,'reverse_axis_union_matches':True,'actual_atoms_checked':atoms_checked,'complete_C_difference_upper_negative':True})
 out={'status':'PASS_saved_two_candidate_independent_reconstruction','rows':rows,'result_sha256':sha(HERE/'exact_results.json'),'audit_source_sha256':sha(Path(__file__)),'wall_seconds':time.monotonic()-start,'wall_cap_seconds':120}
 (HERE/'audit_results.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
if __name__=='__main__':main()
