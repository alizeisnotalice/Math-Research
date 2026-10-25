#!/usr/bin/env python3
"""One rational rotation, exact algebraic-spacing finite-family enclosure."""
from fractions import Fraction as Q
from functools import lru_cache
from itertools import product
from math import prod
from pathlib import Path
from decimal import Decimal,localcontext
import hashlib,json,time
HERE=Path(__file__).resolve().parent;ROOT=HERE.parent
DEN=2**80;GRID=2**40;P,R,QC=4,3,5;CAP_SECONDS=120;CAP_NODES=200000
@lru_cache(None)
def side_bounds(n):
 lo,hi=0,8*DEN
 while lo+1<hi:
  mid=(lo+hi)//2
  if mid**5<=n*DEN**5:lo=mid
  else:hi=mid
 if lo**5==n*DEN**5:lo-=1
 a,b=Q(lo,DEN),Q(lo+1,DEN);assert a**5<n<=b**5
 return a,b

def indices(l):return tuple(range(-(l-1)//2,(l-1)//2+1)) if l%2 else tuple(range(1-l//2,1+l//2))
def interval(c,ab):return (c*ab[0],c*ab[1]) if c>=0 else (c*ab[1],c*ab[0])
def down(x):return Q((x*GRID).numerator//(x*GRID).denominator,GRID)
def up(x):return -down(-x)

def run():
 started=time.monotonic();nodes=0;stage='setup';configs=0
 def budget(node=False):
  nonlocal nodes
  if node:nodes+=1
  if nodes>CAP_NODES:raise RuntimeError('node cap at '+stage)
  if time.monotonic()-started>CAP_SECONDS:raise TimeoutError('wall cap at '+stage)
 assert P*P+R*R==QC*QC
 reps=[]
 for U,V in product(range(QC*QC),repeat=2):
  if (P*U+R*V)%(QC*QC)==0 and (-R*U+P*V)%(QC*QC)==0:
   reps.append((U,V,(P*U+R*V)//(QC*QC),(-R*U+P*V)//(QC*QC)))
 assert len(reps)==QC*QC==len({(U,V) for U,V,_,_ in reps})
 starts=[]
 for U,V,m,n in reps:
  for h,k in product((-1,0,1),repeat=2):starts.append((m+P*h+R*k,n-R*h+P*k))
 assert len(starts)==len(set(starts))==QC*QC*9
 ul,uh=side_bounds(2);aa=(1+ul**3/2,1+uh**3/2);ss=((1+ul)/2,(1+uh)/2)
 domains=[(aa[0]/2,ss[0]/2,ss[0]/2,QC*ss[0],QC*ss[0]),(aa[1]/2,ss[1]/2,ss[1]/2,QC*ss[1],QC*ss[1])]
 raw=[{},{}];stage='geometry'
 for ls in product(range(1,5),repeat=5):
  budget();n=prod(ls);slo,shi=side_bounds(n);rlo,rhi=slo/2,shi/2
  l1,l2,l3,l4,l5=ls
  fixed=[]
  for ell,space in ((l1,aa),(l4,ss),(l5,ss)):
   ii=indices(ell);fixed.append((interval(Q(min(ii)),space),interval(Q(max(ii)),space)))
  for m,k in starts:
   configs+=1
   cxmin=Q(P*m-R*(k+l3-1),QC);cxmax=Q(P*(m+l2-1)-R*k,QC)
   cymin=Q(R*m+P*k,QC);cymax=Q(R*(m+l2-1)+P*(k+l3-1),QC)
   extremes=fixed+[(interval(cxmin,ss),interval(cxmax,ss)),(interval(cymin,ss),interval(cymax,ss))]
   for mode in (0,1):
    box=[]
    for q,(mn,mx) in enumerate(extremes):
     if mode==0:lo=up(max(Q(0),mx[1]-rlo));hi=down(min(domains[0][q],mn[0]+rlo))
     else:lo=down(max(Q(0),mx[0]-rhi));hi=up(min(domains[1][q],mn[1]+rhi))
     if hi<=lo:break
     box.extend((lo,hi))
    else:raw[mode][tuple(box)]=(ls,m,k,n)
  budget()
 stage='dominance';retained=[]
 for data in raw:
  boxes=list(data);kept=[]
  for b in boxes:
   budget()
   if not any(o!=b and all(o[2*q]<=b[2*q] and b[2*q+1]<=o[2*q+1] for q in range(5)) for o in boxes):kept.append(b)
  retained.append(tuple(sorted(kept)))
 stage='direct_geometry_guard';guards=0;max_n=0
 for b in retained[0]:
  ls,m,k,n=raw[0][b];l1,l2,l3,l4,l5=ls;rlo=side_bounds(n)[0]/2
  pts=list(product(indices(l1),range(m,m+l2),range(k,k+l3),indices(l4),indices(l5)))
  assert len(pts)==len(set(pts))==n;max_n=max(max_n,n)
  coeffs=[(Q(j1),Q(j4),Q(j5),Q(P*mm-R*kk,QC),Q(R*mm+P*kk,QC)) for j1,mm,kk,j4,j5 in pts]
  for q,space in enumerate((aa,ss,ss,ss,ss)):
   mn=interval(min(p[q] for p in coeffs),space);mx=interval(max(p[q] for p in coeffs),space)
   assert mx[1]-rlo<=b[2*q] and b[2*q+1]<=mn[0]+rlo
   guards+=2
  assert (2*rlo)**5<n;budget()
 def union(boxes,label):
  nonlocal stage
  stage=label+'_sweep'
  @lru_cache(None)
  def sweep(bs):
   budget(True)
   if not bs:return Q(0)
   if len(bs[0])==2:
    out=Q(0);lo=hi=None
    for x,y in sorted(bs):
     if lo is None:lo,hi=x,y
     elif hi<x:out+=hi-lo;lo,hi=x,y
     else:hi=max(hi,y)
    return out+hi-lo
   cuts=sorted({x for b in bs for x in b[:2]});out=Q(0)
   for lo,hi in zip(cuts,cuts[1:]):out+=(hi-lo)*sweep(tuple(sorted({b[2:] for b in bs if b[0]<=lo and hi<=b[1]})))
   return out
  answer=sweep(boxes);segments=[]
  for q in range(5):
   budget();cuts=sorted({b[2*q+k] for b in boxes for k in (0,1)});axis=[]
   for lo,hi in zip(cuts,cuts[1:]):
    mask=sum(1<<i for i,b in enumerate(boxes) if b[2*q]<=lo and hi<=b[2*q+1])
    if mask:axis.append((hi-lo,mask))
   segments.append(axis)
  stage=label+'_mask'
  @lru_cache(None)
  def grid(q,mask):
   budget(True)
   if q==5:return Q(1)
   out=Q(0)
   for d,other in segments[q]:
    nxt=other if mask<0 else mask&other
    if nxt:out+=d*grid(q+1,nxt)
   return out
  assert grid(0,-1)==answer
  return answer,[len(axis) for axis in segments]
 vlo,slo=union(retained[0],'inner');vhi,shi=union(retained[1],'outer');lo,hi=Q(8,QC*QC)*vlo,Q(8,QC*QC)*vhi;assert lo<=hi
 cp=ROOT/'lower_round003/a2_equals_a3_results.json';c=json.loads(cp.read_text());cl,ch=map(Q,c['ratio_interval'])
 with localcontext() as ctx:
  ctx.prec=60;display=[str(Decimal(q.numerator)/Decimal(q.denominator)) for q in (lo,hi)]
 return {'status':'PASS_strict_inner_above_full_C_pending_audit' if lo>ch else ('NO_IMPROVEMENT_exact_finite_family' if hi<cl else 'UNRESOLVED_interval_overlap'),
  'rotation':{'cos':str(Q(P,QC)),'sin':str(Q(R,QC)),'supercell_atoms':QC*QC},'spacings':'exact a=1+u^3/2, s=(1+u)/2, u^5=2; no spacing substitution','representative_count':len(reps),'distinct_representatives_and_starts':True,'start_translations':[-1,0,1],'configurations':configs,'raw_boxes_inner_outer':[len(x) for x in raw],'retained_boxes_inner_outer':[len(x) for x in retained],'direct_selected_atom_distinct_guard':True,'direct_geometry_endpoint_inequalities':guards,'max_atoms_retained':max_n,'normalization':'8*absolute_union_volume/25','union_engines_agree_inner_and_outer':True,'ratio_interval':[str(lo),str(hi)],'decimal_display_only':display,'inner_above_full_C_upper':lo>ch,'outer_below_full_C_lower':hi<cl,'coordinate_segments_inner_outer':[slo,shi],'root_denominator':str(DEN),'coordinate_rounding_grid':str(GRID),'nodes':nodes,'wall_seconds':time.monotonic()-started,'limits':{'seconds':CAP_SECONDS,'nodes':CAP_NODES},'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'C_source_sha256':hashlib.sha256(cp.read_bytes()).hexdigest()}
def main():
 started=time.monotonic()
 try:out=run()
 except (RuntimeError,TimeoutError) as e:out={'status':'COMPUTE_LIMIT','reason':str(e),'wall_seconds':time.monotonic()-started,'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
 (HERE/'rotation_4_3_5_results.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
if __name__=='__main__':main()
