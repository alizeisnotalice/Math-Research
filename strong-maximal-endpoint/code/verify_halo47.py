"""Exact rational witness checker; independent mpmath budget diagnostic."""
from fractions import Fraction as F
from itertools import combinations
import json,sys,platform
from pathlib import Path
import mpmath as mp
mp.mp.dps=70
src=Path(sys.argv[1]) if len(sys.argv)>1 else Path(__file__).resolve().parents[1]/'rounds/r047/fullhalo47.json'
d=json.loads(src.read_text(),parse_float=str)
xs=list(map(F,d['receiver_endpoint_grid']['x']));ys=list(map(F,d['receiver_endpoint_grid']['y']))
cells=[(a,b,c,e) for a,b in zip(xs,xs[1:]) for c,e in zip(ys,ys[1:])]
def overlap(a,b,c,e):return max(F(0),min(b,e)-max(a,c))
def prim(s,eps,y):return F(35,32)*s-s**3/8+eps*y*(s/2+s*s)
def mass(rect,eps):
 a,b,c,e=rect;z=F(0)
 for l,r,y in [(F(0),F(1,2),1),(F(1,2),F(1),-1)]:
  height=overlap(c,e,l,r)
  hm=overlap(a,b,0,1)*(F(31,32)-eps*y/2)
  lo=max(a,F(-1));hi=min(b,F(0))
  if hi>lo:hm+=prim(-lo,eps,y)-prim(-hi,eps,y)
  z+=height*hm
 return z
out=[]
for row in d['rows']:
 eps=F(row['epsilon']); good=[];uc=set();margins=[]
 for (a,b),(c,e) in ((x,y) for x in combinations(xs,2) for y in combinations(ys,2)):
  area=(b-a)*(e-c);avg=mass((a,b,c,e),eps)/area
  if avg>1:
   good.append((a,b,c,e));margins.append(avg-1)
   uc.update(cell for cell in cells if a<=cell[0]<cell[1]<=b and c<=cell[2]<cell[3]<=e)
 area=sum(((b-a)*(e-c) for a,b,c,e in uc),F(0))
 claimed={tuple(map(F,cell)) for cell in row['witness_cell_bounds']}
 assert uc==claimed and len(good)==row['valid_rectangle_count']
 assert area==F(row['witness_union_area_exact']) and min(margins)==F(row['least_rectangle_margin_exact'])
 em=mp.mpf(eps.numerator)/eps.denominator
 extra=mp.mpf(0)
 for y in [1,-1]:
  p=lambda s:mp.mpf(35)/32-mp.mpf(3)/8*s*s+em*y*(mp.mpf(1)/2+2*s)
  roots=mp.polyroots([-mp.mpf(3)/8,2*em*y,mp.mpf(3)/32+em*y/2])
  cuts=[mp.mpf(0)]+sorted(r for r in roots if not isinstance(r,mp.mpc) and 0<r<1)+[mp.mpf(1)]
  for a,b in zip(cuts,cuts[1:]):
   if p((a+b)/2)>1:extra+=mp.quad(lambda s:p(s)*mp.log(p(s)),[a,b])/2
  fp=mp.mpf(31)/32-em*y/2
  if fp>1:extra+=fp*mp.log(fp)/2
 budget=mp.mpf(31)/16+extra;diff=abs(budget-mp.mpf(str(row['actual_Phi2_budget'])))
 assert diff<mp.mpf('1e-28')
 out.append({'epsilon':str(eps),'exact_witness_area':str(area),'strict_margin':str(min(margins)),'rectangles':len(good),'cells':len(uc),'rational_witness_verified':True,'independent_budget':mp.nstr(budget,55),'budget_difference_from_Wolfram':mp.nstr(diff,8)})
result={'python':platform.python_version(),'mpmath':mp.__version__,'working_decimal_precision':mp.mp.dps,'rows':out,'scope':'Exact rational witness verification; budget quadrature is independent floating evidence, no interval ratio certificate'}
print(json.dumps(result,ensure_ascii=False,indent=2))
