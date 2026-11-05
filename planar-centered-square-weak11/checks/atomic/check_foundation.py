from fractions import Fraction as F
from itertools import combinations,product
import json

def union_area(rects):
 xs=sorted({z for a,b,c,d in rects for z in [a,b]});total=F(0)
 for a,b in zip(xs,xs[1:]):
  x=(a+b)/2;ys=sorted((c,d) for l,r,c,d in rects if l<x<r and c<d)
  if not ys:continue
  length=F(0);l,r=ys[0]
  for c,d in ys[1:]:
   if c>r:length+=r-l;l,r=c,d
   else:r=max(r,d)
  length+=r-l;total+=(b-a)*length
 return total

def intersection_area(rs):
 return max(F(0),min(r[1] for r in rs)-max(r[0] for r in rs))*max(F(0),min(r[3] for r in rs)-max(r[2] for r in rs))
def inclusion_area(rs):
 return sum((1 if k%2 else -1)*intersection_area(sel) for k in range(1,len(rs)+1) for sel in combinations(rs,k))
R1=(-F(1),F(1),-F(1),F(1));R2=(F(9,10),F(29,10),-F(1),F(1));Rp=(F(1,2),F(7,5),-F(7,5),F(7,5))
singles=union_area([R1,R2]);certificate=union_area([R1,R2,Rp]);extra=certificate-singles
assert singles==F(39,5) and extra==F(18,25) and certificate/8==F(213,200)
assert certificate==inclusion_area([R1,R2,Rp]);assert 4*F(7,5)**2<8

def subset_halo(x,points,masses):
 for mask in range(1,1<<len(points)):
  inds=[i for i in range(len(points)) if mask>>i&1];mass=sum(masses[i] for i in inds)
  dist=max(max(abs(x[k]-points[i][k]) for k in [0,1]) for i in inds)
  if mass>4*dist**2:return True
 return False

def box_halo(x,points,masses):
 xs=sorted(set(z[0] for z in points));ys=sorted(set(z[1] for z in points))
 for a,b in combinations(xs+xs,2):
  if a>b:continue
  for c,d in combinations(ys+ys,2):
   if c>d:continue
   mass=sum(m for (u,v),m in zip(points,masses) if a<=u<=b and c<=v<=d)
   dist=max(abs(x[0]-a),abs(x[0]-b),abs(x[1]-c),abs(x[1]-d))
   if mass>4*dist**2:return True
 return False
configs=[
 ([(F(0),F(0))],[F(1)]),
 ([(F(0),F(0)),(F(19,10),F(0))],[F(4),F(4)]),
 ([(F(0),F(0)),(F(3,4),F(0)),(F(2),F(2))],[F(1,2),F(1,2),F(0)]),
 ([(F(0),F(0)),(F(0),F(0)),(F(3,5),F(1,7))],[F(1,5),F(3,10),F(1,2)]),
 ([(F(0),F(0)),(F(5,4),F(1,9)),(F(1,3),F(7,8))],[F(3,8),F(1,8),F(1,2)])]
count=0
for points,masses in configs:
 for xy in product([F(k,5) for k in range(-5,16)],repeat=2):
  assert subset_halo(xy,points,masses)==box_halo(xy,points,masses),(points,masses,xy)
  count+=1
out={'scope':'Exact finite arithmetic checks, not proofs of all-input equivalences, compactness, algebraicity or the global sharp constant.','certificate':{'singles_union':str(singles),'pair_inner_radius':'7/5','four_radius_squared':str(4*F(7,5)**2),'pair_mass':'8','extra_area':str(extra),'certified_union_area':str(certificate),'certified_ratio':str(certificate/8),'is_full_halo_area':False},'independent_area_algorithms_agree':True,'subset_vs_coordinate_box_membership_checks':count,'configurations':len(configs),'zero_mass_and_collision_included':True,'all_passed':True}
print(json.dumps(out,indent=2))
