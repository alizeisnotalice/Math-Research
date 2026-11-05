"""Independent rational sanity checks for the planar atomic foundation.

These are finite exact checks, not proofs of discretization, compactness, or QE.
"""
from fractions import Fraction as F
from itertools import combinations
from pathlib import Path
import hashlib
import json

def area(rect):
    l,r,b,t=rect
    return max(F(0),r-l)*max(F(0),t-b)

def intersect(rects):
    return (max(r[0] for r in rects),min(r[1] for r in rects),
            max(r[2] for r in rects),min(r[3] for r in rects))

def ie_area(rects):
    ans=F(0)
    for k in range(1,len(rects)+1):
        for group in combinations(rects,k):
            ans+=(-1)**(k+1)*area(intersect(group))
    return ans

def union_area(rects):
    rects=[r for r in rects if area(r)>0]
    if not rects:return F(0)
    xs=sorted(set(x for r in rects for x in r[:2]))
    total=F(0)
    for xl,xr in zip(xs,xs[1:]):
        xm=(xl+xr)/2
        intervals=sorted((r[2],r[3]) for r in rects if r[0]<xm<r[1])
        height=F(0)
        if intervals:
            lo,hi=intervals[0]
            for b,t in intervals[1:]:
                if b<=hi:hi=max(hi,t)
                else:height+=hi-lo;lo,hi=b,t
            height+=hi-lo
        total+=(xr-xl)*height
    assert total==ie_area(rects)
    return total

def singleton(x,y,s):return (x-s,x+s,y-s,y+s)
def pair_receiver(x,y,s):return (x-s,s,y-s,s)

# The prompt's rational certificate: singleton radii 1, pair witness r=7/5.
single=[singleton(F(0),F(0),F(1)),singleton(F(19,10),F(0),F(1))]
pair=pair_receiver(F(19,10),F(0),F(7,5))
single_area=union_area(single)
combined=union_area(single+[pair])
assert single_area==F(39,5)
assert combined-single_area==F(18,25)
assert combined/F(8)==F(213,200)
assert F(8)-4*F(7,5)**2==F(4,25)>0

# Positive masses 9/25 and 16/25: exact single radii 3/10 and 2/5;
# total mass one gives pair radius 1/2, including collisions.
collision=[]
for eps in (F(0),F(1,1000),F(1,100),F(1,10)):
    rects=[singleton(F(0),F(0),F(3,10)),
           singleton(eps,F(0),F(2,5)),pair_receiver(eps,F(0),F(1,2))]
    value=union_area(rects)
    assert value==1-eps
    collision.append({"separation":str(eps),"area":str(value)})

# Rational parametrization of square roots of masses, with sum masses one.
zero=[]
for u in (F(0),F(1,64),F(1,32),F(1,16),F(1,8)):
    root1=(1-u*u)/(1+u*u)
    root2=2*u/(1+u*u)
    m1,m2=root1*root1,root2*root2
    assert m1+m2==1
    rects=[singleton(F(0),F(0),root1/2),
           singleton(F(3,4),F(1,4),root2/2),
           pair_receiver(F(3,4),F(1,4),F(1,2))]
    value=union_area(rects)
    if u==0:assert value==1
    zero.append({"u":str(u),"mass2":str(m2),"area":str(value)})

# Separate components may be moved closer while keeping intercomponent gap >1.
component_areas=[]
for separation in (F(2),F(100)):
    rects=[singleton(F(0),F(0),F(3,10)),
           singleton(separation,F(0),F(2,5)),
           pair_receiver(separation,F(0),F(1,2))]
    value=union_area(rects)
    assert value==1
    component_areas.append({"separation":str(separation),"area":str(value)})

# A subset's bounding box may contain another atom: its mass only helps.
subset_mass,actual_box_mass,r=F(2),F(3),F(7,10)
assert subset_mass-4*r*r==F(1,25)>0
assert actual_box_mass>=subset_mass

out={"prompt_obstruction":{"singleton_union_area":str(single_area),
                            "additional_pair_area":str(combined-single_area),
                            "ratio":str(combined/F(8)),"strict_mass_margin":str(F(4,25))},
     "collision_tests":collision,"zero_mass_tests":zero,
     "component_translation_tests":component_areas,
     "extra_atom_inside_subset_box":{"subset_mass":str(subset_mass),
                                      "actual_box_mass":str(actual_box_mass),
                                      "strict_subset_margin":str(subset_mass-4*r*r)},
     "all_assertions_passed":True,
     "proofs_certified":False,
     "script_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
Path(__file__).with_suffix(".json").write_text(json.dumps(out,indent=2)+"\n")
print(json.dumps(out,indent=2))
