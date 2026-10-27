"""Minimal R037 algebra/strict-box checks, not a production rerun or LP solve."""
from fractions import Fraction as F
from pathlib import Path
import hashlib
import json
import sys

HERE=Path(__file__).resolve().parent
y=[(4,0),(-4,0),(0,4),(0,-4)]
def box(c,r):
    return [(F(t-r),F(t+r)) for t in c]
def intersects(v,w):
    return all(max(a,c)<min(b,d) for (a,b),(c,d) in zip(v,w))
def member(c,v):
    return all(a<x<b for x,(a,b) in zip(c,v))
V=[box(c,4) for c in y]
edges=[(i,j) for i in range(4) for j in range(i+1,4) if intersects(V[i],V[j])]
assert edges==[(0,2),(0,3),(1,2),(1,3)]
assert all(not intersects(box(y[i],1),box(y[j],1)) for i in range(4) for j in range(i+1,4))
# Eligibility loads normalized by p; this is the single exact prescribed order.
selected=[]
loads=[]
for i in (0,2,3,1):
    load=sum(member(y[i],box(y[j],7)) for j in selected)
    loads.append(load)
    if load<=F(3,2):
        selected.append(i)
assert loads==[0,1,1,2] and selected==[0,2,3]
E={j:box(y[j],8) for j in selected}
assert {j for j in selected if member((F(-9,2),F(0)),E[j])}=={2,3}
assert {j for j in selected if member((F(-7,2),F(0)),E[j])}=={0,2,3}
# The strips have identical width1 and common other coordinates, so each is v/2.
c=[F(1,6),F(0),F(5,12),F(5,12)]
assert sum(c,F(0))==1
assert F(1,3)<=F(2,3) and F(1,2)<=F(2,3)
q=[F(7,6),F(0),F(17,12),F(17,12)]
new=[F(7,12),F(7,12),F(17,12),F(17,12)]
assert sum(q,F(0))==sum(new,F(0))==4
assert max(q[i]+q[j] for i,j in edges)==F(31,12)
assert {new[i]+new[j] for i,j in edges}=={F(2)}
# Normalized mass primal H/2 on each edge and dual (0,0,1,1).
m=[F(1,2)]*4
z=[0,0,1,1]
assert all(sum(m[k] for k,e in enumerate(edges) if i in e)==1 for i in range(4))
assert all(z[i]+z[j]<=1 for i,j in edges)
assert sum(m)==sum(z)==2
assert (F(31,12)-2)/F(31,12)==F(7,31)>F(1,10000)

# Three exact representative pair cases: crossing interior, clipped left/right.
# These only check algebra and clipping; the general proof is monotonic-affine.
pair_checks=[]
for hi,hj,P,Ai,Aj,Ar in [(2,3,10,1,2,0),(2,3,10,20,0,0),(2,3,10,0,20,50)]:
    hi,hj,P,Ai,Aj,Ar=map(F,(hi,hj,P,Ai,Aj,Ar))
    crossing=(P+hj*(Aj-Ai))/(hi+hj)
    x=max(F(0),min(P/hi,crossing))
    def phi(v):
        return max(Ar,Ai+v,Aj+(P-hi*v)/hj)
    candidates=[F(0),P/hi]
    if 0<=crossing<=P/hi:
        candidates.append(crossing)
        assert Ai+crossing==Aj+(P-hi*crossing)/hj
    assert phi(x)==min(map(phi,candidates))
    assert x>=0 and (P-hi*x)/hj>=0
    pair_checks.append({'crossing':str(crossing),'xstar':str(x),'minimum':str(phi(x))})

receipt={'status':'passed minimal algebra checks; no production cases or LP rerun',
         'source_sha256':hashlib.sha256(Path(sys.argv[1]).read_bytes()).hexdigest(),
         'normalized_four_cube':{'strict_edges':[list(e) for e in edges],
             'eligibility_loads_over_p':loads,'selected':selected,
             'q_over_p':[str(x) for x in q],'split_over_p':[str(x) for x in new],
             'mass_primal_and_dual_over_H':'2','gap_fraction':'7/31'},
         'pair_clip_checks':pair_checks,
         'verifier_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
         'interpretation':'General-n proof uses these two active coordinates and passive equal intervals. This is not an operator bound or registered batch statistic.'}
(HERE/'verification.json').write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps(receipt,indent=2))
