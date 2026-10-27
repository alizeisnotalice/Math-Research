"""R039 minimal new rational checks; no production pair scan or LP solve.

Usage: python verify.py PROOF.md [round038/compute/cache-directory]
Optional caches are read only for source hashes and current peak intersections.
"""
from fractions import Fraction as F
from pathlib import Path
from itertools import combinations
import hashlib,json,sys
sys.set_int_max_str_digits(100000)

HERE=Path(__file__).resolve().parent
def merit(q,h,cliques):
    return (max(sum(q[k] for k in c) for c in cliques),sum(q),sum(a*b*b for a,b in zip(h,q)))

def pair(q,h,cliques,i,j):
    P=h[i]*q[i]+h[j]*q[j]
    prices={}
    for c in cliques:
        key=(i in c,j in c)
        outside=sum((q[k] for k in c if k not in (i,j)),F(0))
        prices[key]=max(prices.get(key,outside),outside)
    branches=[]
    for (ii,jj),A in prices.items():
        branches.append((F(ii)-F(jj)*h[i]/h[j],A+F(jj)*P/h[j]))
    cap=P/h[i]
    candidates={F(0),cap}
    for (a,b),(c,d) in combinations(branches,2):
        if a!=c:
            x=(d-b)/(a-c)
            if 0<=x<=cap:
                candidates.add(x)
    assert len(candidates)<=8
    phi=lambda x:max(a*x+b for a,b in branches)
    m=min(map(phi,candidates))
    lo,hi=F(0),cap
    for a,b in branches:
        if a>0: hi=min(hi,(m-b)/a)
        elif a<0: lo=max(lo,(m-b)/a)
        else: assert b<=m
    assert lo<=hi
    x=hi if h[i]>h[j] else lo if h[i]<h[j] else max(lo,min(hi,P/(h[i]+h[j])))
    new=q.copy();new[i]=x;new[j]=(P-h[i]*x)/h[j]
    assert all(v>=0 for v in new)
    assert sum(a*b for a,b in zip(h,new))==sum(a*b for a,b in zip(h,q))
    assert merit(new,h,cliques)[0]==m==phi(x)
    assert merit(new,h,cliques)<=merit(q,h,cliques)
    return new,{'pair':[i,j],'classes':len(branches),'candidate_count':len(candidates),
                'interval':[str(lo),str(hi)],'new_merit':[str(v) for v in merit(new,h,cliques)]}

# New P3 example only, normalized by h=p=1 for coefficient arithmetic.
cliques=[(0,1),(0,2)]
centers=[0,-5,5]
edges=[(i,j) for i,j in combinations(range(3),2) if max(centers[i]-4,centers[j]-4)<min(centers[i]+4,centers[j]+4)]
assert edges==cliques
assert all(abs(centers[i]-centers[j])>2 for i,j in combinations(range(3),2))
h=[F(1)]*3;q=[F(1)]*3
p3=[]
for i,j in combinations(range(3),2):
    new,row=pair(q,h,cliques,i,j)
    assert new==q and merit(new,h,cliques)==(2,3,3)
    p3.append(row)
mid=[F(1,2),F(3,2),F(1)]
end=[F(0),F(3,2),F(3,2)]
assert merit(mid,h,cliques)==(2,3,F(7,2))
assert merit(end,h,cliques)==(F(3,2),3,F(9,2))
# Explicit normalized mass primal and dual; no solver.
mass=[F(1),F(1)];dual=[0,1,1]
assert sum(mass)==sum(dual)==2
assert all(sum(mass[k] for k,c in enumerate(cliques) if i in c)>=1 for i in range(3))
assert all(sum(dual[i] for i in c)<=1 for c in cliques)
assert (F(2)-F(3,2))/2==F(1,4)

# Unequal-budget plateau: S chooses the large-budget endpoint despite E rising.
extra=[]
for hs,expected in [([2,1,1],[F(3,2),F(0),F(100)]),([1,2,1],[F(0),F(3,2),F(100)])]:
    hh=list(map(F,hs));qq=list(map(F,[1,1,100]));cs=[(0,1),(2,)]
    new,row=pair(qq,hh,cs,0,1)
    assert new==expected
    before,after=merit(qq,hh,cs),merit(new,hh,cs)
    assert before[0]==after[0] and before[1]>after[1] and before[2]<after[2]
    extra.append(row)
# All four classes present, and an adjacent pair with only an 11 branch.
new,row=pair(list(map(F,[1,2,3,4,5,6])),list(map(F,[2,3,1,1,1,1])),[(0,2),(1,3),(0,1,4),(5,)],0,1)
assert row['classes']==4 and new[0]==0 and new[1]==F(8,3) and row['new_merit'][0]=='23/3'
extra.append(row)
new,row=pair([F(1),F(1)],[F(2),F(1)],[(0,1)],0,1)
assert row['classes']==1 and new==[F(3,2),F(0)]
extra.append(row)
new,row=pair([F(0),F(0)],[F(1),F(2)],[(0,1)],0,1)
assert new==[0,0] and row['candidate_count']==1
extra.append(row)

observations=[]
if len(sys.argv)>2:
    for name,expected in [('phasea_C_2_05_full_disjoint.json',('aa259d09c5cee0e974f8d2a8f52b3bfd8b0fe1355e44007551749dc46f3726de',87,3,True)),('phasea_C_1_06_full_disjoint.json',('2c755add961ccddc1e71d482c0d0491844a6358c272094bf53b9ab5c7f1a3649',71,24,False)),('phasea_A_3_15_partial.json',('5582abceeb789961dd24714b10c336304b1a34f81f7817b8a612ce1b42c0b814',37,24,False))]:
        p=Path(sys.argv[2])/name
        digest=hashlib.sha256(p.read_bytes()).hexdigest();data=json.loads(p.read_text())
        assert digest==expected[0]
        qs=list(map(F,data['q_after']));cs=data['maximal_cliques']
        peak=max(sum(qs[i] for i in c) for c in cs)
        active=[set(c) for c in cs if sum(qs[i] for i in c)==peak]
        common=set.intersection(*active)
        assert (len(cs),len(active),len(set(qs))==1)==expected[1:]
        assert not common
        observations.append({'source':name,'sha256':digest,'complete_cached_cliques':len(cs),'active_peaks':len(active),'peak_intersection':[],'q_constant':len(set(qs))==1})

receipt={'status':'passed minimal new rational checks',
         'proof_sha256':hashlib.sha256(Path(sys.argv[1]).read_bytes()).hexdigest(),
         'P3_pair_checks':p3,'P3_mass_primal_dual':'2h','P3_gap_fraction':'1/4',
         'P3_plateau_merits':[['2','3','3'],['2','3','7/2'],['3/2','3','9/2']],
         'new_affine_cases':extra,'cached_peak_observations':observations,
         'verifier_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
         'scope':'No production pair scan, old LP/union/selection recomputation, n7, or Git mutation. Cached clique completeness is inherited, not re-enumerated.'}
(HERE/'verification.json').write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps(receipt,indent=2))
