from fractions import Fraction as F
from pathlib import Path
import hashlib
import json


def intervals(mu):
    out = []
    for i in range(len(mu)):
        mass = F(0)
        for j in range(i, len(mu)):
            mass += mu[j][1]
            l, r = mu[j][0]-mass, mu[i][0]+mass
            if l <= r:
                out.append((l, r))
    return out


def union_intervals(iv):
    out = []
    for l,r in sorted(iv):
        if not out or l > out[-1][1]:
            out.append([l,r])
        else:
            out[-1][1] = max(out[-1][1],r)
    return out


def data(mu):
    mu = sorted(mu)
    merged = union_intervals(intervals(mu))
    mass = sum((m for _,m in mu),F(0))
    length = sum((r-l for l,r in merged),F(0))
    sep = all(mu[i+1][0]-mu[i][0] > mu[i+1][1]+mu[i][1]
              for i in range(len(mu)-1))
    connected = merged == [[mu[0][0]-mu[0][1],
                            mu[-1][0]+mu[-1][1]]]
    return mass, length, length/(2*mass), sep, connected


def concat(mu, count):
    base = mu[0][0]
    z = [(x-base,m) for x,m in mu]
    D = z[-1][0]
    all_atoms = {}
    for j in range(count):
        for x,m in z:
            xx = j*D+(x if j%2 == 0 else D-x)
            if xx in all_atoms:
                assert all_atoms[xx] == m
            else:
                all_atoms[xx] = m
    return sorted(all_atoms.items())


def rect_area(rects):
    xs = sorted({x for l,r,b,t in rects for x in (l,r)})
    ys = sorted({y for l,r,b,t in rects for y in (b,t)})
    total = F(0)
    for l,r in zip(xs,xs[1:]):
        for b,t in zip(ys,ys[1:]):
            x,y=(l+r)/2,(b+t)/2
            if any(ll<x<rr and bb<y<tt for ll,rr,bb,tt in rects):
                total += (r-l)*(t-b)
    return total


seed = [(F(0),F(1)),(F(3),F(1))]
M,_,R,sep,connected = data(seed)
assert sep and connected and R == F(5,4)
alpha = 2*R
A = (alpha**2-alpha)*M+(alpha-1)*seed[0][1]
nu = seed+[(alpha*x+A,alpha*m) for x,m in seed]
last, s, r = len(nu)-1,len(nu)-1,len(seed)
y0 = 2*nu[0][0]-nu[s][0]-2*nu[0][1]-nu[s][1]
k0 = nu[s][0]-nu[0][0]-sum((m for _,m in nu[1:s]),F(0))
yn = 2*nu[last][0]-nu[r][0]+2*nu[last][1]+nu[r][1]
kn = nu[last][0]-nu[r][0]-sum((m for _,m in nu[r+1:last]),F(0))
transformed = [(y0,k0)]+nu+[(yn,kn)]
MT,ET,RT,sep,connected = data(transformed)
assert sep and connected
D = yn-y0
den = 2*MT-k0-kn
star = D/den
fR = (20*R*R-4*R)/(12*R*R-2*R+1)
assert star == fR == F(35,23)
rows = []
for L in [1,2,3,4,8,16,32]:
    atoms = concat(transformed,L)
    ML,EL,RL,sep,connected = data(atoms)
    B = 2*k0 if L%2 == 0 else k0+kn
    assert len(atoms) == L*(len(transformed)-1)+1
    assert sep and connected
    assert RL == (L*D+B)/(L*den+B)
    assert star-RL == B*(D-den)/(den*(L*den+B))
    rows.append({'copies':L,'atoms':len(atoms),'mass':str(ML),
                 'receiver_length':str(EL),'ratio':str(RL),
                 'gap_to_modified_ratio':str(star-RL)})

d = F(19,10)
r = F(7,5)
singles = [(F(-1),F(1),F(-1),F(1)),
           (d-1,d+1,F(-1),F(1))]
inner_pair = (d-r,r,-r,r)
assert 8 > 4*r*r
single_union = rect_area(singles)
cert_area = rect_area(singles+[inner_pair])
assert single_union == F(39,5)
assert cert_area == F(213,25)
assert cert_area/8 == F(213,200)

def q2mul(x,y):
    return (x[0]*y[0]+2*x[1]*y[1],x[0]*y[1]+x[1]*y[0])

pair_extra = q2mul((F(-19,10),F(2)),(F(-2),F(2)))
exact_area = (F(39,5)+pair_extra[0],pair_extra[1])
exact_ratio = tuple(c/8 for c in exact_area)
assert exact_area == (F(98,5),F(-39,5))
assert exact_ratio == (F(49,20),F(-39,40))
assert exact_ratio[1] != 0  # sqrt(2) irrational, hence this ratio irrational.

def pmul(a,b):
    c=[0]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b):
            c[i+j]+=x*y
    return c

def psub(a,b):
    return [(a[i] if i<len(a) else 0)-(b[i] if i<len(b) else 0)
            for i in range(max(len(a),len(b)))]

NN,DD=[0,-4,20],[1,-2,12]
assert psub(NN,[0]+DD)==[0,-5,22,-12]
assert psub(pmul([-4,40],DD),pmul(NN,[-2,24]))==[-4,40,8,0]

out = {
    'scope':'finite exact checks; no planar sharp upper bound or fixed-N optimization',
    'seed_ratio':str(R),'alpha':str(alpha),'translation_A':str(A),
    'T_seed':[[str(x),str(m)] for x,m in transformed],
    'T_seed_plain_ratio':str(RT),'T_seed_modified_ratio':str(star),
    'reflection_concatenation':rows,
    'planar_prompt':{'single_union':str(single_union),
                     'inner_radius':str(r),'strict_pair_witness':True,
                     'certificate_area':str(cert_area),
                     'certificate_ratio':str(cert_area/8),
                     'full_area_Qsqrt2':[str(c) for c in exact_area],
                     'full_ratio_Qsqrt2':[str(c) for c in exact_ratio],
                     'rational_input_has_irrational_ratio':True},
    'common_radius':[
        {'R':str(q),'product_of_1d_maxima':str(F(1,4*q)),
         'same_radius_square_maximum':str(F(1,4*q*q)),
         'ratio':str(q)} for q in [2,3,10,100]],
    'all_assertions_passed':True,
    'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
}
Path(__file__).with_suffix('.json').write_text(
    json.dumps(out,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(out,ensure_ascii=False,indent=2))
