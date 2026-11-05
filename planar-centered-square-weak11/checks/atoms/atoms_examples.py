from itertools import combinations, product
from functools import cmp_to_key
from math import isqrt
import sympy as s

Q=s.Rational

def compare(a,b):
    d=s.simplify(a-b)
    if d==0:return 0
    if bool(d>0):return 1
    if bool(d<0):return -1
    raise ValueError(('undecidable comparison',a,b))

def lo(xs): return min(xs,key=cmp_to_key(compare))
def hi(xs): return max(xs,key=cmp_to_key(compare))
def rect(points, radius):
    x=[p[0] for p in points]; y=[p[1] for p in points]
    r=(hi(x)-radius,lo(x)+radius,hi(y)-radius,lo(y)+radius)
    return r if compare(r[0],r[1])<0 and compare(r[2],r[3])<0 else None

def subset_rects(points,masses,radius=None):
    out=[]
    for k in range(1,len(points)+1):
        for S in combinations(range(len(points)),k):
            mass=sum(masses[i] for i in S)
            r=s.sqrt(mass)/2 if radius is None else radius(mass)
            b=rect([points[i] for i in S],r)
            if b is not None: out.append(b)
    return out

def inclusion_area(boxes):
    total=0
    for k in range(1,len(boxes)+1):
        for bs in combinations(boxes,k):
            x0=hi([b[0] for b in bs]);x1=lo([b[1] for b in bs])
            y0=hi([b[2] for b in bs]);y1=lo([b[3] for b in bs])
            if compare(x0,x1)<0 and compare(y0,y1)<0:
                total+=(-1)**(k+1)*(x1-x0)*(y1-y0)
    return s.simplify(total)

def coordinate_rects(points,masses):
    # Independent enumeration: four atom-coordinate boundaries define each closed box.
    xs=sorted(set(p[0] for p in points));ys=sorted(set(p[1] for p in points))
    out=set()
    for xl,xr in product(xs,repeat=2):
        if xl>xr:continue
        for yl,yr in product(ys,repeat=2):
            if yl>yr:continue
            ids=[i for i,(x,y) in enumerate(points) if xl<=x<=xr and yl<=y<=yr]
            if not ids:continue
            mass=sum(masses[i] for i in ids);radius=s.sqrt(mass)/2
            # Use the chosen box, as in the source formula.
            b=(xr-radius,xl+radius,yr-radius,yl+radius)
            if compare(b[0],b[1])<0 and compare(b[2],b[3])<0:out.add(b)
    return list(out)

def sweep_area(boxes):
    xs=sorted(set(t for b in boxes for t in b[:2]),key=cmp_to_key(compare))
    total=0
    for x0,x1 in zip(xs,xs[1:]):
        ys=sorted([(b[2],b[3]) for b in boxes if compare(b[0],x0)<=0 and compare(b[1],x1)>=0],key=cmp_to_key(lambda a,b: compare(a[0],b[0])))
        if not ys:continue
        start,end=ys[0]; length=0
        for a,b in ys[1:]:
            if compare(a,end)>0:length+=end-start;start,end=a,b
            elif compare(b,end)>0:end=b
        length+=end-start;total+=(x1-x0)*length
    return s.simplify(total)

def inner(m):
    a,b=int(s.numer(m)),int(s.denom(m)); k=isqrt((a*1000**2-1)//(4*b));r=Q(k,1000)
    assert 4*r*r<m
    return r

cases=[('three',[(Q(0),Q(0)),(Q(7,10),Q(0)),(Q(7,5),Q(0))],[Q(1,3)]*3),('prompt',[(Q(0),Q(0)),(Q(19,10),Q(0))],[Q(4),Q(4)])]
for name,points,masses in cases:
    first=inclusion_area(subset_rects(points,masses));second=sweep_area(coordinate_rects(points,masses))
    expected=Q(7,3)-Q(7,5)*s.sqrt(Q(2,3)) if name=='three' else (Q(98)-39*s.sqrt(2))/5
    print(name,'subset inclusion-exclusion =',first,'coordinate sweep =',second,'expected =',expected,flush=True)
    print('both exact equal =',s.simplify(first-second)==0,'expected exact =',s.simplify(first-expected)==0,'decimal=',s.N(first,20),flush=True)
    assert s.simplify(first-second)==0 and s.simplify(first-expected)==0
    if name=='three':
        for k in range(1,4):
            for ids in combinations(range(3),k):
                m=sum(masses[i] for i in ids);r=inner(m)
                print('inner subset',ids,'mass',m,'radius',r,'strict',4*r*r<m,'box',rect([points[i] for i in ids],r))
        lower=inclusion_area(subset_rects(points,masses,inner))
        print('rational inner area =',lower,'above 7/6 =',lower>Q(7,6),'rational sweep =',sweep_area(subset_rects(points,masses,inner)))
        assert lower==Q(3702,3125) and lower>Q(7,6)
    else:
        lower=inclusion_area(subset_rects(points,masses,lambda m: Q(1) if m==4 else Q(7,5)))
        print('prompt normalized exact =',s.simplify(first/8),'r=7/5 strict =',4*Q(7,5)**2<8,'inner normalized =',lower/8,'prompt bound verified =',lower/8==Q(213,200))
        assert lower/8==Q(213,200)
