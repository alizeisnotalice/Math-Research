#!/usr/bin/env python3
"""Strict rational reconstruction of two closed off-lattice diagnostics."""
from fractions import Fraction as Q
from functools import lru_cache
from itertools import product
from pathlib import Path
from math import prod
from decimal import Decimal,localcontext
import json,time,hashlib,sys
HERE=Path(__file__).resolve().parent
DEN=2**64;EPS=Q(1,1024);SHIFT=(0,1,0,0,0)
def floorroot(n,k=5):
    lo,hi=0,1<<((n.bit_length()+k-1)//k)
    while lo+1<hi:
        m=(lo+hi)//2
        if m**k<=n:lo=m
        else:hi=m
    return lo
def root_interval(q,den=2**100):
    q=Q(q);v=floorroot(q.numerator*den**5//q.denominator)
    return Q(v,den),Q(v+1,den)
@lru_cache(None)
def side(mass):
    mass=Q(mass);v=floorroot(mass.numerator*DEN**5//mass.denominator)
    if Q(v,DEN)**5==mass:v-=1
    s=Q(v,DEN);assert s**5<mass and (s+Q(1,DEN))**5>=mass
    return s
def nearest(shift,length):
    # Independent explicit consecutive sets; one definition, checked per atom below.
    if shift:
        start=-(length//2)
    else:
        start=-((length-1)//2)
    return tuple(range(start,start+length))
def rational_spacings(codilate):
    ul,uh=root_interval(2);vl,vh=root_interval(Q(1025,1024)) if codilate else (Q(1),Q(1))
    bounds=[((1+ul**3/2)*vl,(1+uh**3/2)*vh)]+[((1+ul)/2*vl,(1+uh)/2*vh)]*4
    a=[]
    for low,high in bounds:
        n=(low*DEN).numerator//(low*DEN).denominator
        assert Q(n+1,DEN)>high
        if Q(n,DEN)==low:n-=1
        assert Q(n,DEN)<high
        a.append(Q(n,DEN))
    return tuple(a)
def run(codilate):
    begin=time.monotonic();nodes=0
    def budget():
        if time.monotonic()-begin>120 or nodes>200000:raise TimeoutError('registered cap')
    a=rational_spacings(codilate);raw={};configs=0
    def add(l0,l1):
        nonlocal configs
        configs+=1
        mass=(prod(l0) if l0 else 0)+EPS*(prod(l1) if l1 else 0)
        r=side(mass)/2;b=[]
        for j in range(5):
            mins=[];maxs=[]
            for orbit,ls in ((0,l0),(1,l1)):
                if ls:
                    shift=SHIFT[j] if orbit else 0
                    idx=nearest(shift,ls[j])
                    mins.append((min(idx)+Q(shift,2))*a[j]);maxs.append((max(idx)+Q(shift,2))*a[j])
            lo=max(Q(0),max(maxs)-r);hi=min(a[j]/2,min(mins)+r)
            if hi<=lo:return
            b.extend((lo,hi))
        raw[tuple(b)]=(l0,l1,mass,r)
    for ls in product(range(1,5),repeat=5):add(ls,None);add(None,ls)
    opts=[[(p,q) for p in range(1,5) for q in range(1,5) if (abs(p-q)<=1 if s else p==q)] for s in SHIFT]
    for pq in product(*opts):add(tuple(t[0] for t in pq),tuple(t[1] for t in pq))
    kept=tuple(sorted(b for b in raw if not any(o!=b and all(o[2*j]<=b[2*j] and b[2*j+1]<=o[2*j+1] for j in range(5)) for o in raw)))
    guards=0;max_indices=[0]*5
    for b in kept:
        l0,l1,mass,r=raw[b];atoms={}
        for orbit,ls,weight in ((0,l0,Q(1)),(1,l1,EPS)):
            if ls is None:continue
            for idx in product(*(nearest(SHIFT[j] if orbit else 0,l) for j,l in enumerate(ls))):
                point=tuple(2*z+(SHIFT[j] if orbit else 0) for j,z in enumerate(idx))
                assert point not in atoms
                atoms[point]=weight
                max_indices=[max(old,abs(z)) for old,z in zip(max_indices,idx)]
        assert sum(atoms.values())==mass and (2*r)**5<mass
        for j in range(5):
            co=[p[j]*a[j]/2 for p in atoms]
            assert max(co)-r<=b[2*j] and b[2*j+1]<=min(co)+r
            guards+=2
    @lru_cache(None)
    def sweep(boxes):
        nonlocal nodes
        nodes+=1;budget()
        if not boxes:return Q(0)
        if len(boxes[0])==2:
            total=Q(0);left=right=None
            for lo,hi in sorted(boxes):
                if left is None:left,right=lo,hi
                elif right<lo:total+=right-left;left,right=lo,hi
                elif hi>right:right=hi
            return total+right-left
        cuts=sorted({v for b in boxes for v in b[:2]});total=Q(0)
        for lo,hi in zip(cuts,cuts[1:]):
            sub=tuple(sorted({b[2:] for b in boxes if b[0]<=lo and hi<=b[1]}))
            total+=(hi-lo)*sweep(sub)
        return total
    volume=sweep(kept);segments=[]
    for j in range(5):
        cuts=sorted({b[2*j+k] for b in kept for k in (0,1)});axis=[]
        for lo,hi in zip(cuts,cuts[1:]):
            mask=sum(1<<q for q,b in enumerate(kept) if b[2*j]<=lo and hi<=b[2*j+1])
            if mask:axis.append((hi-lo,mask))
        segments.append(axis)
    @lru_cache(None)
    def grid(j,mask):
        if j==5:return Q(1)
        return sum((d*grid(j+1,mask&m) for d,m in segments[j] if mask&m),Q(0))
    assert grid(0,(1<<len(kept))-1)==volume
    ratio=32*volume/(1+EPS)
    cp=HERE.parent/'lower_round003/a2_equals_a3_results.json';baseline=json.loads(cp.read_text());cl,ch=map(Q,baseline['ratio_interval'])
    with localcontext() as ctx:
        ctx.prec=55;display=str(Decimal(ratio.numerator)/Decimal(ratio.denominator))
    out=dict(status='NO_IMPROVEMENT_exact' if ratio<cl else ('PASS_above_C_pending_independent_audit' if ratio>ch else 'UNRESOLVED'),codilation=codilate,dimension=5,spacing_rationals=[str(x) for x in a],epsilon=str(EPS),half_shift=list(SHIFT),period_mass=str(1+EPS),configs=configs,raw_boxes=len(raw),retained_boxes=len(kept),ratio_exact=str(ratio),ratio_decimal_display_only=display,strictly_below_C_lower=ratio<cl,strictly_above_C_upper=ratio>ch,difference_from_C_interval=[str(ratio-ch),str(ratio-cl)],union_engines_agree=True,coordinate_segments=[len(x) for x in segments],nodes=nodes,endpoint_inequalities=guards,distinct_atom_mass_guard=True,strict_common_radius_guard=True,max_period_indices=max_indices,wall_seconds=time.monotonic()-begin,baseline_source_sha256=hashlib.sha256(cp.read_bytes()).hexdigest(),source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    out['witnesses']=[dict(lengths_primary=raw[b][0],lengths_satellite=raw[b][1],mass=str(raw[b][2]),radius=str(raw[b][3]),box=[str(v) for v in b]) for b in kept]
    name='exact_codilation' if codilate else 'exact_fixed'
    (HERE/(name+'_results.json')).write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({k:v for k,v in out.items() if k not in ('witnesses','ratio_exact','difference_from_C_interval')},indent=2))
if __name__=='__main__':run(sys.argv[1]=='codilation')
