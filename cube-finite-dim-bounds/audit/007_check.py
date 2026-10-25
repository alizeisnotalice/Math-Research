#!/usr/bin/env python3
"""Independent round005 satellite certificate: reused Newton/slab engine, new full-T3 target."""
import argparse,hashlib,itertools,json,random,time
from collections import defaultdict
from fractions import Fraction as F
from functools import lru_cache
from pathlib import Path
HERE=Path(__file__).resolve().parent
DEFAULT_INPUT=HERE.parent/'lower'/'n3'/'round005_certificate_input.json'
TARGET=HERE/'006-satellite_n3_results.json'
SEED=HERE.parent/'lower'/'n3'/'round004_certificate_input.json'
EXPECTED='5ebe44bd14b85cb8483ac648a690f9226a07da35081eda848cd0ca38785f8cd7'
D=10**30;DEN=2*D

def cube_floor(v):
    if v==0:return 0
    x=1<<((v.bit_length()+2)//3)
    while True:
        y=(2*x+v//(x*x))//3
        if y>=x:break
        x=y
    while x**3>v:x-=1
    while (x+1)**3<=v:x+=1
    assert x**3<=v<(x+1)**3
    return x
@lru_cache(None)
def radius(m):
    q=cube_floor(m.numerator*D**3//m.denominator)
    assert q**3*m.denominator<=m.numerator*D**3<(q+1)**3*m.denominator
    return q,q+1

def intcoord(z):
    z*=DEN
    assert z.denominator==1
    return z.numerator

def slabs(boxes,axis):
    events=defaultdict(int)
    for i,b in enumerate(boxes):
        events[b[2*axis]]^=1<<i;events[b[2*axis+1]]^=1<<i
    cuts=sorted(events);mask=0;out=defaultdict(int)
    for k,x in enumerate(cuts):
        mask^=events[x]
        if k+1<len(cuts) and mask:out[mask]+=cuts[k+1]-x
    assert mask==0
    return out,len(cuts)

def grid_volume(boxes):
    dims=[slabs(boxes,j) for j in range(3)]
    sx,sy,sz=[d[0] for d in dims]
    @lru_cache(None)
    def zcovered(mask):return sum(length for mz,length in sz.items() if mz&mask)
    total=0;xy=0
    for mx,lx in sx.items():
        for my,ly in sy.items():
            mask=mx&my
            if mask:total+=lx*ly*zcovered(mask);xy+=1
    return total,{'cuts':[z[1] for z in dims],'distinct_axis_masks':[len(z[0]) for z in dims],'nonzero_xy_masks':xy,'cached_z_masks':zcovered.cache_info().currsize}

def selfcheck():
    rng=random.Random(20261025006)
    for case in range(30):
        boxes=[]
        for _ in range(1+case%5):
            b=[]
            for j in range(3):b.extend(sorted(rng.sample(range(6),2)))
            boxes.append(tuple(b))
        v,_=grid_volume(boxes)
        brute=sum(any(all(b[2*j]<=x[j]<b[2*j+1] for j in range(3)) for b in boxes) for x in itertools.product(range(5),repeat=3))
        assert v==brute
    return {'seed':20261025006,'cases':30,'pass':True}

def run(input_path,author_path=None,target_path=TARGET,seed_path=SEED):
    start=time.monotonic();assert hashlib.sha256(input_path.read_bytes()).hexdigest()==EXPECTED
    inp=json.loads(input_path.read_text());periods=list(map(F,inp['periods']));phases=[list(map(F,z)) for z in inp['phases']];weights=list(map(F,inp['weights']))
    assert len(periods)==3 and len(phases)==len(weights)==2 and all(z>0 for z in periods+weights)
    assert all(0<=z<1 for phase in phases for z in phase)
    assert any((phases[1][j]-phases[0][j]).denominator!=1 for j in range(3))
    inner=[];outer=[];records=[];indices=[];point_checks=0;critical=0
    periodints=list(map(intcoord,periods))
    for key in inp['witness_keys']:
        assert len(key)==12 and all(type(z) is int for z in key)
        atoms=[];mass=F(0);counts=[]
        for orbit in range(2):
            low=key[orbit*3:orbit*3+3];high=key[6+orbit*3:9+orbit*3]
            points=list(itertools.product(*(range(a,b+1) for a,b in zip(low,high))))
            counts.append(len(points));mass+=weights[orbit]*len(points)
            for k in points:
                indices.append(k)
                atoms.append(tuple(intcoord(periods[j]*(k[j]+phases[orbit][j])) for j in range(3)))
        assert atoms and mass>0
        rlo,rhi=radius(mass)
        critical+=rlo**3*mass.denominator==mass.numerator*D**3
        vmin=[min(z[j] for z in atoms) for j in range(3)];vmax=[max(z[j] for z in atoms) for j in range(3)]
        ib=tuple(z for j in range(3) for z in (max(0,vmax[j]-rlo),min(periodints[j],vmin[j]+rlo)))
        ob=tuple(z for j in range(3) for z in (max(0,vmax[j]-rhi),min(periodints[j],vmin[j]+rhi)))
        assert all(ib[2*j]<ib[2*j+1] for j in range(3))
        assert all(ob[2*j]<=ib[2*j] and ib[2*j+1]<=ob[2*j+1] for j in range(3))
        for atom in atoms:
            assert all(max(abs(ib[2*j]-atom[j]),abs(ib[2*j+1]-atom[j]))<=rlo for j in range(3))
            point_checks+=1
        inner.append(ib);outer.append(ob)
        records.append({'mass':str(mass),'orbit_counts':counts,'radius_lower':str(F(rlo,DEN)),'radius_upper':str(F(rhi,DEN))})
    index_min=[min(z[j] for z in indices) for j in range(3)];index_max=[max(z[j] for z in indices) for j in range(3)];b=max(abs(z) for z in index_min+index_max)
    assert b==2 and all(-b<=z<=b for k in indices for z in k)
    # This is the exact symbolic collar implication for every integer N>=0.
    # -b<=k_j<=b and -N<=q_j<=N imply -N-b<=q_j+k_j<=N+b.
    assert weights[0]+weights[1]==F('1.000561416422')
    ilow,istats=grid_volume(inner);ihigh,ostats=grid_volume(outer)
    totalmass=sum(weights);low=F(ilow,DEN**3)/totalmass;high=F(ihigh,DEN**3)/totalmass
    previous=json.loads(target_path.read_text());oldlo,oldhi=map(F,previous['ratio_interval'])
    assert previous['input_sha256']=='5f7c4e028a4b3d34d119410c24db9a569f96146427fdb60bd7eefed0e5515970'
    seed=json.loads(seed_path.read_text());assert inp['witness_keys']==seed['witness_keys'] and [[F(z) for z in phase] for phase in inp['phases']]==[[F(z) for z in phase] for phase in seed['phases']] and list(map(F,inp['weights']))==list(map(F,seed['weights']))
    assert list(map(F,inp['periods'][1:]))==list(map(F,seed['periods'][1:])) and inp['periods'][0]=='1.630361663339249'
    assert low>oldhi
    author_overlap=None
    if author_path is not None:
        author=json.loads(author_path.read_text());authorlo,authorhi=map(F,author['ratio_interval'])
        assert max(low,authorlo)<min(high,authorhi)
        author_overlap=True
    assert low-oldhi>F('0.00000008467763')
    # Compare interval widths and new exact lower to a short decimal rational independently.
    shown=F('1.7939726890816552');assert low>shown
    return {'status':'exact independent finite-witness lower certificate PASS','input_sha256':EXPECTED,'checker_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'target_certificate_sha256':hashlib.sha256(target_path.read_bytes()).hexdigest(),'target_exact_upper':str(oldhi),'seed_input_sha256':hashlib.sha256(seed_path.read_bytes()).hexdigest(),'only_P1_changed_and_170_keys_preserved':True,'root_denominator':str(D),'coordinate_denominator':str(DEN),'witness_count':len(inner),'enumerated_atom_endpoint_checks':point_checks,'critical_equal_radius_count':critical,'index_min':index_min,'index_max':index_max,'b':b,'mass_per_index_pair':str(totalmass),'ratio_interval':[str(low),str(high)],'decimal_display_only':[float(low),float(high)],'exact_gain_above_T3_full_upper':str(low-oldhi),'lower_above_T3_full_upper':low>oldhi,'author_interval_overlap':author_overlap,'inner_grid':istats,'outer_grid':ostats,'selfcheck':selfcheck(),'minimum_inner_side':str(F(min(ib[2*j+1]-ib[2*j] for ib in inner for j in range(3)),DEN)),'radius_records':records,'elapsed_seconds':time.monotonic()-start}
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--input',type=Path,default=DEFAULT_INPUT);ap.add_argument('--author-result',type=Path);ap.add_argument('--target-certificate',type=Path,default=TARGET);ap.add_argument('--seed-input',type=Path,default=SEED);ap.add_argument('--save',action='store_true');args=ap.parse_args()
    out=run(args.input,args.author_result,args.target_certificate,args.seed_input)
    if args.save:HERE.joinpath('007_results.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({k:v for k,v in out.items() if k!='radius_records'},indent=2))
