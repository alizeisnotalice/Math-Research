#!/usr/bin/env python3
"""Independent atom-extrema audit; never imports a producer module.

Run bounded by 120 seconds. Reconstruct both the defective published boxes
and the genuine all-atom intersections, then use audit-owned exact union.
"""
from fractions import Fraction as Q
from functools import lru_cache
from importlib.util import spec_from_file_location, module_from_spec
from itertools import product
from pathlib import Path
from math import prod
import hashlib, json, time

HERE = Path(__file__).resolve().parent
CORE_PATH = HERE/'010_independent_new_lower.py'
spec = spec_from_file_location('audit010', CORE_PATH)
core = module_from_spec(spec); spec.loader.exec_module(core)
Alg, u, v, root, union = core.Alg, core.u, core.v, core.root, core.bitmask_union
# Scalars are immutable; cached rational interval bounds do not change arithmetic.
Alg.interval = lru_cache(None)(Alg.interval)
BASE = Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube5-lower-continue/cube-dim-5-bounds')
PROD = BASE/'lower_round010'
a1, s = 1+u**3/2, (1+u)/2
h = a1/8
domain = (Alg(0), a1, Alg(0), 2*s, Alg(0), s/2, Alg(0), s/2, Alg(0), s/2)

def coeff(x):
    return {f'{i//5},{i%5}':str(c) for i,c in enumerate(x.coeff) if c}

def display(x):
    return [float(t) for t in x.interval()]

def retain(raw):
    unique = {box:tag for tag,box in raw}
    return tuple((tag,box) for box,tag in unique.items() if not any(
        other != box and all(other[2*i]<=box[2*i] and box[2*i+1]<=other[2*i+1] for i in range(5))
        for other in unique))

def construct(deadline):
    """Atom lists, not min/max parity formulas, determine first two extrema."""
    raw_bad, raw_good = [], []
    first = {}
    for l1,l2,j,k in product(range(1,5),range(1,5),range(-6,7),range(-6,7)):
        atoms = tuple((column*a1+(h if row%2 else Alg(0)), row*s)
                      for row in range(k,k+l2) for column in range(j,j+l1))
        assert len(set(atoms))==l1*l2
        first[l1,l2,j,k] = tuple((min(p[i] for p in atoms),max(p[i] for p in atoms)) for i in range(2))
    radius = {n:root(n)/2 for n in {prod(ls) for ls in product(range(1,5),repeat=5)}}
    for ls in product(range(1,5),repeat=5):
        if time.monotonic()>deadline: raise TimeoutError('120 second wall cap in geometry')
        r = radius[prod(ls)]
        # Last-three coordinate starts are fixed explicitly at consecutive grid atoms.
        tail = []
        for length in ls[2:]:
            start = -(length-1)//2 if length%2 else 1-length//2
            points = tuple(q*s for q in range(start,start+length))
            lo, hi = max(Alg(0),max(points)-r), min(s/2,min(points)+r)
            if hi<=lo: break
            tail.extend((lo,hi))
        else:
            for j,k in product(range(-6,7),repeat=2):
                extrema = first[ls[0],ls[1],j,k]
                xlo,xhi = max(Alg(0),extrema[0][1]-r), min(a1,extrema[0][0]+r)
                if xhi<=xlo: continue
                ylo,yhi = max(Alg(0),extrema[1][1]-r), min(2*s,extrema[1][0]+r)
                blo,bhi = max(Alg(0),extrema[1][0]-r), min(2*s,extrema[1][1]+r)
                tag = (j,k,ls)
                if bhi>blo: raw_bad.append((tag,(xlo,xhi,blo,bhi,*tail)))
                if yhi>ylo: raw_good.append((tag,(xlo,xhi,ylo,yhi,*tail)))
    return raw_bad,raw_good

def counterexample(kept):
    for tag,box in kept:
        j,k,ls=tag; r=root(prod(ls))/2
        if ls[1]==1: continue
        correct_lo,correct_hi=(k+ls[1]-1)*s-r,k*s+r
        if box[2]<correct_lo:
            y=(box[2]+min(correct_lo,box[3]))/2
            row=k+ls[1]-1
        elif correct_hi<box[3]:
            y=(max(correct_hi,box[2])+box[3])/2
            row=k
        else: continue
        point=tuple((box[2*i]+box[2*i+1])/2 if i!=1 else y for i in range(5))
        gap = row*s-y if row*s>y else y-row*s
        assert all(box[2*i]<point[i]<box[2*i+1] for i in range(5))
        assert gap>r
        return {'retained_tag':[j,k,list(ls)],'N':prod(ls),'radius_coefficients':coeff(r),
                'point_coefficients':[coeff(p) for p in point],
                'point_display':[display(p) for p in point],
                'box_endpoint_coefficients':[coeff(p) for p in box],
                'box_endpoint_display':[display(p) for p in box],
                'missed_selected_row':row,'missed_atom_y_coefficients':coeff(row*s),
                'y_distance_minus_radius_coefficients':coeff(gap-r),
                'y_distance_minus_radius_interval':[str(t) for t in (gap-r).interval()],
                'point_is_strictly_in_retained_box':True,'selected_atom_is_outside_common_radius_cube':True}
    raise AssertionError('expected defective retained box')

def main():
    started=time.monotonic(); deadline=started+120
    bad,good=construct(deadline); kb,kg=retain(bad),retain(good)
    example=counterexample(kb)
    # This particular retained box has a still simpler complete atom count.
    # At x=(r/2,s,s/4,s/4,s/4), exactly the row-1, column-0,
    # last-three-indices-0 atom is present in the entire periodic cube.
    r=root(4)/2
    assert example['retained_tag']==[0,0,[1,4,1,1,1]]
    assert s>r and 3*s/4>r
    assert -r/2<h<3*r/2 and h-a1 < -r/2 and h+a1>3*r/2
    assert 3*s>2*r
    example['entire_periodic_measure_mass_at_this_radius']=1
    example['actual_average_at_this_radius']='1/4'
    example['selected_four_row_y_span_exceeds_cube_diameter']=True
    vb,sb,nb=union(tuple(b for _,b in kb),deadline)
    vg,sg,ng=union(tuple(b for _,b in kg),deadline)
    rb,rg=4*vb,4*vg
    producer=json.loads((PROD/'two_row_tau_1_8_results.json').read_text())
    assert coeff(rb)==producer['ratio_coefficients']
    assert len(bad)==producer['raw_boxes'] and len(kb)==producer['retained_boxes']
    cdata=json.loads((BASE/'lower_round003/a2_equals_a3_results.json').read_text())
    C=Alg(coeff=[Q(cdata['coefficients'].get(f'{i//5},{i%5}','0')) for i in range(25)])
    delta=rg-C
    out={'status':'REJECT_published_candidate_wrong_second_coordinate_intersection',
         'independence':'No producer imports. Atom extrema rebuilt; audit-owned Q(u,v) scalar and bitmask union reused.',
         'raw_boxes':{'published_geometry':len(bad),'correct_geometry':len(good)},
         'retained_boxes':{'published_geometry':len(kb),'correct_geometry':len(kg)},
         'published_ratio_coefficients_reproduced':True,'published_ratio_display':display(rb),
         'correct_ratio_coefficients':coeff(rg),'correct_ratio_interval':[str(t) for t in rg.interval()],
         'correct_ratio_display':display(rg),
         'correct_minus_complete_C_interval':[str(t) for t in delta.interval()],
         'correct_minus_complete_C_display':display(delta),
         'correct_strictly_above_complete_C':delta.interval()[0]>0,
         'coordinate_segments':{'published':sb,'correct':sg},'sweep_nodes':{'published':nb,'correct':ng},
         'counterexample':example,
         'elapsed_seconds':time.monotonic()-started,'wall_cap_seconds':120,
         'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
         'producer_source_sha256':hashlib.sha256((PROD/'two_row_staggered.py').read_bytes()).hexdigest(),
         'producer_result_sha256':hashlib.sha256((PROD/'two_row_tau_1_8_results.json').read_bytes()).hexdigest(),
         'complete_C_source_sha256':hashlib.sha256((BASE/'lower_round003/a2_equals_a3_results.json').read_bytes()).hexdigest()}
    (HERE/'016_staggered_independent_results.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out,indent=2))

if __name__=='__main__': main()
