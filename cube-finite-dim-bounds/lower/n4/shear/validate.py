"""Exact indexing/toy references plus independent finite visibility checks."""
from fractions import Fraction as F
from itertools import product,combinations
import numpy as np
import engine

def fraction_union2(boxes):
    out=F(0)
    for n in range(1,len(boxes)+1):
        for subset in combinations(boxes,n):
            area=F(1)
            for j in range(2):area*=max(F(0),min(b[2*j+1] for b in subset)-max(b[2*j] for b in subset))
            out+=area*(1 if n%2 else -1)
    return out

def toy_area(t,singleton):
    lo=[F(-1,2)]*2 if singleton else [max(F(0),t),F(0)]
    hi=[F(1,2)]*2 if singleton else [1+min(F(0),t),F(1)]
    boxes=[]
    for k1,k2 in product(range(-3,4),repeat=2):
        shifts=[k1+t*k2,F(k2)];b=[]
        for j in range(2):b.extend([max(F(0),lo[j]+shifts[j]),min(F(1),hi[j]+shifts[j])])
        if all(b[2*j]<b[2*j+1] for j in range(2)):boxes.append(tuple(b))
    return fraction_union2(boxes)

def preflight(t):
    p,q=t.numerator,t.denominator;indices=0
    for k1,k2 in product(range(-6,7),repeat=2):
        m,j=divmod(k2,q);c=(p*j)//q;phi=F(p*j,q)-c;n1=k1+p*m+c
        assert F(n1)+phi==F(k1)+t*k2 and q*m+j==k2
        assert n1-p*m-c==k1;indices+=1
    assert len({(F(p*j,q)%1,F(j,q)) for j in range(q)})==q
    saved=engine.A;references=[]
    try:
        engine.A=np.ones(4)
        for shear in [F(-1,16),F(-1,32),F(1,32),F(1,16)]:
            for singleton,key in [(True,(1,1,0,0,0)),(False,(2,2,0,0,1,1,0,1))]:
                expected=toy_area(shear,singleton);assert expected==(F(1) if singleton else 1-abs(shear))
                b,meta=engine.periodic_boxes({key},shear);got=engine.volume(b)[0];assert abs(got-float(expected))<1e-13
                references.append(dict(t=str(shear),singleton=singleton,exact_area=str(expected),numeric=got))
    finally:engine.A=saved
    larger=np.array(list(product(range(-6,7),repeat=4)),int);xyz=engine.coordinates(larger,t);small={tuple(k) for k in engine.IDS};cases=0;shapes=set()
    # Fixed test centers; none is a new input parameter.
    centers=np.array(list(product([.17,.73],repeat=4)))*engine.A
    for x in centers:
        dist=np.max(np.abs(xyz-x),axis=1)
        for r in [.5,1,2,3]:
            ids=larger[dist<=r];assert all(tuple(k) in small for k in ids)
            if len(ids):
                key=engine.canonical(ids);_,_,n=engine.raw_box(key,t);assert n==len(ids);shapes.add(key)
            cases+=1
    return dict(exact_negative_index_supercell_checks=indices,supercell_orbits=q,supercell_mass=q,determinant='q*prod(P0)',exact_toy_reference_cases=references,larger_grid_capture_cases=cases,distinct_checked_shapes=len(shapes),proven_small_grid='[-3,4]^4',preflight_passed=True)

def row_check(t,keys):
    q=t.denominator;rows=[]
    for j in [0,1,q-1]:
        b,meta=engine.periodic_boxes(keys,t,row=j);v=engine.volume(b)[0];rows.append(dict(row=j,volume=v,boxes=meta))
    assert max(v['volume'] for v in rows)-min(v['volume'] for v in rows)<2e-12
    return dict(full_lattice_closure_checked=True,row_values=rows,axisparallel_supercell_volume=q*rows[0]['volume'],supercell_mass=q,normalized_supercell_ratio=rows[0]['volume'],all_passed=True)
