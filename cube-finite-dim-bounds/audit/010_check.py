#!/usr/bin/env python3
"""Independent actual-volume checks for the N0=0 4D satellite term."""
from fractions import Fraction as F
from itertools import product
from pathlib import Path
import hashlib,json

def volume(b):
    z=F(1)
    for lo,hi in b:z*=max(F(0),hi-lo)
    return z

def uncovered(cube,boxes):
    # Exact cut-cell midpoint indicators; no orthant-cone production code used.
    cuts=[]
    for j,(lo,hi) in enumerate(cube):
        cuts.append(sorted({lo,hi}|{min(hi,max(lo,b[j][side])) for b in boxes for side in (0,1)}))
    out=F(0)
    for ids in product(*(range(len(z)-1) for z in cuts)):
        b=[(cuts[j][ids[j]],cuts[j][ids[j]+1]) for j in range(4)]
        p=[(a+c)/2 for a,c in b]
        if not any(all(old[j][0]<p[j]<old[j][1] for j in range(4)) for old in boxes):out+=volume(b)
    return out

rows=[]
for t in [F(1,8),F(1,16),F(1,32),F(1,64)]:
    eps=(2*t)**4;center=[eps/3,-eps/4,eps/5,-eps/6]
    cube=[(p-t,p+t) for p in center];old=[[(eps,F(1))]*4]
    net=uncovered(cube,old);error=net-eps*F(15,16)
    bound=F(233,60)*eps*t**3
    assert F(0)<=error<=bound
    origin_cube=[(-t,t)]*4;thin=[(-eps,eps)]+[(-F(1),F(1))]*3
    overlap=eps-uncovered(origin_cube,[thin]);assert overlap==16*eps*t**3
    rows.append({'r':str(t),'epsilon':str(eps),'moving_orthant_uncovered':str(net),'leading_term':str(eps*F(15,16)),'error':str(error),'strip_error_bound':str(bound),'thin_box_overlap':str(overlap)})
t=F(1,16);eps=(2*t)**4;cube=[(-t,t)]*4
clipped=[[(F(0),F(1))]*4]
images=[[(F(k),F(k+1)) for k in ks] for ks in product((-1,0),repeat=4)]
wrong=uncovered(cube,clipped);right=uncovered(cube,images)
assert wrong==eps*F(15,16) and right==0
scale=[F(2),F(3),F(5),F(7)];J=F(210)
normalized=[(-t/a,t/a) for a in scale]
norm_uncovered=uncovered(normalized,clipped);assert J*norm_uncovered==eps*F(15,16)
out={'status':'independent exact cut-cell volume checks PASS','arithmetic':'stdlib Fraction; no author imports','moving_rows':rows,'periodic_boundary_example':{'clipped_cone_uncovered':str(wrong),'full_16_image_uncovered':str(right),'epsilon':str(eps)},'normalization_example':{'jacobian':str(J),'normalized_uncovered':str(norm_uncovered),'physical_uncovered':str(J*norm_uncovered)},'checker_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
Path(__file__).with_name('010_results.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
