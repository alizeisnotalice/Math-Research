#!/usr/bin/env python3
"""Five fixed independent geometry reconstructions; no optimizer or trajectory.

The independent implementation rebuilds atoms, all coordinate-rectangular
closures, and torus projections without calling shared geometry generators.
The production path is called only for comparison at the same five points.
Uses floating-point checks, not an interval/rational lower-bound certificate.
"""
from pathlib import Path
import importlib.util
import sys
import math
import json
from hashlib import sha256
import numpy as np

sys.dont_write_bytecode = True
ROOT = Path(__file__).resolve().parent

def independent(config):
    # Directly construct atoms and masses. No shared atom/witness/pruning/wrap code.
    hr,vr,wr,x,y,m=config
    u=(2+math.sqrt(22))/3
    h=(1+u)/2*hr;v=((1+u)/2+1)*vr;w=(u*u-1)*wr
    px,py=2*h,v
    seeds=[(0.,0.,1.),(h,0.,w)]
    if m>0:seeds.append((px*x,py*y,m))
    atoms=[(sx+i*px,sy+j*py,z) for i in range(-1,2) for j in range(-1,2) for sx,sy,z in seeds]
    xs=sorted({a[0] for a in atoms});ys=sorted({a[1] for a in atoms})
    # Bin coincident atoms at coordinate ranks; keep all their mass.
    bins={}
    for a,b,z in atoms:
        k=(xs.index(a),ys.index(b));bins[k]=bins.get(k,0.)+z
    rectangles=[];nonempty=0;positive=0
    def project(a,b,p):
        out=[]
        for k in range(math.floor(a/p)-1,math.ceil(b/p)+2):
            lo=max(0.,a-k*p);hi=min(p,b-k*p)
            if lo<hi:out.append((lo,hi))
        return out
    for l in range(len(xs)):
      for r in range(l,len(xs)):
       for b in range(len(ys)):
        for t in range(b,len(ys)):
         selected=[(i,j,z) for (i,j),z in bins.items() if l<=i<=r and b<=j<=t]
         if not selected:continue
         nonempty+=1;M=sum(z for i,j,z in selected);rad=math.sqrt(M)/2
         ax=xs[max(i for i,j,z in selected)]-rad;bx=xs[min(i for i,j,z in selected)]+rad
         ay=ys[max(j for i,j,z in selected)]-rad;by=ys[min(j for i,j,z in selected)]+rad
         if ax>=bx or ay>=by:continue
         positive+=1
         rectangles += [(a,b,c,d) for a,b in project(ax,bx,px) for c,d in project(ay,by,py)]
    return rectangles,1+w+m,{'atom_count':len(atoms),'distinct_physical_sites':len(bins),'nonempty_endpoint_rectangles':nonempty,'positive_center_rectangles':positive,'raw_wrapped_rectangles':len(rectangles),'periods':[px,py]}

def compare_masks(a,b):
    xs=sorted({z for box in a+b for z in box[:2]});ys=sorted({z for box in a+b for z in box[2:]})
    xi={v:i for i,v in enumerate(xs)};yi={v:i for i,v in enumerate(ys)}
    def cover(boxes):
        out=np.zeros((len(xs)-1,len(ys)-1),dtype=bool)
        for l,r,bot,top in boxes:out[xi[l]:xi[r],yi[bot]:yi[top]]=True
        return out
    ca,cb=cover(a),cover(b);areas=np.diff(xs)[:,None]*np.diff(ys)[None,:]
    return float(areas[ca].sum()),float(areas[cb].sum()),float(areas[ca^cb].sum())

def main():
    spec = importlib.util.spec_from_file_location('round015_geometry_audit_target', ROOT/'search.py')
    target = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = target
    spec.loader.exec_module(target)
    core, source = target.load_geometry()
    configs = [
        ('zero_baseline', (1, 1, 1, 0, 0, 0)),
        ('saved_nonzero_control', (1, 1, 1, .5, 0, 1e-4)),
        ('all_lower_geometry', (.85, .85, .7, .37, .61, .2)),
        ('all_upper_geometry', (1.15, 1.15, 1.3, 1, 1, 1)),
        ('mixed_geometry', (1.07, .91, 1.12, .5, 0, .03)),
    ]
    rows = []
    for name, config in configs:
        hr, vr, wr, x, y, mass_input = config
        p = target.parameters(core, (hr, vr, wr), x, y, mass_input)
        value, boxes, mass, area = target.geometry(core, p)
        independent_boxes, independent_mass, meta = independent(config)
        ia, pa, symmetric_difference = compare_masks(independent_boxes, boxes)
        if not (abs(mass-independent_mass) < 1e-13 and abs(ia-area) < 1e-11
                and abs(pa-area) < 1e-11 and symmetric_difference < 1e-11):
            raise ValueError('fixed geometry comparison failed: '+name)
        rows.append({'label': name, 'config': list(config), **meta,
                     'production_pruned_boxes': len(boxes), 'period_mass': mass,
                     'independent_period_mass': independent_mass,
                     'production_sweep_area': area, 'independent_cell_area': ia,
                     'same_grid_production_area': pa,
                     'symmetric_difference_area': symmetric_difference,
                     'absolute_area_difference': abs(ia-area),
                     'independent_score': ia/mass, 'production_score': value})
    saved = json.loads((ROOT/'result.json').read_text())
    if sha256((ROOT/'search.py').read_bytes()).hexdigest() != saved['search_sha256']:
        raise ValueError('search source hash differs from recorded run')
    if sha256(source.read_bytes()).hexdigest() != saved['geometry_sha256']:
        raise ValueError('round009 geometry source hash differs from recorded run')
    result = {'status': 'five_fixed_independent_geometry_checks_passed',
              'independent_fixed_configurations': 5, 'optimizer_or_search_calls': 0,
              'production_fixed_rebuilds': 5, 'independent_fixed_rebuilds': 5,
              'independence': 'independent atoms, coordinate-rank closures, integer-translate torus projection; common endpoint-cell coverage comparison',
              'checks': rows, 'tolerances': {'complete_mass': 1e-13, 'area_and_symmetric_difference': 1e-11},
              'search_sha256': saved['search_sha256'], 'geometry_sha256': saved['geometry_sha256'],
              'verifier_sha256': sha256(Path(__file__).read_bytes()).hexdigest(),
              'scope': 'five fixed floating-point checks only; no optimizer rerun, no all-trajectory geometry claim, no rigorous lower-bound certificate'}
    (ROOT/'geometry_verification.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
