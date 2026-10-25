"""Supplement only: exact periodic tangent-cone theta for stored critical phases.

No optimization or finite-epsilon search. Q(u,v) coordinates are mapped to
rational ranks only after rigorous interval-separated ordering. Theta uses
only coordinate order/equality, so the Fraction interface receives exact data.
"""
from pathlib import Path
from functools import lru_cache
from fractions import Fraction as F
from itertools import product
import importlib.util,json,hashlib,time
import numpy as np
from round005_exposure import AV,rvec
P=Path(__file__).resolve().parent
INTERFACE=P.parent/'variation_geometry/normal_intervals.py'
if not INTERFACE.exists():INTERFACE=Path('/Users/zhengzhihao/Desktop/T/Math-Research-finite-variation/cube-finite-dim-bounds/lower/variation_geometry/normal_intervals.py')
spec=importlib.util.spec_from_file_location('variation_normal_intervals',INTERFACE);normal=importlib.util.module_from_spec(spec);spec.loader.exec_module(normal)
cert=json.loads((P.parents[1]/'baseline/original_certificates.json').read_text())[1]
D=cert['root_denominator'];nums=cert['root_lower_numerators']
for k,target in zip(nums,[2,3]):assert k**4<target*D**4<(k+1)**4
u=(F(nums[0],D),F(nums[0]+1,D));v=(F(nums[1],D),F(nums[1]+1,D))
B=[(u[0]**i*v[0]**j/4,u[1]**i*v[1]**j/4) for i,j in product(range(4),repeat=2)]
ZERO=(0,)*16;AA=[tuple(map(int,row)) for row in AV]
def add(a,b):return tuple(x+y for x,y in zip(a,b))
def scale(a,c):return tuple(c*x for x in a)
@lru_cache(None)
def interval(vec):
    lo=hi=F(0)
    for c,(a,b) in zip(vec,B):
        lo+=c*(a if c>=0 else b);hi+=c*(b if c>=0 else a)
    return lo,hi
@lru_cache(None)
def compare(a,b):
    z=add(a,scale(b,-1))
    if z==ZERO:return 0
    lo,hi=interval(z)
    if lo>0:return 1
    if hi<0:return -1
    raise AssertionError('unresolved exact algebraic order')

def baseline_boxes():
    out=set()
    for ell in product(range(1,5),repeat=4):
        rv=tuple(map(int,rvec(int(np.prod(ell)))));choices=[]
        for j in range(4):
            opts=[]
            for k in range(-2,2):
                lo=add(scale(AA[j],k+ell[j]-1),scale(rv,-1));hi=add(scale(AA[j],k),rv)
                if compare(lo,ZERO)<0:lo=ZERO
                if compare(hi,AA[j])>0:hi=AA[j]
                if compare(lo,hi)<0:opts.append((lo,hi))
            choices.append(opts)
        for box in product(*choices):out.add((tuple(x[0] for x in box),tuple(x[1] for x in box)))
    return sorted(out)

def center(formula):
    if formula['kind']=='coincident':return (ZERO,)*4
    root=scale(tuple(map(int,rvec(formula['primary_mass']))),2);cc=[]
    for j,sign in enumerate(formula['signs']):
        raw=scale(root,sign);ll,hh=interval(raw);al,ah=interval(AA[j]);ratios=[x/y for x,y in product([ll,hh],[al,ah])]
        m0=min(ratios).numerator//min(ratios).denominator;m1=max(ratios).numerator//max(ratios).denominator
        assert m0==m1
        z=add(raw,scale(AA[j],-m0));assert compare(z,ZERO)>=0 and compare(z,AA[j])<0;cc.append(z)
    return tuple(cc)

def main():
    st=time.monotonic();data=json.loads((P/'round005_exposure_results.json').read_text());boxes=baseline_boxes()
    periodic=[]
    for shift in product([-1,0,1],repeat=4):
        displacement=[scale(AA[j],shift[j]) for j in range(4)]
        for lo,hi in boxes:periodic.append((tuple(add(lo[j],displacement[j]) for j in range(4)),tuple(add(hi[j],displacement[j]) for j in range(4))))
    centers=[center(row['phase_formula']) for row in data['phases']]
    ranks=[]
    for j in range(4):
        values={c[j] for c in centers}
        for lo,hi in periodic:values.update([lo[j],hi[j]])
        values=sorted(values,key=lambda x:sum(interval(x)))
        for a,b in zip(values,values[1:]):assert interval(a)[1]<interval(b)[0]
        ranks.append({x:i for i,x in enumerate(values)})
    lows=np.array([[ranks[j][lo[j]] for j in range(4)] for lo,hi in periodic]);highs=np.array([[ranks[j][hi[j]] for j in range(4)] for lo,hi in periodic])
    out=[]
    for row,c in zip(data['phases'],centers):
        z=np.array([ranks[j][c[j]] for j in range(4)]);active=np.all((lows<=z)&(highs>=z),axis=1)
        exact_boxes=[(tuple(map(int,lo)),tuple(map(int,hi))) for lo,hi in zip(lows[active],highs[active])]
        theta=normal.pure_satellite_term(tuple(map(int,z)),1,exact_boxes)
        interior=bool(np.any(np.all((lows<z)&(highs>z),axis=1)))
        out.append(dict(phase_index=row['phase_index'],phase_formula=row['phase_formula'],theta=str(theta),active_periodic_images=int(sum(active)),strictly_inside_one_old_box=interior,old_best_Oepsilon_prediction=row['best_prediction'],corrected_best_prediction=row['best_prediction']+float(theta)))
    runs=[]
    for row in data['runs']:
        idx=next(i for i,x in enumerate(data['phases']) if x['phase_formula']==row['phase_formula']);theta=F(out[idx]['theta'])
        runs.append(dict(phase_formula=row['phase_formula'],direction=row['direction'],theta=str(theta),old_Oepsilon_prediction=row['prediction'],corrected_prediction=row['prediction']+float(theta)))
    result=dict(status='supplemental exact theta only; no optimization or finite-epsilon rerun',old_model_scope='N0>0 fixed-radius witnesses with O(epsilon) endpoints; pure N0=0 singleton not explicitly retained',physical_term='mass coefficient1 times theta, added once to ratio derivative; no determinant factor',periodic_image_shifts=[-1,0,1],baseline_unique_boxes=len(boxes),periodic_boxes=len(periodic),algebraic_order_certified=True,rank_counts=list(map(len,ranks)),phase_count=len(out),theta_histogram={q:sum(r['theta']==q for r in out) for q in sorted(set(r['theta'] for r in out))},phases=out,optimized_runs=runs,positive_stored_corrected_prediction=any(r['corrected_best_prediction']>1e-9 for r in out) or any(r['corrected_prediction']>1e-9 for r in runs),elapsed=time.monotonic()-st,input_sha256={str(q):hashlib.sha256(q.read_bytes()).hexdigest() for q in [P/'round005_exposure_results.json',P.parents[1]/'baseline/original_certificates.json',INTERFACE]})
    (P/'round005_pure_satellite_results.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:v for k,v in result.items() if k!='phases'},indent=2))
if __name__=='__main__':main()
