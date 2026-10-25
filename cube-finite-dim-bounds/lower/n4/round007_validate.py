"""New three-orbit interface regressions and independent point geometry."""
from pathlib import Path
from itertools import product
import json,time,sys,hashlib,numpy as np
from fractions import Fraction as F
from scipy.stats import qmc
from round007_three_orbit import A,PRIMARY,PURE,PAIRS,boxes,capture,volume
from round007_joint_model import base_phases
P=Path(__file__).resolve().parent

def main():
    st=time.monotonic();synced=[]
    for key in PRIMARY:
        ll=np.array(key[:12]).reshape(3,4).copy();hh=np.array(key[12:]).reshape(3,4).copy();ll[1]=ll[2]=ll[0];hh[1]=hh[2]=hh[0];synced.append(tuple(map(int,np.r_[ll.ravel(),hh.ravel()])))
    sigma=.001;weights=np.array([1,.8*sigma,.2*sigma]);a=A*(1+sigma)**.25
    b,kk=boxes(a,np.zeros((3,4)),weights,synced);ratio=float(volume(b)[0]/sum(weights));assert abs(ratio-1.9017426111009554)<1e-12
    phi,phiB,_=base_phases();phases=np.array([np.zeros(4),phi,phi]);pureboxes,_=boxes(A,phases,weights,PURE|PAIRS);purevolume=float(volume(pureboxes)[0]);assert abs(purevolume-sigma)<1e-12
    phases=np.array([np.zeros(4),phi,1-phi]);keys=capture(A,phases,weights,qmc.Sobol(4,scramble=True,seed=20261025709).random_base2(7));checks=0;maxabs=0
    for key in keys:
        ll=np.array(key[:12]).reshape(3,4);hh=np.array(key[12:]).reshape(3,4);points=[];mass=0.
        for orb in range(3):
            inds=list(product(*(range(int(ll[orb,j]),int(hh[orb,j])+1) for j in range(4))))
            if inds:
                arr=np.array(inds);points.extend((arr+phases[orb])*A);mass+=len(inds)*weights[orb];maxabs=max(maxabs,int(abs(arr).max()))
        points=np.array(points);r=mass**.25/2;lo=np.maximum(0,points.max(axis=0)-r);hi=np.minimum(A,points.min(axis=0)+r);bb,_=boxes(A,phases,weights,[key])
        if np.all(hi>lo+1e-14):assert len(bb)==1 and np.allclose(bb[0],np.stack([lo,hi],axis=1).ravel(),atol=1e-13,rtol=1e-13)
        else:assert not len(bb)
        checks+=1
    variation=P.parent/'variation_geometry'
    if not (variation/'vg7_joint.py').exists():variation=Path('/Users/zhengzhihao/Desktop/T/Math-Research-finite-variation/cube-finite-dim-bounds/lower/variation_geometry')
    sys.path.insert(0,str(variation))
    from vg7_joint import descriptor
    audited=PRIMARY|PURE|PAIRS|keys
    for key in audited:
        info=descriptor(key,[F(4,5),F(1,5)]);assert info['complete_cell_mass_slope']==1
        n,k1,k2=info['counts'];assert info['captured_mass_slope']==F(4*k1+k2,5)
        if n:assert info['radius_velocity_factor']*4*n==info['captured_mass_slope']
    descriptor_hash=hashlib.sha256((variation/'vg7_joint.py').read_bytes()).hexdigest()
    out=dict(exact_descriptor_key_checks=len(audited),descriptor_sha256=descriptor_hash,complete_cell_mass_slope="1",coincident_three_orbits_scaled_period_ratio=ratio,coincident_satellite_only_volume=purevolume,expected_satellite_only_volume=sigma,independent_explicit_point_keys=checks,max_abs_test_index=maxabs,pure_singleton_keys=len(PURE),pure_pair_keys=len(PAIRS),all_passed=True,elapsed=time.monotonic()-st)
    (P/'round007_validation_results.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
if __name__=='__main__':main()
