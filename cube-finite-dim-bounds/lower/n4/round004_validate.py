"""Coincident-orbit scaling regression and explicit point geometry checks."""
from pathlib import Path
import json
import numpy as np
from scipy.stats import qmc
from round004_satellite_search import A,PRIMARY,capture,witness_boxes,volume
P=Path(__file__).resolve().parent

def main():
    synchronized=[]
    for k in PRIMARY:
        q=np.array(k).copy();q[4:8]=q[:4];q[12:16]=q[8:12];synchronized.append(tuple(q))
    calibration=[]
    for eps in [.001,.01,.05]:
        a=A*(1+eps)**.25;b=witness_boxes(a,np.zeros(4),eps,synchronized);r=float(volume(b)[0]/(1+eps))
        assert abs(r-1.9017426111009554)<1e-12
        calibration.append(dict(epsilon=eps,scale=(1+eps)**.25,ratio=r,box_count=len(b)))
    # Regenerate independent rectangles at a genuinely free phase, then explicit atom extrema/mass.
    phi=np.array([.17,.73,.39,.91]);eps=.01;t=qmc.Sobol(4,scramble=True,seed=20261025406).random_base2(6)
    keys=capture(A,phi,eps,t);checks=0;maxabs=0
    from itertools import product
    for key in keys:
        k=np.array(key);ll=k[:8].reshape(2,4);hh=k[8:].reshape(2,4);pts=[];mass=0
        for orb in range(2):
            indices=list(product(*(range(int(ll[orb,j]),int(hh[orb,j])+1) for j in range(4))))
            if not indices:continue
            ii=np.array(indices);maxabs=max(maxabs,int(abs(ii).max()));pts.extend((ii+(phi if orb else 0))*A);mass+=len(ii)*(eps if orb else 1)
        pts=np.array(pts);r=mass**.25/2;lo=np.maximum(0,pts.max(axis=0)-r);hi=np.minimum(A,pts.min(axis=0)+r)
        b=witness_boxes(A,phi,eps,[key]);expected=np.stack([lo,hi],axis=1).ravel()
        if np.all(hi>lo+1e-13):assert len(b)==1 and np.allclose(b[0],expected,rtol=1e-13,atol=1e-13)
        else:assert len(b)==0
        checks+=1
    out=dict(coincident_orbit_scaled_period=calibration,explicit_two_orbit_key_checks=checks,all_passed=True,test_phase=phi.tolist(),test_epsilon=eps,test_seed=20261025406,max_abs_captured_index=maxabs,transfer_margin=8,transfer_mass='(1+epsilon)*(2N+17)^4',transfer_factor='((2N+1)/(2N+17))^4')
    (P/'round004_validation.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
if __name__=='__main__':main()
