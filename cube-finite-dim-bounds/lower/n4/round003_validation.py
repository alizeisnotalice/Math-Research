"""Finite independent explicit-atom checks of the phase/stair formulas."""
from pathlib import Path
from itertools import product
import numpy as np
import json,time
from round003_stair_search import A,templates,evaluate
P=Path(__file__).resolve().parent

def main():
    rng=np.random.default_rng(2026102504);checked=nonempty=attempts=0
    while checked<64:
        i=attempts;attempts+=1
        assert attempts<4096
        a=A*rng.uniform(.9,1.12,4);f=rng.uniform(0,.5);h=rng.uniform(.3,.7);w=rng.uniform(-.6,.6);z=np.r_[a,f,h,w]
        parity=int(rng.integers(2))
        if i%2:
            m,n=rng.integers(1,6,2);d=int(rng.integers(-2,3));l3,l4=rng.integers(1,5,2)
            ids=[(parity,k2,k3,k4) for k2,k3,k4 in product(range(m),range(l3),range(l4))]+[(parity+1,k2,k3,k4) for k2,k3,k4 in product(range(d,d+n),range(l3),range(l4))]
        else:
            ell=rng.integers(1,5,4);ids=list(product(range(parity,parity+ell[0]),range(ell[1]),range(ell[2]),range(ell[3])))
        ids=np.array(ids);sgn=1-2*(ids[:,0]%2);points=ids*a
        points[:,0]=2*a[0]*(ids[:,0]//2)+(ids[:,0]%2)*2*a[0]*h;points[:,1]+=sgn*f*a[1]/2
        mass=np.sum(1+sgn*w);r=mass**.25/2;lo=np.max(points,axis=0)-r;hi=np.min(points,axis=0)+r
        if np.all(hi>lo):
            box=np.stack([lo,hi],axis=1).ravel();bb=templates(z)
            assert np.any(np.all(np.isclose(bb,box,rtol=1e-12,atol=1e-12),axis=1));nonempty+=1
            checked+=1
    rows=[evaluate(np.r_[A,f,.5,0]) for f in [.03,.15,.25]]
    # The old 30-atom example is stair parity1,m=3,n=2,d=0,l3=2,l4=3;
    # translate its k4 interval by -1 without changing its orbit.
    result=dict(random_template_attempts=attempts,explicit_atom_template_checks=checked,nonempty_checks=nonempty,all_passed=True,old_30_atom_stair=dict(parity=1,m=3,n=2,d=0,l3=2,l4=3),diagnostics=rows,finite_index_margin=dict(b=16,rectangular_index_range=[[-8,12],[-4,7],[-4,7],[-4,7]],stair_index_range=[[-8,10],[-6,10],[-4,7],[-4,7]],label_map=['floor(k1/2)','k2','k3','k4'],mass='2*(2N+33)^4',weak_ratio_factor='((2N+1)/(2N+33))^4'))
    (P/'round003_validation_results.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
if __name__=='__main__':main()
