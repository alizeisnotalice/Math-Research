"""Separate explicit interval-placement check of the analytical safe bound."""
import importlib.util,itertools,json,math
from pathlib import Path
import numpy as np
HERE=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('probe',HERE/'probe.py');p=importlib.util.module_from_spec(spec);spec.loader.exec_module(p)
rng=np.random.default_rng(202610076333);records=[]
for n in [1,2]:
    for k in range(20):
        L=2;Z,omega,labels=p.draws(n,L,1,202610079000+n*100+k);r=float(rng.uniform(1,2));dx=r*omega/2
        scores,rr,lm,lr,ix,valid=p.candidates(dx,Z,L);R=float(rr[0,ix[0]])
        masses={}
        for q,prior in zip(p.Q,p.P):
            for key in itertools.product(range(L*q),repeat=n):
                zz=tuple(v*(4//q) for v in key);masses[zz]=masses.get(zz,0.)+float(prior)/(L*q)**n
        atoms=np.array(list(masses))/4;w=np.array(list(masses.values()));local=atoms-Z[0]/4
        # The explicit check works in local physical coordinates to retain tagged
        # source faces. 1e-12 only joins numerically identified boundary arrivals.
        inside=np.max(2*np.abs(local-dx[0]),axis=1)<=R+1e-12
        local=local[inside];w=w[inside];h=2/n
        possible=[np.unique(np.concatenate((local[:,i],local[:,i]-h))) for i in range(n)]
        best=0.
        for lower in itertools.product(*possible):
            low=np.array(lower);capt=np.all((local>=low-1e-12)&(local<=low+h+1e-12),axis=1)
            best=max(best,float(w[capt].sum()))
        upper=float(math.exp(lr[0,ix[0]]));ratio=best/float(w.sum())
        assert ratio<=upper+1e-11,(n,k,ratio,upper)
        records.append(dict(n=n,case=k,winner_R=R,explicit_sup_captured_ratio=ratio,safe_component_product_ratio_upper=upper))
face=[]
for n in [16,32,64]:
    for sign in [-1,1]:
        for r in [1.,math.sqrt(2),2.]:
            dx=np.zeros((1,n));dx[0,0]=sign*r/2
            # common quarter-unit coordinate as high as support permits.
            Z=np.full((1,n),252);sc,rr,lm,lr,ix,valid=p.candidates(dx,Z,64)
            assert 2*abs(dx[0,0])==r
            face.append(dict(n=n,sign=sign,r=r,tagged_source_arrival=2*abs(dx[0,0]),winner_R=float(rr[0,ix[0]])))
(HERE/'microbox_validation.json').write_text(json.dumps(dict(status='PASS_FLOAT_EXPLICIT_MICROBOX_PLACEMENTS',cases=40,records=records,local_face_tests=face,tolerance=1e-12,scope='Safe per-component product inequality is analytic; numerical comparisons and winner cross-check are not interval proof.'),indent=2))
print('MICROBOX PASS',len(records))
