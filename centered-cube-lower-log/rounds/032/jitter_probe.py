from pathlib import Path
import numpy as np,json
rng=np.random.default_rng(2026102932);rows=[];Q=131072;x=(np.arange(Q)+.5)/Q
for K in [15,63]:
    v=4;delta=1/(4*K*v);step=1/(K*v)
    for mode in ['uniform_jitter','paired_boundary_clusters']:
        jitter=rng.random(K)
        if mode=='paired_boundary_clusters':
            jitter[::2]=1-1e-7;jitter[1::2]=1e-7
        atoms=(np.arange(K)+jitter)/K
        def primitive(t):
            periods=np.floor(t);u=t-periods;z=np.zeros_like(t)
            for a in atoms:
                for k in [-1,0,1]:
                    b=a+k
                    z+=np.clip((u-b+delta/2)/delta,0,1)-np.clip((-b+delta/2)/delta,0,1)
            return z+K*periods-K*t
        def H(L):return primitive(x+L/2)-primitive(x-L/2)
        overlap=sum((abs((x-a+.5)%1-.5)<=delta/2).astype(int) for a in atoms)
        for L in [10/9,1.16,4/3]:
            h=H(L);diff=H(L+step)-h;e=float(np.mean(diff*diff))
            assert abs(h).max()<=2+1e-8
            assert e<=4*K*step+1e-5
            assert overlap.max()<=3
            rows.append(dict(K=K,mode=mode,L=L,delta=delta,step=step,
                   H_abs_max=float(abs(h).max()),increment_energy=e,
                   short_increment_upper=4*K*step,
                   maximum_box_overlap=int(overlap.max())))
(Path(__file__).resolve().parent/'jitter_checks.json').write_text(json.dumps(dict(seed=2026102932,
          numpy_version=np.__version__,query_grid=Q,rows=rows),indent=2))
print(json.dumps(rows,indent=2))
