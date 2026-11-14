"""Exact-rational primal/dual certificates for finite rectangle catalogues.
These are fixed-catalogue physical indicator certificates, not general H_n upper bounds.
"""
import os
for key in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS']:os.environ[key]='1'
from pathlib import Path
from fractions import Fraction as F
from concurrent.futures import ProcessPoolExecutor,as_completed
import json,itertools,time
import numpy as np
from scipy import sparse as ss
from scipy.optimize import linprog
HERE=Path(__file__).resolve().parent

def rational(x):return F(float(x)).limit_denominator(10**8)
def encode(x):return str(x.numerator) if x.denominator==1 else str(x.numerator)+'/'+str(x.denominator)

def geometry(boxes,common):
    n=len(boxes[0]);axes=[sorted(set(v for box in boxes for v in box[j])) for j in range(n)]
    if common:
        axes=[sorted(set(a+[0])) for a in axes]
    lengths=[len(a)-1 for a in axes];tuples=list(itertools.product(*[range(k) for k in lengths]));ind=np.array(tuples,dtype=np.int16)
    low=np.array([[axes[j][a[j]] for j in range(n)] for a in tuples],dtype=np.int64)
    high=np.array([[axes[j][a[j]+1] for j in range(n)] for a in tuples],dtype=np.int64)
    B=np.array([np.all((low>=np.array([b[0] for b in box]))&(high<=np.array([b[1] for b in box])),axis=1) for box in boxes])
    keep=np.any(B,axis=0);B=B[:,keep];ind=ind[keep];low=low[keep];high=high[keep]
    masses=np.prod(high-low,axis=1).astype(np.int64);vol=np.array([np.prod([r-l for l,r in box]) for box in boxes],dtype=np.int64)
    assert np.all(B@masses==vol)
    edges=[]
    if common:
        lookup={tuple(a):j for j,a in enumerate(ind)}
        for far,j in lookup.items():
            for axis in range(n):
                c=list(far);a=axes[axis];idx=far[axis]
                if a[idx]>=0 and idx>0 and a[idx-1]>=0:c[axis]-=1
                elif a[idx+1]<=0 and idx+2<len(a) and a[idx+2]<=0:c[axis]+=1
                else:continue
                near=lookup.get(tuple(c))
                if near is not None:edges.append((near,j))
    return B,masses,vol,edges

def solve(boxes,common,isotone):
    B,masses,vol,edges=geometry(boxes,common);M=len(masses);V=int(sum(masses));r=len(boxes)
    constraints=[-ss.csr_matrix(B*masses[None,:]/vol[:,None])];rhs=[-np.ones(r)/2]
    if isotone and edges:
        rr=np.repeat(np.arange(len(edges)),2);cc=np.array(edges).ravel();dd=np.tile([1.,-1.],len(edges));O=ss.csr_matrix((dd,(rr,cc)),shape=(len(edges),M))
        constraints.append(-O);rhs.append(np.zeros(len(edges)))
    result=linprog(masses/V,A_ub=ss.vstack(constraints),b_ub=np.concatenate(rhs),bounds=(0,1),method='highs',
       options={'primal_feasibility_tolerance':1e-9,'dual_feasibility_tolerance':1e-9})
    assert result.success,result.message
    x=[max(F(0),min(F(1),rational(v))) for v in result.x]
    if isotone:
        # Monotone majorant raises only densities. Repeated sweeps terminate on the acyclic partial order.
        for _ in range(sum(len(box) for box in [boxes[0]])*10+20):
            changed=False
            for near,far in edges:
                if x[near]<x[far]:x[near]=x[far];changed=True
            if not changed:break
        assert all(x[a]>=x[b] for a,b in edges)
    base_primal=x.copy()
    epsilon=F(0);rowmass=[]
    for i in range(r):
        mass=sum((int(masses[j])*x[j] for j in np.flatnonzero(B[i])),F(0));rowmass.append(mass);target=F(int(vol[i]),2)
        if mass<target:epsilon=max(epsilon,(target-mass)/(int(vol[i])-mass))
    epsilon+=F(1,10**10);assert epsilon<1
    x=[(1-epsilon)*a+epsilon for a in x]
    primal=sum((int(m)*a for m,a in zip(masses,x)),F(0))
    margins=[]
    for i in range(r):
        mass=sum((int(masses[j])*x[j] for j in np.flatnonzero(B[i])),F(0));margin=mass/F(int(vol[i]))-F(1,2);assert margin>0;margins.append(margin)
    multipliers=-result.ineqlin.marginals
    lam=[max(F(0),rational(V*multipliers[i]/int(vol[i]))) for i in range(r)]
    coeff=[F(int(m)) for m in masses]
    for i,l in enumerate(lam):
        if l:
            for j in np.flatnonzero(B[i]):coeff[j]-=int(masses[j])*l
    eta_values=[]
    if isotone and edges:
        for q,(near,far) in enumerate(edges):
            z=max(F(0),rational(V*multipliers[r+q]));coeff[near]-=z;coeff[far]+=z
            if z:eta_values.append([q,encode(z)])
    dual=sum((l*int(v)/2 for l,v in zip(lam,vol)),F(0))+sum((min(F(0),v) for v in coeff),F(0))
    assert 0<dual<=primal
    groups={}
    for j,v in enumerate(base_primal):
        if v:groups.setdefault(encode(v),[]).append(j)
    certificate=dict(base_primal_groups=groups,strict_mix=encode(epsilon),
        dual_box_multipliers=[[i,encode(v)] for i,v in enumerate(lam) if v],dual_order_multipliers=eta_values)
    return dict(certificate=certificate,atoms=M,rows=r,union_volume=V,order_edges=len(edges) if isotone else 0,
       lp_objective_float=float(result.fun*V),source_interval=[encode(dual),encode(primal)],
       strict_source_mean_margin=encode(min(margins)),strict_mixture=encode(epsilon),
       fixed_catalogue_score_interval=[encode(F(V)/primal),encode(F(V)/dual)],
       relative_certificate_gap=float((primal-dual)/dual),solver_success=bool(result.success),
       solver_message=str(result.message))

def run(job):
    n,K,seed,common=job;start=time.time();rng=np.random.default_rng(seed);boxes=[]
    count=12 if n<=3 else (20 if n<=4 else 28)
    for _ in range(count):
        box=[]
        for j in range(n):
            if common:
                a,b=map(int,rng.choice([0,1,K],2));
                if a==b==0:b=1
                box.append([-a,b])
            else:
                a,b=sorted(map(int,rng.choice([0,1,2,3],2,replace=False)));box.append([a,b])
        boxes.append(box)
    boxes=list({tuple(map(tuple,b)):b for b in boxes}.values())
    plain=solve(boxes,common,False);ordered=solve(boxes,common,True) if common else None
    if ordered:
        a,b=map(F,plain['source_interval']);c,d=map(F,ordered['source_interval']);assert max(a,c)<=min(b,d)
    return dict(n=n,K=K,seed=seed,common_point=common,boxes=boxes,plain=plain,isotone=ordered,seconds=time.time()-start)
if __name__=='__main__':
    jobs=[(n,K,2026111400+100*n+10*K+j,True) for n in [2,3,4,5,6] for K in [2,4,16] for j in range(2)]
    jobs +=[(n,2,2026111400+100*n+j,False) for n in [2,3,4,5,6] for j in range(2)]
    out=[]
    with ProcessPoolExecutor(12) as pool:
        for future in as_completed([pool.submit(run,j) for j in jobs]):
            r=future.result();out.append(r);print(r['n'],r['K'],r['common_point'],r['plain']['atoms'],r['plain']['relative_certificate_gap'],flush=True)
            (HERE/'catalogues.json').write_text(json.dumps(sorted(out,key=lambda r:(r['n'],not r['common_point'],r['K'],r['seed'])),indent=2)+'\n')
