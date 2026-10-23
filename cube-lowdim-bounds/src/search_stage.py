"""Bounded discovery only. Float scores are NOT proof certificates."""
from pathlib import Path
import sys, json, time, itertools, functools
import numpy as np
from scipy.optimize import minimize, differential_evolution
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'baseline'))
from weighted2 import score as weighted_score
from uniform import score as uniform_score

def union(boxes):
    @functools.lru_cache(None)
    def rec(B):
        if not B:return 0.
        if len(B[0])==2:
            a,b=B[0];v=0.
            for c,e in B[1:]:
                if c>b:v+=b-a;a,b=c,e
                else:b=max(b,e)
            return v+b-a
        cuts=sorted({x for b in B for x in b[:2]})
        return sum((b-a)*rec(tuple(sorted({c[2:] for c in B if c[0]<=a and b<=c[1]}))) for a,b in zip(cuts,cuts[1:]))
    return rec(tuple(sorted(set(boxes))))

def prune(B):
    B=sorted(set(B))
    if not B:return []
    A=np.array(B);keep=[]
    for i,b in enumerate(A):
        con=np.all(A[:,::2]<=b[::2],axis=1)&np.all(A[:,1::2]>=b[1::2],axis=1)
        con[i]=False
        if not con.any():keep.append(tuple(b))
    return keep

def finite_boxes(z,w,closed_only=False):
    n,d=z.shape;out=[]
    for mask in range(1,1<<n):
        ix=[i for i in range(n) if mask>>i&1];zz=z[ix]
        lo=zz.min(0);hi=zz.max(0)
        if closed_only:
            inc=np.all(z>=lo,axis=1)&np.all(z<=hi,axis=1)
            if sum(1<<i for i in np.where(inc)[0])!=mask:continue
        r=w[ix].sum()**(1/d)/2
        a=hi-r;b=lo+r
        if np.all(a<b):out.append(tuple(np.column_stack((a,b)).flat))
    return prune(out)

def finite_score(p,n):
    z=np.r_[np.zeros(2),p[:2*(n-1)]].reshape(n,2)
    w=np.r_[1.,np.exp(p[2*(n-1):])]
    return union(finite_boxes(z,w))/w.sum()

@functools.lru_cache(None)
def templates(d,L):
    out=[]
    for ls in itertools.product(range(1,L+1),repeat=d):
        for q in itertools.product((0,1),repeat=2):
            starts=[-(ls[j]//2)+q[j] for j in range(2)]+[(-l//2+1 if l%2 else 1-l//2) for l in ls[2:]]
            # odd l: -(l-1)/2; even l: 1-l/2
            starts[2:]=[-(l//2) if l%2 else 1-l//2 for l in ls[2:]]
            cnt=[]
            for j in range(2):cnt.append([sum((s%2)==p for s in range(starts[j],starts[j]+ls[j])) for p in (0,1)])
            coeff=np.outer(*cnt).reshape(-1)*np.prod(ls[2:],dtype=int)
            out.append((ls,starts,coeff))
    return np.array([x[0] for x in out]),np.array([x[1] for x in out]),np.array([x[2] for x in out])

def checker_score(p,d,L=4):
    a=np.asarray(p[:d]);t=np.exp(p[d]);W=np.array([1.,t,t,1.])
    if np.any(a<=0):return 0.
    ls,starts,C=templates(d,L);mass=C@W;r=mass**(1/d)/2
    limits=a.copy();limits[2:]/=2
    lo=np.maximum(0,(starts+ls-1)*a-r[:,None]);hi=np.minimum(limits,starts*a+r[:,None])
    ok=np.all(lo<hi,axis=1)
    B=[tuple(np.column_stack((x,y)).flat) for x,y in zip(lo[ok],hi[ok])]
    return 2**d*union(prune(B))/W.sum()

def linparams():
    s=(2+np.sqrt(22))/3;h=(1+s)/2
    return [h,h,h+1,s*s-1]

def baseparams(d):
    u=2**(1/d);v=3**(1/d)
    return [1+u**(d-2)/2,(2+u**(d-1)*v)/4]+[(1+u)/2]*(d-2)

def run(round_no):
    seed=20261023+round_no;rng=np.random.default_rng(seed);records=[];start=time.time()
    def save(label,score,p,**kw):
        row={'label':label,'score':float(score),'parameters':list(map(float,p)),**kw};records.append(row)
        print(json.dumps(row),flush=True)
        (ROOT/'results'/f'round{round_no}_search.json').write_text(json.dumps({'seed':seed,'elapsed_seconds':time.time()-start,'records':records},indent=2))
    lp=linparams();save('Lin_reproduction',weighted_score(lp),lp)
    # Every round changes the discovery model/budget, with a fixed seed.
    if round_no==1:
        for n in (3,4,5):
            bounds=[(-2,2)]*(2*(n-1))+[(-1.5,1.5)]*(n-1)
            res=differential_evolution(lambda p:-finite_score(p,n),bounds,maxiter=24,popsize=5,seed=seed+n,polish=False)
            rr=minimize(lambda p:-finite_score(p,n),res.x,method='Nelder-Mead',options={'maxiter':350})
            save('finite_atoms',-rr.fun,rr.x,n=n,nfev=int(res.nfev+rr.nfev))
    else:
        typ='2x1' if round_no==2 else '2x2'
        p0=lp if round_no==2 else [lp[0],lp[1],lp[2],lp[2],1.,lp[3],lp[3]]
        best=(weighted_score(p0,typ,K=6),p0);evals=0
        for k in range(3):
            x=np.array(p0)*(1+rng.normal(0,.025 if k else 0,len(p0)))
            def obj(p):
                return 1e3 if min(p)<=0 else -weighted_score(p,typ,K=6)
            rr=minimize(obj,x,method='Nelder-Mead',options={'maxiter':400,'xatol':1e-9,'fatol':1e-11});evals+=rr.nfev
            if -rr.fun>best[0]:best=(-rr.fun,rr.x)
        save('joint_periodic_'+typ,*best,nfev=int(evals),K=6)
    for d in (3,4):
        a=baseparams(d);L=4+round_no
        save('expanded_uniform_witnesses',uniform_score(a,d,True,L=L),a,d=d,L=L)
        p0=a+[0.];best=(checker_score(p0,d,L=4),p0);trials=0
        scale=[.01,.035,.1][round_no-1]
        for k in range(24):
            p=np.array(p0);p[:d]*=np.exp(rng.normal(0,scale,d));p[d]=rng.uniform(-2*scale,2*scale)
            z=checker_score(p,d,L=4);trials+=1
            if z>best[0]:best=(z,p)
        save('nonproduct_checkerboard',*best,d=d,L=4,trials=trials)
    print('DONE',round_no,'seconds',time.time()-start,flush=True)

if __name__=='__main__':run(int(sys.argv[1]))
