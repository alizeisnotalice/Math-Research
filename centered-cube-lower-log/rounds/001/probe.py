from pathlib import Path
import numpy as np,json,time,hashlib
P=Path(__file__).resolve().parent
start=time.perf_counter(); rows=[]; seed=20261028

def iso(x,Jlist):
    n=x.shape[1]; z=np.sort(x,axis=1); best=-n*np.log(2*z[:,-1]); out={}
    kk=np.arange(1,n+1)
    if 1 in Jlist:out[1]=best.copy()
    for j in range(1,max(Jlist)):
        d=(j+1)//2 + (z if j%2==0 else -z[:,::-1])
        score=n*np.log(j)+kk*np.log1p(1/j)-n*np.log(2*d)
        best=np.maximum(best,score.max(axis=1))
        if j+1 in Jlist:out[j+1]=best.copy()
    return out

def general(x,a,J):
    n=len(a); m=len(x); events=[];delta=[]
    for j in range(J):
        events.append(a*((j+1)//2 + (x if j%2==0 else -x)))
        delta.append(np.full((m,n),0 if j==0 else np.log1p(1/j)))
    ee=np.concatenate(events,axis=1);dd=np.concatenate(delta,axis=1)
    order=np.argsort(ee,axis=1); rr=np.take_along_axis(ee,order,axis=1)
    ll=np.cumsum(np.take_along_axis(dd,order,axis=1),axis=1)+np.log(a).sum()-n*np.log(2*rr)
    ll[rr < (x*a).max(axis=1)[:,None]]=-np.inf
    return ll.max(axis=1)

def evaluate(vals,name,n,J,sigma):
    mid=len(vals)//2;train=np.sort(vals[:mid])[::-1]; ranks=np.arange(1,mid+1)
    # infinitesimal downward threshold, explicitly frozen before validation
    pos=np.argmax(train+np.log(ranks/mid));h=float(np.nextafter(train[pos],-np.inf))
    hits=int(np.count_nonzero(vals[mid:]>h));q=hits/(len(vals)-mid)
    mean=float(np.exp(h)*q);se=float(np.exp(h)*np.sqrt(q*(1-q)/(len(vals)-mid)))
    row=dict(name=name,n=n,J=J,sigma=sigma,Ntrain=mid,Nvalidation=len(vals)-mid,log_threshold=h,threshold=float(np.exp(h)),hits=hits,estimate=mean,mcse=se,normalized_sqrt=mean/np.sqrt(n),max_log_observed=float(vals.max()))
    rows.append(row);np.save(P/(name+'.npy'),vals)

# independent direct floor-count control at every event
rng=np.random.default_rng(seed);maxerr=0.
for n in [2,3,5]:
    a=np.exp(np.linspace(-.13,.13,n));x=rng.random((9,n))*.5;got=general(x,a,6);want=[]
    for xx in x:
        rr=np.array([a[i]*((j+1)//2+(xx[i] if j%2==0 else -xx[i])) for j in range(6) for i in range(n)])
        rr=np.nextafter(rr,np.inf);center=xx*a
        cnt=np.floor((center[None,:]+rr[:,None])/a+1e-12)-np.ceil((center[None,:]-rr[:,None])/a-1e-12)+1
        valid=np.all(cnt>0,axis=1); val=np.log(cnt[valid]).sum(axis=1)+np.log(a).sum()-n*np.log(2*rr[valid]);want.append(val.max())
    maxerr=max(maxerr,float(np.max(np.abs(got-want))))
assert maxerr<1e-10,maxerr
# identical-source two algorithms
x=rng.random((30,8))*.5
isoerr=float(np.max(np.abs(iso(x,[16])[16]-general(x,np.ones(8),16))))
assert isoerr<1e-10
(P/'control.json').write_text(json.dumps(dict(direct_count_max_error=maxerr,iso_general_max_error=isoerr),indent=2))
for n in [8,32,128,512,2048]:
    rng=np.random.default_rng(np.random.SeedSequence([seed,n,0]));parts={J:[] for J in [2,4,8,16,32,64]}
    for _ in range(32):
        got=iso(rng.random((256,n))*.5,list(parts))
        for J in parts:parts[J].append(got[J])
    for J in parts:evaluate(np.concatenate(parts[J]),f'iso_n{n}_J{J}',n,J,0)
    print('isotropic',n,flush=True)
for n in [8,32,128,512]:
    for si,sigma in enumerate([0,.02,.08,.2]):
        rng=np.random.default_rng(np.random.SeedSequence([seed,n,si+1]));a=np.exp(sigma*np.r_[-np.ones(n//2),np.ones(n//2)]);parts=[]
        for _ in range(64):parts.append(general(rng.random((64,n))*.5,a,16))
        evaluate(np.concatenate(parts),f'block_n{n}_s{si}',n,16,sigma)
    print('block',n,flush=True)
(P/'results.json').write_text(json.dumps(dict(seed=seed,numpy=np.__version__,seconds=time.perf_counter()-start,scope='finite-event periodic lattice Monte Carlo; not certified weak constant',rows=rows),indent=2))
print('seconds',time.perf_counter()-start)
