import numpy as np
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
