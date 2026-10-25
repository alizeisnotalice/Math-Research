"""Actual two-coset superperiod, old parity keys always retained."""
from pathlib import Path
import time,numpy as np
from scipy.stats import qmc
from round008_embedding import baseline_keys
from round007_union_fast import volume
from staggered_search import A
BASE=baseline_keys()

def decode(x):
    P=A*np.exp(np.asarray(x[:4]));epsilon=float(x[4]);h=np.array(x[5:9]);Q=P.copy();Q[0]*=2
    phase=np.r_[.5+h[0]/2,h[1:]]
    return P,Q,np.array([np.zeros(4),phase]),np.array([1+epsilon,1-epsilon])

def capture(Q,phases,weights,t):
    events=[]
    for orbit in range(2):
        for j in range(4):events.extend(np.abs((t[:,j,None]-phases[orbit,j]-np.arange(-3,5))*Q[j]).T)
    rr=np.array(events).T;ll=np.ceil(t[:,None,None,:]-phases[None,None,:,:]-rr[:,:,None,None]/Q-1e-11).astype(np.int16);hh=np.floor(t[:,None,None,:]-phases[None,None,:,:]+rr[:,:,None,None]/Q+1e-11).astype(np.int16)
    count=np.prod(np.maximum(0,hh-ll+1),axis=3);mass=count@weights;good=(mass>np.maximum((2*rr)**4,1e-30)*(1+1e-10))&(rr<=3);keys=set()
    for i,j in np.argwhere(good):
        lo=ll[i,j].copy();hi=hh[i,j].copy()
        for orbit in range(2):
            if np.any(hi[orbit]<lo[orbit]):lo[orbit]=0;hi[orbit]=-1
        keys.add(tuple(map(int,np.r_[lo.ravel(),hi.ravel()])))
    return keys

def boxes(Q,phases,weights,keys):
    ordered=sorted(keys);data=np.array(ordered,int);ll=data[:,:8].reshape(-1,2,4);hh=data[:,8:].reshape(-1,2,4);counts=np.prod(np.maximum(0,hh-ll+1),axis=2);mass=countsum=counts@weights;r=mass**.25/2
    low=(ll+phases)*Q;high=(hh+phases)*Q;low[counts==0]=np.inf;high[counts==0]=-np.inf
    lo=np.maximum(0,np.max(high,axis=1)-r[:,None]);hi=np.minimum(Q,np.min(low,axis=1)+r[:,None]);ok=np.all(hi>lo+1e-14,axis=1);ids=np.flatnonzero(ok);b=np.stack([lo[ok],hi[ok]],axis=2).reshape(-1,8);b,ii=np.unique(b,axis=0,return_index=True);ids=ids[ii]
    lows=b[:,::2];highs=b[:,1::2];keep=[]
    for i in range(len(b)):
        inside=np.all(lows<=lows[i],axis=1)&np.all(highs>=highs[i],axis=1);inside[i]=False
        if not inside.any():keep.append(i)
    return b[keep],[ordered[int(ids[i])] for i in keep]

def evaluate(x,pool=(),power=8,seed=20261025808,capturing=True):
    st=time.monotonic();P,Q,ph,w=decode(x);keys=BASE|set(pool)
    if capturing:keys|=capture(Q,ph,w,qmc.Sobol(4,scramble=True,seed=seed).random_base2(power))
    b,selected=boxes(Q,ph,w,keys);value=volume(b)[0]/2
    return dict(x=list(map(float,x)),P=P.tolist(),Q=Q.tolist(),epsilon=float(x[4]),h=list(map(float,x[5:])),physical_d=(P*np.array(x[5:])).tolist(),mass=2,ratio=float(value),key_count=len(keys),selected_boxes=len(b),centers=2**power if capturing else 0,seconds=time.monotonic()-st),keys,selected
