"""Fourcoset chess input; all cosets move by sign*s*d, full mass4."""
import time,numpy as np
from scipy.stats import qmc
from staggered_search import A
from round011_embedding import baseline_keys,COSETS,SIGNS
from round007_union_fast import volume
BASE=baseline_keys();S=np.array(SIGNS);OFFSETS=np.array([[p/2,q/2,0,0] for p,q in COSETS])

def decode(x):
    P=A*np.exp(np.array(x[:4]));Q=P*np.array([2,2,1,1]);epsilon=float(x[4]);h=np.array(x[5:9]);u=h/np.array([2,2,1,1]);phases=OFFSETS+S[:,None]*u;weights=1+S*epsilon
    return P,Q,phases,weights

def capture(Q,phases,weights,t):
    events=[]
    for orbit in range(4):
        for j in range(4):events.extend(np.abs((t[:,j,None]-phases[orbit,j]-np.arange(-3,5))*Q[j]).T)
    rr=np.array(events).T;ll=np.ceil(t[:,None,None,:]-phases[None,None,:,:]-rr[:,:,None,None]/Q-1e-11).astype(np.int16);hh=np.floor(t[:,None,None,:]-phases[None,None,:,:]+rr[:,:,None,None]/Q+1e-11).astype(np.int16)
    counts=np.prod(np.maximum(0,hh-ll+1),axis=3);mass=counts@weights;good=(mass>np.maximum((2*rr)**4,1e-30)*(1+1e-10))&(rr<=3);keys=set()
    for i,j in np.argwhere(good):
        l=ll[i,j].copy();h=hh[i,j].copy()
        for orbit in range(4):
            if np.any(h[orbit]<l[orbit]):l[orbit]=0;h[orbit]=-1
        keys.add(tuple(map(int,np.r_[l.ravel(),h.ravel()])))
    return keys

def boxes(Q,phases,weights,keys):
    order=sorted(keys);data=np.array(order,int);ll=data[:,:16].reshape(-1,4,4);hh=data[:,16:].reshape(-1,4,4);counts=np.prod(np.maximum(0,hh-ll+1),axis=2);mass=counts@weights;r=mass**.25/2
    low=(ll+phases)*Q;high=(hh+phases)*Q;low[counts==0]=np.inf;high[counts==0]=-np.inf
    lo=np.maximum(0,high.max(axis=1)-r[:,None]);hi=np.minimum(Q,low.min(axis=1)+r[:,None]);ok=np.all(hi>lo+1e-14,axis=1);ids=np.flatnonzero(ok);b=np.stack([lo[ok],hi[ok]],axis=2).reshape(-1,8);b,ii=np.unique(b,axis=0,return_index=True);ids=ids[ii];lows=b[:,::2];highs=b[:,1::2];keep=[]
    for i in range(len(b)):
        contained=np.all(lows<=lows[i],axis=1)&np.all(highs>=highs[i],axis=1);contained[i]=False
        if not contained.any():keep.append(i)
    return b[keep],[order[int(ids[i])] for i in keep]

def evaluate(x,pool=(),capturing=True,power=8,seed=20261025111):
    st=time.monotonic();P,Q,phases,weights=decode(x);keys=BASE|set(pool)
    if capturing:keys|=capture(Q,phases,weights,qmc.Sobol(4,scramble=True,seed=seed).random_base2(power))
    b,selected=boxes(Q,phases,weights,keys);v=volume(b)[0]/4
    return dict(x=list(map(float,x)),P=P.tolist(),Q=Q.tolist(),epsilon=float(x[4]),h=list(map(float,x[5:])),physical_d=(P*np.array(x[5:])).tolist(),mass=4,ratio=float(v),key_count=len(keys),selected_boxes=len(b),centers=2**power if capturing else 0,seconds=time.monotonic()-st),keys,selected
