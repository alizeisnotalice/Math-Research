import itertools,functools,time,json,sys
import numpy as np
from scipy.optimize import differential_evolution, minimize

def boxes(d,h,vs,w,K=5,L=4):
    limits=(h,)+tuple(v/2 for v in vs)
    out=[]
    for k in range(1,K+1):
      for s in range(-k,2):
       wx=sum(1 if i%2==0 else w for i in range(s,s+k))
       for ls in itertools.product(range(1,L+1),repeat=d-1):
        mass=wx*np.prod(ls)
        r=mass**(1/d)/2
        lo=max(0,(s+k-1)*h-r); hi=min(h,s*h+r)
        if hi<=lo:continue
        b=[lo,hi]
        for l,v in zip(ls,vs):
         if l%2: a=0;c=min(v/2,r-(l-1)*v/2)
         else:a=max(0,l*v/2-r);c=v/2
         if c<=a:break
         b.extend((a,c))
        else:out.append(tuple(b))
    out=list(set(out))
    # remove boxes contained in another (cost acceptable for small list)
    if out:
      a=np.array(out);keep=[]
      for i,b in enumerate(a):
       contained=np.all(a[:,::2]<=b[::2],axis=1)&np.all(a[:,1::2]>=b[1::2],axis=1)
       contained[i]=False
       if not contained.any():keep.append(tuple(b))
      out=keep
    return out

def union_volume(bs):
    @functools.lru_cache(maxsize=None)
    def rec(B):
      if not B:return 0.
      if len(B[0])==2:
       intervals=sorted(B);z=0.;a,b=intervals[0]
       for c,e in intervals[1:]:
        if c>b:z+=b-a;a=c;b=e
        else:b=max(b,e)
       return z+b-a
      cuts=sorted(set(t for b in B for t in b[:2]));z=0.
      for a,b in zip(cuts,cuts[1:]):
       mid=(a+b)/2
       S=tuple(sorted(set(t[2:] for t in B if t[0]<mid<t[1])))
       if S:z+=(b-a)*rec(S)
      return z
    return rec(tuple(sorted(bs)))

def score(p,d,K=5,L=4,anisotropic=False):
    if anisotropic:h,*vs,w=p
    else:h,v,w=p;vs=[v]*(d-1)
    B=boxes(d,h,vs,w,K,L)
    return 2**d*union_volume(B)/(1+w)
