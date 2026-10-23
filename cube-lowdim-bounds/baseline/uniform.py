import sys,itertools,json,time
import numpy as np
from scipy.optimize import differential_evolution,minimize
from search import union_volume

def boxes(spacings,L=4):
 d=len(spacings);B=[]
 for ls in itertools.product(range(1,L+1),repeat=d):
  r=np.prod(ls)**(1/d)/2;b=[]
  for l,v in zip(ls,spacings):
   if l%2:a=0;c=min(v/2,r-(l-1)*v/2)
   else:a=max(0,l*v/2-r);c=v/2
   if c<=a:break
   b.extend((a,c))
  else:B.append(tuple(b))
 B=list(set(B))
 if B:
  A=np.array(B);keep=[]
  for i,b in enumerate(A):
   con=np.all(A[:,::2]<=b[::2],axis=1)&np.all(A[:,1::2]>=b[1::2],axis=1);con[i]=False
   if not con.any():keep.append(tuple(b))
  B=keep
 return B

def score(p,d,aniso=False,L=4):
 s=list(p) if aniso else [p[0]]+[p[1]]*(d-1)
 return 2**d*union_volume(boxes(s,L))
