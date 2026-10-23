import itertools,sys,json,time
import numpy as np
from scipy.optimize import differential_evolution,minimize
from search import union_volume

def wraps(a,b,P):
 if b<=a:return []
 if b-a>=P:return [(0.,P)]
 a=a%P;b=a+b-(b-(b-a)) if False else b
 # reconstruct length before mod in caller version

def split(a,b,P):
 l=b-a
 if l<=0:return []
 if l>=P:return [(0.,P)]
 a=a%P
 return [(a,a+l)] if a+l<=P else [(a,P),(0.,a+l-P)]

def getboxes(p,typ='2x1',K=5):
 if typ=='2x1':
  a,b,v,w=p;xs=[0.,a];ys=[0.];P=a+b;Q=v;W=np.array([[1.],[w]])
 elif typ=='2x2':
  a,b,c,e,w1,w2,w3=p;xs=[0.,a];ys=[0.,c];P=a+b;Q=c+e;W=np.array([[1.,w1],[w2,w3]])
 if min(p)<=0:return [],1
 nx=len(xs);ny=len(ys);B=[]
 for s,t in itertools.product(range(nx),range(ny)):
  for k,l in itertools.product(range(1,K+1),repeat=2):
   x1=xs[s];x2=xs[(s+k-1)%nx]+P*((s+k-1)//nx)
   y1=ys[t];y2=ys[(t+l-1)%ny]+Q*((t+l-1)//ny)
   mass=sum(W[(s+i)%nx,(t+j)%ny] for i in range(k) for j in range(l))
   r=mass**.5/2
   if x2-x1>=2*r or y2-y1>=2*r:continue
   for X,Y in itertools.product(split(x2-r,x1+r,P),split(y2-r,y1+r,Q)):B.append(X+Y)
 B=list(set(B))
 if B:
  A=np.array(B);keep=[]
  for i,b in enumerate(A):
   con=np.all(A[:,::2]<=b[::2],axis=1)&np.all(A[:,1::2]>=b[1::2],axis=1);con[i]=False
   if not con.any():keep.append(tuple(b))
  B=keep
 return B,W.sum()

def score(p,typ='2x1',K=5):
 B,m=getboxes(p,typ,K)
 return union_volume(B)/m
