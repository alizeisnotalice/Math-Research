"""Numerical full joint union for finite rectangular lattice witnesses."""
import numpy as np
from itertools import product
from functools import lru_cache
from scipy.optimize import minimize
import json,time
from pathlib import Path

def boxes(a,cap):
 d=len(a);l=np.array(list(product(range(1,cap+1),repeat=d)))
 r=np.prod(l,axis=1)**(1/d)/2
 odd=l%2==1
 lo=np.where(odd,0,np.maximum(0,l*a/2-r[:,None]))
 hi=np.where(odd,np.minimum(a/2,r[:,None]-(l-1)*a/2),a/2)
 ok=np.all(hi>lo,axis=1);l=l[ok];lo=lo[ok];hi=hi[ok]
 b=np.stack([lo,hi],axis=2).reshape(-1,2*d)
 b,idx=np.unique(b,axis=0,return_index=True);l=l[idx];lo=b[:,::2];hi=b[:,1::2]
 keep=[]
 for i,x in enumerate(b):
  contains=np.all(lo<=lo[i],axis=1)&np.all(hi>=hi[i],axis=1)
  contains[i]=False
  if not contains.any():keep.append(i)
 return b[keep],l[keep]

def volume(box):
 @lru_cache(None)
 def rec(B):
  if not B:return 0.
  if len(B[0])==2:
   a,b=B[0];s=0.
   for c,e in B[1:]:
    if c>b:s+=b-a;a,b=c,e
    else:b=max(b,e)
   return s+b-a
  cuts=sorted({x for b in B for x in b[:2]})
  return sum((v-u)*rec(tuple(sorted({b[2:] for b in B if b[0]<=u and b[1]>=v}))) for u,v in zip(cuts,cuts[1:]))
 return rec(tuple(sorted(set(map(tuple,box)))))

def value(a,cap=4):
 if min(a)<.8 or max(a)>2.2:return 0.
 b,l=boxes(np.asarray(a),cap)
 return 2**len(a)*volume(b)

def main():
 out=[];start=time.monotonic()
 for n in [3,4]:
  u=2**(1/n);v=3**(1/n)
  a=np.array([1+u**(n-2)/2,(2+u**(n-1)*v)/4]+[(1+u)/2]*(n-2))
  ext=[{'cap':c,'value':value(a,c),'boxes':len(boxes(a,c)[0])} for c in [4,6,8]]
  print(n,ext,flush=True)
  rng=np.random.default_rng(20261025+n);runs=[]
  for i in range(4):
   x=a if i==0 else a*(1+rng.normal(0,.025,n))
   r=minimize(lambda z:-value(z),x,method='Nelder-Mead',options={'maxfev':400,'xatol':1e-10,'fatol':1e-12})
   row={'start':x.tolist(),'a':r.x.tolist(),'value':-r.fun,'nfev':r.nfev,'message':r.message}
   runs.append(row);print(n,i,row,flush=True)
  out.append({'dimension':n,'seed':a.tolist(),'extension':ext,'runs':runs})
 result={'elapsed':time.monotonic()-start,'status':'numerical_only','dimensions':out}
 (Path(__file__).resolve().parents[1]/'results/round002_search.json').write_text(json.dumps(result,indent=2)+'\n')
if __name__=='__main__':main()
