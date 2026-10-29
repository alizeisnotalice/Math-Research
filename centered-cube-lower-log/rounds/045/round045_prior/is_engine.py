"""Replaceable positive-source whole-vector IS interface, independent of the box bank.
source.sample_1d(rng, shape, bank_id, n): iid source coordinates in [0,1).
source.average_1d(x, L, bank_id, n): actual positive periodic 1d window density.
These two methods must describe the same normalized Lebesgue source.
"""
import numpy as np
from scipy.special import logsumexp
class WholeVectorIS:
 def __init__(self,source,bank_ids,n,lengths,uniform_weight=.2,anchor_weights=None):
  self.source=source;self.ids=list(bank_ids);self.n=n;self.lengths=list(lengths);self.gamma=float(uniform_weight)
  assert 0<self.gamma<1 and min(self.lengths)>0
  R=len(self.ids);J=len(self.lengths)
  self.weights=np.full((R,J),(1-self.gamma)/(R*J))if anchor_weights is None else np.asarray(anchor_weights,float)
  assert self.weights.shape==(R,J) and np.min(self.weights)>0 and abs(self.weights.sum()+self.gamma-1)<1e-12
 def draw(self,rng,N):
  R,J=self.weights.shape;c=rng.choice(1+R*J,size=N,p=np.r_[self.gamma,self.weights.ravel()]);x=rng.random((N,self.n))
  for cc in np.unique(c):
   if cc==0:continue
   ids=np.flatnonzero(c==cc);b,j=divmod(int(cc)-1,J);shape=(len(ids),self.n);y=self.source.sample_1d(rng,shape,self.ids[b],self.n);assert y.shape==shape
   x[ids]=(y+self.lengths[j]*(rng.random(shape)-.5))%1
  return x,c
 def scores(self,x):
  scores=np.empty((len(x),len(self.ids),len(self.lengths)))
  for b,bankid in enumerate(self.ids):
   for j,L in enumerate(self.lengths):
    aa=self.source.average_1d(x,L,bankid,self.n);assert aa.shape==x.shape and np.min(aa)>0 and np.isfinite(aa).all();scores[:,b,j]=np.log(aa).sum(1)
  return scores
 def log_density(self,scores):return np.logaddexp(np.log(self.gamma),logsumexp(scores+np.log(self.weights)[None,:,:],axis=(1,2)))
 def integrand_cap(self,pi):
  pi=np.asarray(pi,float);assert len(pi)==len(self.ids)and np.min(pi)>=0 and abs(pi.sum()-1)<1e-12
  return float(np.max(pi[:,None]/self.weights))
