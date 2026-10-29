"""Unchanged selected helpers from round038 frozen driver; no cache imports."""
import math
from fractions import Fraction as F
import numpy as np

def lse(x,axis=-1):
 m=np.max(x,axis=axis,keepdims=True);safe=np.where(np.isfinite(m),m,0.)
 with np.errstate(divide='ignore',invalid='ignore',under='ignore'):out=safe+np.log(np.exp(x-safe).sum(axis=axis,keepdims=True))
 return np.squeeze(out,axis)

def phi(z):
 z=np.asarray(z);a=np.fromiter((.5*math.erfc(-float(v)/math.sqrt(2))for v in z.ravel()),dtype=float,count=z.size)
 return a.reshape(z.shape)

def cdf_enclosure(cut,mean,var,error,eps):
 vm=var-error;vp=var+error
 assert np.min(vm)>0
 sm=np.sqrt(vm);sp=np.sqrt(vp)
 nm=cut+.5-(mean+error);np_=cut+.5-(mean-error)
 zl=np.minimum(nm/sm,nm/sp);zu=np.maximum(np_/sm,np_/sp)
 return np.maximum(0,phi(zl)-eps),np.minimum(1,phi(zu)+eps)

def log_prob(p):
 with np.errstate(divide='ignore'):return np.log(p)

def log_mixture(logl,logu,logr):
 lv=2*lse(logu+.5*logr,axis=1)
 return lse(logl,axis=1),lse(logu,axis=1),lv

def tail_enclosures(logl,logu,logv,h):
 with np.errstate(over='ignore',under='ignore',divide='ignore',invalid='ignore'):
  d=np.maximum(0,-np.expm1(h-logl));r=np.exp(logv-2*logl)
  low=np.where(logl>h,d*d/(r+d*d),0.)
  t=np.maximum(0,-np.expm1(logu-h));rv=np.exp(logv-2*h)
  up=np.where(logu<h,rv/(rv+t*t),1.)
  up=np.minimum(up,np.minimum(1,np.exp(logu-h)))
 low=np.nan_to_num(low,nan=0.);up=np.nan_to_num(up,nan=1.)
 return low.max(axis=1),np.minimum(1,up.sum(axis=1))

def bounded_ci(mean,cap,n,alpha,upper):
 x=max(0,min(1,mean/cap));target=math.log(1/alpha)/n
 if (upper and x==1)or(not upper and x==0):return cap if upper else 0.
 a,b=(x,1.)if upper else(0.,x)
 for _ in range(80):
  p=(a+b)/2
  div=(x*math.log(x/p)if x else 0)+((1-x)*math.log((1-x)/(1-p))if x<1 else 0)
  if upper:
   if div>target:b=p
   else:a=p
  else:
   if div>target:a=p
   else:b=p
 return cap*((a+b)/2)

def summarize(values,cap,alpha,upper):
 n=len(values);s=float(values.sum());ss=float(np.dot(values,values))
 if s==0:ess=0.
 else:
  positive=values[values>0];ll=np.log(positive);ess=math.exp(2*float(lse(ll))-float(lse(2*ll)))
 return {'mean':float(values.mean()),'SE_descriptive':float(values.std(ddof=1)/math.sqrt(n)),
         'ESS':ess,'CI':bounded_ci(float(values.mean()),cap,n,alpha,upper),
         'max_contribution':float(values.max())}
