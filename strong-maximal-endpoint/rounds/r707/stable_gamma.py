"""Stable lower Gamma tail from a positive Laplace integral.
Gauss-Laguerre is floating-point quadrature; certify_gamma.py supplies point certificates.
"""
import numpy as np
from scipy import special as sp
_ORIGINAL_P=sp.gammainc;_ORIGINAL_Q=sp.gammaincc
_RULES={order:sp.roots_laguerre(order) for order in [32,64]}
def log1pmx(x):
    x=np.asarray(x);small=np.abs(x)<.05
    y=np.empty_like(x);y[~small]=np.log1p(x[~small])-x[~small]
    u=x[small];acc=np.zeros_like(u)
    for k in range(16,1,-1):acc=acc*u+(-1)**(k+1)/k
    y[small]=acc*u*u
    return y
def lower(n,x,order=64):
    x=np.asarray(x,dtype=float);mask=(n>=100000)&(x<n-4*np.sqrt(n))&(x>0)
    out=np.empty_like(x);out[~mask]=_ORIGINAL_P(n,x[~mask])
    if np.any(mask):
        selected=x[mask];values=np.empty_like(selected);nodes,weights=_RULES[order]
        for start in range(0,len(selected),4096):
            xx=selected[start:start+4096];lam=(n-1-xx)/n
            u=-nodes[None,:]/(n*lam[:,None])
            correction=(n-1)*log1pmx(u)
            integral=np.exp(correction)@weights
            logpref=n*log1pmx((xx-n)/n)-.5*np.log(2*np.pi*n)-1/(12*n)+1/(360*float(n)**3)
            values[start:start+4096]=np.exp(logpref)/lam*integral
        out[mask]=values
    return float(out) if out.ndim==0 else out
def upper(n,x,order=64):
    x=np.asarray(x,dtype=float);mask=(n>=100000)&(x<n-4*np.sqrt(n))&(x>0)
    out=np.empty_like(x);out[~mask]=_ORIGINAL_Q(n,x[~mask]);out[mask]=1-lower(n,x[mask],order)
    return float(out) if out.ndim==0 else out
class Special:
    def __init__(self,order=64):self.order=order
    def __getattr__(self,name):return getattr(sp,name)
    def gammainc(self,n,x):return lower(n,x,self.order)
    def gammaincc(self,n,x):return upper(n,x,self.order)
