"""True 3D body-diagonal Cauchy heat/Green formulas, numerical interface.

Fixed V=(1,+/-1,+/-1), standard Fourier inverse (2*pi)^-3.
Default Green uses a strictly positive rational integrand, not termwise
integration of divergent cone terms. Adaptive quadrature errors are numerical,
not certified intervals. All nonzero body-diagonal lines return +infinity.
"""
import math
import numpy as np
from scipy.integrate import quad

def cones():
    out=[]
    for i in range(3):
        j,k=[v for v in range(3) if v!=i]
        for sign in (-1,1):
            def ray(axis,ss):
                r=[0,0,0];r[i]=sign;r[axis]=ss;return tuple(r)
            a,c=ray(j,1),ray(j,-1);b,d=ray(k,1),ray(k,-1)
            out.extend([(a,b,c),(a,c,d)])
    from itertools import product
    for e in product((-1,1),repeat=3):out.append(((e[0],e[1],0),(e[0],0,e[2]),(0,e[1],e[2])))
    return tuple(out)
CONES=cones()

def _data(x):
    x=np.asarray(x,dtype=float)
    if x.shape!=(3,) or not np.all(np.isfinite(x)):raise ValueError('x must be a finite real 3-vector')
    scale=float(np.max(np.abs(x)))
    if scale==0:return scale,None
    x=np.sort(np.abs(x/scale))[::-1]
    X,Y,Z=x*x;S1=X+Y+Z;S2=X*Y+Y*Z+Z*X
    # SOS identity S1*S2-9*X*Y*Z avoids catastrophic cancellation near a line.
    C=X*(Y-Z)**2+Y*(Z-X)**2+Z*(X-Y)**2
    A=np.array([(x[1]-x[2])**2,(x[1]+x[2])**2,(x[0]-x[1])**2,(x[0]-x[2])**2,(x[0]+x[2])**2,(x[0]+x[1])**2])
    return scale,(S1,S1*S1+5*S2,C,A)

def _integrand(w,S,T,C,A):
    v=1-w;zeros=int(np.count_nonzero(A==0));aa=A[A!=0]
    den=math.prod(w+float(a)*v for a in aa)
    if zeros==0:num=w*(20*w**3+24*S*w*w*v+4*T*w*v*v+4*C*v**3)
    elif zeros==1:num=20*w**3+24*S*w*w*v+4*T*w*v*v+4*C*v**3
    elif zeros==2:num=20*w*w+24*S*w*v+4*T*v*v # coordinate axis; C=0 exactly
    else:return math.inf
    return num/den

def green(x, *, method='auto', dps=50, epsrel=2e-10, epsabs=2e-10, return_error=False):
    """G(x)=integral_0^infty p_t(x)dt, float by default; no certified error.

    method='mp' uses mpmath and returns mpf (optional error is None).
    Origin and all |x1|=|x2|=|x3| nonzero directions are true +inf.
    """
    if method=='mp':return _green_mp(x,dps=dps,return_error=return_error)
    if method not in ('auto','quad'):raise ValueError('method must be auto, quad or mp')
    scale,data=_data(x)
    if scale==0 or int(np.count_nonzero(data[3]==0))>=3:
        return (math.inf,math.inf) if return_error else math.inf
    S,T,C,A=data
    points=sorted(set(float(a/(1+a)) for a in A if a>0 and a/(1+a)<1))
    value,error=quad(_integrand,0.,1.,args=(S,T,C,A),points=points,epsrel=epsrel,epsabs=epsabs,limit=250)
    normal=1/(32*math.pi**3*scale**2)
    return (value*normal,error*normal) if return_error else value*normal

def green_face(y,z,**kwargs):
    """G(1,y,z); scalar interface for a face-cost driver, supports arrays too."""
    if np.ndim(y)==0 and np.ndim(z)==0:return green((1.,y,z),**kwargs)
    y,z=np.broadcast_arrays(y,z);return green_array(np.stack([np.ones_like(y),y,z],axis=-1),**kwargs)

def green_array(x,**kwargs):
    x=np.asarray(x,dtype=float)
    if x.shape[-1]!=3:raise ValueError('last dimension must be 3')
    return np.array([green(v,**kwargs) for v in x.reshape(-1,3)]).reshape(x.shape[:-1])

def heat_kernel(t,x, *, method='rational'):
    """p_t for t>0; rational form is positive, cones retains all complex terms."""
    t=float(t);x=np.asarray(x,dtype=float)
    if t<=0 or x.shape!=(3,):raise ValueError('t>0 and real 3-vector required')
    if method=='cones':
        return (sum(1/math.prod(complex(4*t,-float(np.dot(x,r))) for r in c) for c in CONES)/(4*math.pi**3)).real
    scale=max(4*t,float(np.max(np.abs(x))));xx=x/scale;a=4*t/scale;u=a*a
    X,Y,Z=xx*xx;S=X+Y+Z;S2=X*Y+Y*Z+Z*X;C=X*(Y-Z)**2+Y*(Z-X)**2+Z*(X-Y)**2
    L=[xx[1]-xx[2],xx[1]+xx[2],xx[0]-xx[1],xx[0]-xx[2],xx[0]+xx[2],xx[0]+xx[1]]
    N=u*(20*u**3+24*S*u*u+4*(S*S+5*S2)*u+4*C)
    return a*N/math.prod(u+l*l for l in L)/(4*math.pi**3*scale**3)

def _mp(v,mp):
    if hasattr(v,'numerator') and hasattr(v,'denominator'):return mp.mpf(int(v.numerator))/int(v.denominator)
    return mp.mpf(str(v))

def _green_mp(x, *, dps=50, return_error=False):
    import mpmath as mp
    with mp.workdps(dps):
        x=sorted([abs(_mp(v,mp)) for v in x],reverse=True);scale=x[0]
        if not scale or x[0]==x[1]==x[2]:return (mp.inf,None) if return_error else mp.inf
        x=[v/scale for v in x];X,Y,Z=[v*v for v in x];S=X+Y+Z;T=S*S+5*(X*Y+Y*Z+Z*X);C=X*(Y-Z)**2+Y*(Z-X)**2+Z*(X-Y)**2
        A=[(x[1]-x[2])**2,(x[1]+x[2])**2,(x[0]-x[1])**2,(x[0]-x[2])**2,(x[0]+x[2])**2,(x[0]+x[1])**2];zero=sum(a==0 for a in A);aa=[a for a in A if a]
        def f(w):
            v=1-w;den=mp.fprod(w+a*v for a in aa)
            if zero==0:num=w*(20*w**3+24*S*w*w*v+4*T*w*v*v+4*C*v**3)
            elif zero==1:num=20*w**3+24*S*w*w*v+4*T*w*v*v+4*C*v**3
            else:num=20*w*w+24*S*w*v+4*T*v*v
            return num/den
        pts=[mp.mpf(0)]+sorted(set(a/(1+a) for a in A if a))+[mp.mpf(1)]
        value=mp.quad(f,pts)/(32*mp.pi**3*scale**2)
        return (value,None) if return_error else value
