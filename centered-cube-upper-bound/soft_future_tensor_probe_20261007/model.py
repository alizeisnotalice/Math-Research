"""Exact finite-grid geometric aggregation, bounded scalar quadrature.
Bounds include an analytic quadrature remainder and a heuristic floating guard;
there is NO rigorous floating-point interval certificate.
"""
import math
import numpy as np
Q=(1,2,4);P=(.25,.5,.25)

def geo(a,d,m,v):
    # Sum exp(-(a+kd)/v), k=0,...,m-1.
    return np.where(m>0,np.exp(-a/v)*(-np.expm1(-m*d/v))/(-np.expm1(-d/v)),0.)

def primedsum(a,d,m):
    # Upper on sum (1+a+kd) exp(-(a+kd)), using exact geometric formula.
    r=np.exp(-d);rm=np.exp(-m*d);den=-np.expm1(-d)
    g=np.exp(-a)*(-np.expm1(-m*d))/den
    rank=np.exp(-a)*(r-m*rm+(m-1)*rm*r)/(den*den)
    return np.where(m>0,(1+a)*g+d*rank,0.)

def coordinates(Z,dx,Ls,N,block=8):
    n=len(Z);Ls=np.array(Ls);hc=[];soft=[];err=[];curv=[]
    for q in Q:
        count=64*q;step=4//q
        diff=dx[:,None]+(Z[:,None]-step*np.arange(count)[None,:])/4
        hard=np.empty((len(Ls),n));phi=np.empty_like(hard);delta=np.empty_like(hard);bounds=np.empty_like(hard)
        for begin in range(0,len(Ls),block):
            L=Ls[begin:begin+block,None]
            inside=2*np.abs(diff)[None,:,:]<=L[:,:,None]
            c=inside.sum(axis=2);left=(diff[None,:,:]>L[:,:,None]/2).sum(axis=2);right=(diff[None,:,:]<-L[:,:,None]/2).sum(axis=2)
            assert np.all(c+left+right==count)
            # Components are aligned arithmetic grids; all hard boundary groups
            # are counted closed, including exact local tagged source faces.
            lo=left;hi=count-right-1;d=1/(q*L)
            def at(k):return (dx[None,:]+(Z[None,:]-step*k)/4)/L
            ai=np.maximum(0,.5+at(hi));bi=np.maximum(0,.5-at(lo))
            al=np.maximum(0,at(left-1)-.5);ar=np.maximum(0,-at(count-right)-.5)
            curvature=primedsum(ai,d,c)+primedsum(bi,d,c)+primedsum(al,d,left)+primedsum(al+1,d,left)+primedsum(ar,d,right)+primedsum(ar+1,d,right)
            integral=np.zeros(c.shape)
            for vb in range(1,N+1,32):
                vv=np.arange(vb,min(N+1,vb+32))/N
                v=vv[None,None,:];weight=np.ones(len(vv));weight[vv==1]=.5
                C=c[:,:,None];D=d[:,:,None]
                inner=2*C-geo(ai[:,:,None],D,C,v)-geo(bi[:,:,None],D,C,v)
                outer=(-np.expm1(-1/v))*(geo(al[:,:,None],D,left[:,:,None],v)+geo(ar[:,:,None],D,right[:,:,None],v))
                integral+=np.sum(v*(inner+outer)*weight,axis=-1)/N
            # f(0)=0; trapezoid Peano bound h²/8 int|f''|. At a=0,
            # g_a(v)=v and curvature=0; (1+a)e^-a=1 remains a safe upper.
            quad=curvature/(8*N*N)
            hard[begin:begin+block]=c/(q*L);phi[begin:begin+block]=integral/(q*L)
            # Numerical guard is documented, not a verified rounding bound.
            delta[begin:begin+block]=(quad+1e-10)/(q*L);bounds[begin:begin+block]=curvature/(q*L)
        hc.append(hard);soft.append(phi);err.append(delta);curv.append(bounds)
    return np.array(hc),np.array(soft),np.array(err),np.array(curv)

def response(h,b,s):
    a=(1-s)*h+s*b
    with np.errstate(divide='ignore'):
        logs=np.log(np.maximum(a,0.)).sum(axis=-1)+np.log(np.array(P))[:,None]
    return np.exp(np.logaddexp.reduce(logs,axis=0))

def sigma_cell_upper(h,blo,bhi,a,b):
    # Scalar coordinate responses are affine in sigma; endpoint maximum
    # gives a positive product upper for an entire sigma interval.
    aa=(1-a)*h+a*bhi;bb=(1-b)*h+b*bhi
    with np.errstate(divide='ignore'):
        logs=np.log(np.maximum(np.maximum(aa,bb),0.)).sum(axis=-1)+np.log(np.array(P))[:,None]
    return np.exp(np.logaddexp.reduce(logs,axis=0))

def weighted_future_average(h,b,s):
    """Positive exact Bernstein integration, alpha integer, in exact arithmetic.
    One physical L: h,b arrays (components,n). No v quadrature is used.
    """
    n=h.shape[-1];alpha=math.ceil(math.sqrt(n));out=0.
    weights=np.full(n+1,alpha/(n+alpha))
    k=np.arange(n+1)
    for j in range(1,alpha):weights*= (k+j)/(n+j)
    for comp in range(3):
        A=(1-s)*h[comp]+s*b[comp];B=b[comp]
        coeff=np.array([1.])
        for degree in range(1,n+1):
            kk=np.arange(degree+1);new=np.zeros(degree+1)
            new[:-1]+= (degree-kk[:-1])/degree*A[degree-1]*coeff
            new[1:]+= kk[1:]/degree*B[degree-1]*coeff
            coeff=new
        out+=P[comp]*float(coeff@weights)
    return out
