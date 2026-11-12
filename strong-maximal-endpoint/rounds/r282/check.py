from math import exp,log,log1p

def saddle(theta):
    def eq(h):
        return (1-theta)*exp(theta*h)-1-theta*exp(-(1-theta)*h)
    lo,hi=0.,1.
    while eq(hi)<0:
        hi*=2
    for _ in range(65):
        mid=(lo+hi)/2
        if eq(mid)<0: lo=mid
        else: hi=mid
    h=(lo+hi)/2
    mean=1/(1-theta)+h/(1+exp((1-theta)*h))
    return theta*(mean-h),mean,h

for b in (.1,.5,1.,2.):
    left,right=.00001,.99
    for _ in range(60):
        t=(left+right)/2
        if saddle(t)[0]<b:left=t
        else:right=t
    t=(left+right)/2
    beta,A,h=saddle(t)
    psi=log1p(exp(-(1-t)*h))-log1p(-t)-log1p(exp(-h))
    assert abs(beta-b)<1e-12 and abs(psi-t*h)<1e-12
    assert abs((b+psi)/t-A)<1e-12
    print(b,1/t,A)
