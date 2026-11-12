from math import erfc, exp, log, pi, sqrt

mu = 1.2784645427610738
sigma = sqrt(mu)
def phi(z):
    return exp(-z*z/2)/sqrt(2*pi)
def Phi(z):
    return erfc(-z/sqrt(2))/2
def root(fn):
    l,r=-10.,10.
    assert fn(l)>0>fn(r)
    for _ in range(70):
        mid=(l+r)/2
        if fn(mid)>0:l=mid
        else:r=mid
    return (l+r)/2
for p,c in ((.5,1.),(.25,1.),(.75,2.)):
    a=sigma/c
    def g(z):
        return Phi(-z)+p*exp(-a*z+a*a/2)*Phi(z-a)
    zp=root(lambda z:g(z)-p)
    z=max(zp,root(lambda z:a*Phi(-z)-(1-p)*phi(z)))
    step=(z+12)/20000
    partial=sum(phi(-12+(i+.5)*step)*exp(a*(-12+(i+.5)*step-z))
                for i in range(20000))*step
    assert abs(g(z)-Phi(-z)-p*partial)<1e-7
    val=lambda t:sigma*t+c*log(min(1,g(t)/p))
    assert all(val(z)>=val(z+j/100)-1e-12 for j in range(-100,101))
    print(p,c,zp,z,val(z))
