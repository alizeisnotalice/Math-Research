from math import comb, log1p, expm1, exp
for n in (32,64,128,256,512):
    a=n//2
    z=3*n//4
    base=comb(z,a)
    logmiss=log1p(-1/base)
    union=private=0.0
    for j in range(a,n+1):
        prob=comb(n-a,j-a)/2**(n-a)
        m=comb(j,a)
        union+=prob*(-expm1(m*logmiss))/(m/base)
        private+=prob*exp((m-1)*logmiss)
    print(n,round(union,8),round(private,8))
