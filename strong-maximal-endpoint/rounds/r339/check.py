from math import comb,log,exp,log1p,expm1
a,l,x=.5,.5,.625
h=x*log(x)-a*log(a)-(x-a)*log(x-a)
rate=x*log(x/l)+(1-x)*log((1-x)/(1-l))
coeff=.5+log(x*(1-l)/(l*(1-x)))/(2*log(x/(x-a)))
for n in (32,64,128,256,512):
    theta=exp(-n*h)
    mean=sum((comb(n,j)/2**n)*(-expm1(comb(j,n//2)*log1p(-theta))) for j in range(n//2,n+1))
    print(n,round(mean,9),round(log(mean)+n*rate+coeff*log(n),5))
