from fractions import Fraction as Q
import json
N=128
def L(x):
 s=2*sum((x**(2*j+1)/Q(2*j+1)for j in range(N)),Q())
 t=2*x**(2*N+1)/(Q(2*N+1)*(1-x*x))
 return s,s+t
l2,u2=L(Q(1,3));lr,ur=L(Q(31,95));lo=5*l2+lr;hi=5*u2+ur
A=1+Q(65,63)**7;B=Q(64,63)**7
lower=A*(1+B*(1+7*lo));upper=A*(1+B*(1+7*hi));U=Q(77433340322253,10**12)
assert lower<upper<U<Q(77436179562403,10**12)
print(json.dumps(dict(status='PASS independent exact rational log63 enclosure',a=63,N=N,expression_lower=str(lower),expression_upper=str(upper),Cn_strict_upper=str(U),old_upper='77436179562403/1000000000000'),indent=2))
