from fractions import Fraction as F
from pathlib import Path
import json
v,b,c=F(1,128),F(3,5),F(1,4)
tau=b*v/c
checks=[]
for N in [9,10,13]:
 w={h:tau**(N-h) for h in range(1,N+1)};omega=sum(w.values());a=F(1,64)
 f={h:F(h*h+1,h+3) for h in range(N+1)}
 lhs=sum(w[h]*(c*f[h]-b*v*f[h-1]) for h in w);rhs=c*f[N]-b*v*tau**(N-1)*f[0]
 assert lhs==rhs
 T,q=F(7,11),F(3,7);E={h:F(2*h+5,h+9) for h in w}
 residual=sum(w[h]*(E[h]-q*T) for h in w)-(sum(w[h]*E[h] for h in w)-omega*q*T)
 assert residual==0
 A0,M=F(5,3),F(7,13)
 endpoint=b*v*tau**(N-1)*(A0*v*M)/(a*omega)
 target=A0*b**N*v**(N+1)*M/(c**(N-1)*a*omega)
 assert endpoint==target
 D={h:F(h+7,N+1) for h in w}
 pos=sum(w[h]*D[h]/(4*a) for h in w)/omega
 two=(D[N]+tau*D[N-1])/(8*a)
 assert tau<F(1,2) and omega<2 and pos>=two
 checks.append({'N':N,'telescope_residual':str(lhs-rhs),'center_residual':str(residual),'endpoint_residual':str(endpoint-target),'duplicate_reserve_slack':str(pos-two),'initial_order':N+1})
orders={'m9':9+10-9,'m8':8+10-8,'m7':8+10-7,'m6_qualitative':6+10-6,'early_max5':6+10-5,'initial':10+1}
assert orders=={'m9':10,'m8':10,'m7':11,'m6_qualitative':10,'early_max5':11,'initial':11}
print(json.dumps({'scope':'Finite rational checks of scalar identities, not geometry or an arbitrary-N theorem. No timestamps/PIDs invented.','v':str(v),'b':str(b),'c':str(c),'checks':checks,'N10_orders':orders},indent=2))
