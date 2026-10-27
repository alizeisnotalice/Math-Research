from fractions import Fraction as Q
from pathlib import Path
import json,hashlib,time
start=time.perf_counter();N=128
s=2*sum((Q(1,(2*j+1)*3**(2*j+1))for j in range(N)),Q())
tail=Q(9,4*(2*N+1)*3**(2*N+1))
A=1+Q(33,32)**7;B=Q(65,64)**7
lo=A*(1+B*(1+42*s));hi=A*(1+B*(1+42*(s+tail)))
U=Q(77436179562403,10**12)
assert lo<hi<U<81<128
assert 1+Q(7,10)+Q(7,10)**2/2+Q(7,10)**3/6==Q(12013,6000)>2
r=A*(1+B*(1+42*Q(7,10)))
assert hi<r<Q(10013,125)<81
print(json.dumps(dict(status='PASS exact rational scalar enclosure; selection proof separately audited',N=N,log2_lower=str(s),log2_upper=str(s+tail),expression_lower=str(lo),expression_upper=str(hi),strict_C7_rational_upper_if_theorem_accepted=str(U),old_upper='128',exact_comparison=True,elapsed_seconds=time.perf_counter()-start),indent=2))
