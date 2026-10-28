from pathlib import Path
from fractions import Fraction as F
import math,json
P=Path(__file__).resolve().parent
rows=[]
for w in [4,16,64,256,1024,10**6,10**12]:
 j=w.bit_length()-1;S1=F(j+1)+F(w,2**j);S2=F((j+1)*(j+2),2)+F(w*(j+3),2**j)
 finite1=sum((min(F(1),F(w,2**i)) for i in range(200)),F(0));finite2=sum(((i+1)*min(F(1),F(w,2**i)) for i in range(200)),F(0))
 assert 0<=S1-finite1<F(1,10**40) and 0<=S2-finite2<F(1,10**38)
 A=1+math.log(w);rows.append(dict(w=str(w),S1=str(S1),S2=str(S2),normalized_S1=float(S1)/A,normalized_S2=float(S2)/(A*A)))
growth=[]
for n in [10**8,10**16,10**24,10**40]:
 degree=1
 while 9*.25*n/(29*(2**(degree+1)-1)**2)>=16*math.log(n):degree+=1
 K=2**degree-1;eps=.5/K
 for power in [1,2]:
  A=1+power*math.log(n);growth.append(dict(n=str(n),K=K,power=power,epsilon_logw2=eps*A*A,epsilon_logw=eps*A))
(P/'series_checks.json').write_text(json.dumps(dict(series=rows,growth=growth),indent=2));print(json.dumps(dict(series=rows,growth=growth),indent=2))
