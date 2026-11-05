from fractions import Fraction as F
from math import comb
import json

def op(p):
 out=[F(0)]*(len(p)+1)
 for i,a in enumerate(p):
  b=a/F(i+1);out[0]+=b
  for j in range(i+2):out[j]+=b*comb(i+1,j)*(-1)**j
 while out and out[-1]==0:out.pop()
 return out
p=[F(1),F(-1)];rows=[]
for j in range(1,4):
 p=op(p);c=2*sum(a/F(i+1) for i,a in enumerate(p));rows.append({'j':j,'coefficients_ascending':[str(x) for x in p],'integral_coefficient':str(c)})
c3=c
for n in (2,3):
 B=F(n*(n-1),2);K=B*n*n/27;k=F(8*n,3)
 L=c3*B*k**6/(512*n**5);feedback=F(2*n,512)*L
 micro=F(2*n,64)*(B*k**3/(8*n*n))
 rows.append({'n':n,'feedback_av9':str(feedback),'feedback_matches':feedback==2*c3*K/27,'micro_av_nplus6':str(micro),'micro_matches':micro==2*K})
print(json.dumps(rows,indent=2))
