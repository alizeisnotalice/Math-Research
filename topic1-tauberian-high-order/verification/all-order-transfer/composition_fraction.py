from fractions import Fraction as F
import json
from itertools import product

def mul(p,q,K):
 out=[F(0)]*(K+1)
 for i,x in enumerate(p):
  for j,y in enumerate(q[:K+1-i]):out[i+j]+=x*y
 return out

def cseq(q,K):
 c=[F(1)]
 for j in range(1,K+1):c.append(c[-1]*(j-q)/F(j+1))
 return c

def recurrence(p,c,K):
 a=[F(1)]
 for k in range(1,K+1):a.append(sum(((p+1)*i-k)*c[i]*a[k-i] for i in range(1,k+1))/k)
 return a

def independent_binomial(p,c,K):
 z=c.copy();z[0]=F(0);power=[F(1)]+[F(0)]*K;out=[F(0)]*(K+1);binom=F(1)
 for j in range(K+1):
  for i in range(K+1):out[i]+=binom*power[i]
  power=mul(power,z,K);binom*= (p-j)/F(j+1)
 return out
rows=[]
for n,r,q in product((2,3,4),(0,1,2,5,8,12),(F(0),F(1,2),F(1))):
 mismatches=[];count=0;K=r//n;c=cseq(q,K)
 for j in range(r+1):
  p=F(n+j,n);a=recurrence(p,c,K);direct=independent_binomial(p,c,K)
  for J in range(r+1):
   k,rem=divmod(J-j,n);expect=a[k] if J>=j and not rem else 0;value=direct[k] if J>=j and not rem else 0
   count+=1
   if expect!=value:mismatches.append([j,J])
 rows.append({'n':n,'r':r,'q':str(q),'assertions':count,'mismatches':mismatches,'pass':not mismatches})
print(json.dumps({'cases':len(rows),'basis_assertions':sum(x['assertions'] for x in rows),'all_pass':all(x['pass'] for x in rows),'rows':rows},indent=2))
