from fractions import Fraction as F
from math import comb
import json

def add(p,q,N):return [p[i]+q[i] for i in range(N+1)]
def mul(p,q,N):
 out=[F(0)]*(N+1)
 for i,x in enumerate(p):
  for j,y in enumerate(q[:N+1-i]):out[i+j]+=x*y
 return out
def inv(p,N):
 out=[F(0)]*(N+1);out[0]=1/p[0]
 for k in range(1,N+1):out[k]=-sum(p[j]*out[k-j] for j in range(1,k+1))/p[0]
 return out
rows=[]
for n in range(2,7):
 N=n+20
 # v^n/(1+v)^n = v^n sum_j (-1)^j binom(n+j-1,j) v^j
 t=[F(0)]*(N+1)
 for j in range(N-n+1):t[n+j]=F((-1)**j*comb(n+j-1,j))
 direct=inv(add([F(1)]+[F(0)]*N,[-x for x in t],N),N)
 # ell coefficients: (1+v)^n - v^n times ell = (1+v)^n.
 denom=[F(comb(n,k)) if k<=n else F(0) for k in range(N+1)];denom[n]-=1
 rr=inv(denom,N); numer=[F(comb(n,k)) if k<=n else F(0) for k in range(N+1)];rec=mul(numer,rr,N)
 ell=[F(1)]
 for j in range(1,21):ell.append(-sum(F(comb(n,k))*ell[j-k] for k in range(1,min(j,n-1)+1)))
 for j in range(21):
  floor=sum((-1)**(j-(k-1)*n)*comb(n+j-1,j-(k-1)*n) for k in range(1,2+j//n))
  rows.append({'n':n,'j':j,'floor':floor,'fraction_series':str(direct[n+j]),'rational_recurrence':str(rec[n+j]),'ell_recurrence':str(ell[j]),'pass':direct[n+j]==rec[n+j]==ell[j]==floor})
print(json.dumps({'cases':len(rows),'all_pass':all(r['pass'] for r in rows),'rows':rows},indent=2))
