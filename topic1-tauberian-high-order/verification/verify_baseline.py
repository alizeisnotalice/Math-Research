#!/usr/bin/env python3
"""Exact formal-series verification; all coefficients are fractions, x=t^(1/n)."""
from fractions import Fraction as F
from math import comb
import json, time

def mul(a,b,N):
 c=[F(0)]*(N+1)
 for i,ai in enumerate(a):
  for j,bj in enumerate(b[:N+1-i]): c[i+j]+=ai*bj
 return c

def power_one(u,p,N):
 out=[F(1)]+[F(0)]*N; v=out[:]; b=F(1)
 for j in range(1,N+1):
  v=mul(v,u,N); b*= (p-j+1)/j
  out=[a+b*z for a,z in zip(out,v)]
 return out

def ds(n,N):
 a=[F(comb(n,k)*2**k) if k<=n else F(0) for k in range(N+1)]
 q=[F(1)]+[F(comb(n,k)*2**k) if k<n else F(0) for k in range(1,N+1)]
 d=[F(1)]
 for k in range(1,N+1): d.append(a[k]-sum(q[j]*d[k-j] for j in range(1,k+1)))
 assert mul(d,q,N)==a
 return d

def composite(n,m,N):
 # s_m(x^n)/x^n = sum_{j>=0} m*(-1)^j*binom(1/m,j+1)*x^(nj)
 u=[F(0)]*(N+1); b=F(1)
 for k in range(1,N//n+2):
  b*= (F(1,m)-(k-1))/k
  idx=n*(k-1)
  if idx<=N:
   coeff= -m*((-1)**k)*b
   if idx==0: assert coeff==1
   else: u[idx]=coeff
 root=power_one(u,F(1,n),N)
 y=[F(0)]+root[:-1]
 out=[F(0)]*(N+1); p=[F(1)]+[F(0)]*N
 for c in ds(n,N):
  out=[a+c*b for a,b in zip(out,p)]; p=mul(p,y,N)
 return out

def string(v): return [str(z) for z in v]
def main():
 start=time.time(); results={}
 for n in (1,2,3,10):
  N=max(2*n+1,n+3)
  results[str(n)]={'order':N,'d':string(ds(n,N)),'compositions':{str(m):string(composite(n,m,N)) for m in (1,2,3,10)}}
  d=ds(n,N)
  assert d[n]==2**n
  if n>=2: assert d[n+1]==-n*2**(n+1)
  else: assert d[1:]==[F(2)]+[F(0)]*(N-1)
  for m in (1,2,3,10):
   c=composite(n,m,N)
   assert c[2*n] == d[2*n]+F(2**n*(m-1),2*m)
 out={'variable':'x=t^(1/n), or x=epsilon^(1/n) for composition','n_values':[1,2,3,10],'m_values':[1,2,3,10],'series':results,'seconds':time.time()-start,'exact':True}
 print(json.dumps(out,indent=2))
if __name__=='__main__': main()
