"""Independent formal power series via rational convolutions; no CAS shared polynomial."""
from fractions import Fraction as Q
from pathlib import Path
import json, math
ROOT=Path(__file__).resolve().parents[1]/'results';N=16

def add(x,y):return [x[i]+y[i] for i in range(N)]
def scale(x,c):return [v*c for v in x]
def mul(x,y):return [sum((x[j]*y[i-j] for j in range(i+1)),Q(0)) for i in range(N)]
def const(c):return [Q(c)]+[Q(0)]*(N-1)
def root1(x):
 out=const(1)
 for i in range(1,N):out[i]=(x[i]-sum(out[j]*out[i-j] for j in range(1,i)))/2
 return out

def divide(x,y):
 z=const(0)
 for i in range(N):z[i]=(x[i]-sum(y[j]*z[i-j] for j in range(1,i+1)))/y[0]
 return z

def check():
 tables=json.loads((ROOT/'direction_common_scale_wolfram.json').read_text())['coefficient_tables'];rows=[]
 for row in tables:
  b=Q(row['b']);a=1-b;al=const(1);al[1]=Q(-1);r=const(0);r[1]=Q(1);r0=r.copy()
  for n in range(2,15):
   r2=mul(r,r);e=add(scale(r,2*b),root1(add(const(1),scale(r2,4*(a-b)**2))))
   F=add(add(const(1),mul(e,e)),add(scale(mul(add(al,r),e),-2),scale(r2,4*a*b)))
   r[n]=F[n]/2
   P0=add(add(scale(mul(r0,r0),4*a*b),scale(mul(al,r0),-2)),add(const(1),scale(mul(al,al),-1)))
   r0[n]=P0[n]/2
  got=list(map(Q,row['radius_coefficients']));assert r[1:15]==got
  # shift each series by one, so division has nonzero constant.
  rat=divide(r[1:]+[Q(0)],r0[1:]+[Q(0)]);sig=[-rat[j+1] for j in range(11)]
  assert sig==list(map(Q,row['sigma_coefficients'])),(b,sig[:3],row['sigma_coefficients'][:3])
  d=1-2*b;K=a*d;r2=mul(r,r);A=add(const(1),scale(r2,2*K));B=add(al,scale(r,d))
  polynomial=add(mul(A,A),add(scale(mul(mul(r,A),B),-4*b),scale(mul(add(const(1),scale(r2,4*(d*d-b*b))),mul(B,B)),-1)))
  assert all(x==0 for x in polynomial[:15])
  rows.append(dict(b=str(b),radius_coefficients=[str(x) for x in r[1:15]],sigma_coefficients=[str(x) for x in sig],polynomial_residual_through_order=14))
 return {'directions':len(rows),'radius_orders':14,'budget_orders':10,'exact_polynomial_residual_order':14,'method':'independent Fraction sqrt-series convolution'}
if __name__=='__main__':print(check())
