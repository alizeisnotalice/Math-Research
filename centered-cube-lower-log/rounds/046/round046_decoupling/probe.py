"""Actual small-field Walsh/Fourier exploration before proposing decoupling bounds."""
from pathlib import Path
from fractions import Fraction as Q
import numpy as np,json,math,hashlib,time
H=Path(__file__).resolve().parent
def mul(a,b,d,poly):
 z=0
 while b:
  if b&1:z^=a
  b>>=1;a<<=1
  if a>>d:a^=poly
 return z

def power(a,n,d,p):
 z=1
 while n:
  if n&1:z=mul(z,a,d,p)
  a=mul(a,a,d,p);n>>=1
 return z

def trace(a,d,p):
 z=0
 for _ in range(d):z^=a;a=mul(a,a,d,p)
 assert z in [0,1];return z

def fwht(a):
 a=np.array(a,dtype=float,copy=True);n=a.shape[-1];h=1
 while h<n:
  sh=a.reshape(*a.shape[:-1],-1,2,h);u=sh[...,0,:].copy();v=sh[...,1,:].copy();sh[...,0,:]=u+v;sh[...,1,:]=u-v;h*=2
 return a/n

def main():
 out=H/'PROBE.json';assert not out.exists();t=time.time();rows=[]
 for b in range(1,12):
  B=2**b;a=np.arange(B);v=(a[None,:]>=a[:,None]).astype(float);c=fwht(v);l1=abs(c).sum(1)
  rows.append(dict(b=b,B=B,max_full_Walsh_l1=float(l1.max()),max_nonconstant_l1=float((abs(c).sum(1)-abs(c[:,0])).max()),dyadic_b_bound_passed=bool(np.all(l1<=b+1e-12))))
 # Source-style two-endpoint score on an actual finite-field orbit.
 d=10;p=(1<<10)|(1<<3)|1;g=2;K=2**d-1;b=4;B=2**b;states=[];x=1
 for j in range(K):states.append(x);x=mul(x,g,d,p)
 assert len(set(states))==K and x==1
 labels=np.array([sum(trace(mul(1<<i,x,d,p),d,p)<<i for i in range(b))for x in states]);theta=Q(1,7);L=Q(10,9);eps=Q(1,128)
 def endpoint(s):
  off=s.numerator//s.denominator;r=s-off;u=(labels[(np.arange(K)+off)%K]+.5)/B
  return np.clip((float(r)-u+float(eps)/2)/float(eps),0,1)-float(r)
 score=endpoint(theta+K*L/2)-endpoint(theta-K*L/2)
 xx=(np.arange(K)+float(theta))/K;fourier=np.cos(2*np.pi*xx)+.3*np.sin(6*np.pi*xx)
 # Joint trace source may have relations, but empirical covariance need not vanish.
 records=[]
 for z in [.001,.01,.1,1.]:
  F=np.exp(z*score);V=np.exp(z*fourier);cross=float(np.mean(F*V)-np.mean(F)*np.mean(V))
  # Conservative suggested A-norm bound, to be proved after this probe.
  bound=math.sqrt(K+1)/K*math.expm1(2*(b+1)*z)*math.expm1(1.3*z)
  records.append(dict(t=z,mean_product=float(np.mean(F*V)),product_means=float(np.mean(F)*np.mean(V)),cross=cross,candidate_bound=bound,bound_passed=abs(cross)<=bound+1e-12))
 # Test all nonzero trace characters and several multiplicative frequencies.
 gauss=[]
 for a in [1,2,3,7,31,1023]:
  ch=np.array([1-2*trace(mul(a,x,d,p),d,p)for x in states])
  for m in [1,2,3,7,31,127]:
   gg=abs(np.sum(ch*np.exp(2j*np.pi*m*np.arange(K)/K)));gauss.append(float(abs(gg-math.sqrt(K+1))))
 result=dict(threshold_Walsh_scan=rows,source_probe=dict(d=d,poly=p,g=g,b=b,B=B,theta=str(theta),L=str(L),epsilon=str(eps),score_mean=float(score.mean()),score_Fourier_covariance=float(np.mean(score*fourier)-score.mean()*fourier.mean()),exponential_cases=records),gauss_checks=len(gauss),max_Gauss_modulus_error=max(gauss),seconds=time.time()-t,numpy=np.__version__,script_SHA256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),scope='Actual finite-field orbit and floating transforms; not a C_n bound or interval certificate.')
 assert all(r['dyadic_b_bound_passed']for r in rows) and all(r['bound_passed']for r in records)
 out.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
if __name__=='__main__':main()
