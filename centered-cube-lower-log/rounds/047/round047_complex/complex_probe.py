"""Finite-field complex moments and anchor-normalized characteristic functions.
No paper/GF artifact is modified. Only local finite sums; no probability lower bound.
"""
from pathlib import Path
from fractions import Fraction as F
from decimal import Decimal as D,localcontext
import numpy as np,math,json,hashlib
HERE=Path(__file__).resolve().parent

def mul(a,b,d,p):
 z=0
 while b:
  if b&1:z^=a
  b>>=1;a<<=1
  if a&(1<<d):a^=p
 return z

def trace(x,d,p):
 z=0
 for _ in range(d):z^=x;x=mul(x,x,d,p)
 assert z in (0,1)
 return z

def bound(K,b,aa,eta,V,anchor=0):
 U=2*(b+1);W=U+V;r=W*max(aa);s=np.abs(eta);s1=float(s@aa);s2=float(s@(aa*aa));Q=math.exp(W*s1/(1-r))
 # Mixed coefficients of (1+a0*(Uu+Vv))*(Q(u,v)-1).
 cap=math.sqrt(K+1)/K*U*V*Q*((1+aa[anchor]*W)*(s1*s1+s2)/(1-r)**2+2*aa[anchor]*s1/(1-r))
 return cap,dict(U=U,V=V,r=r,S1=s1,S2=s2)

def cpair(z):return [float(z.real),float(z.imag)]

def main():
 d=7;poly=0x83;K=(1<<d)-1;b=3;B=1<<b;states=[];x=1
 for j in range(K):states.append(x);x=mul(x,2,d,poly)
 assert x==1 and len(set(states))==K
 labels=np.array([sum(trace(mul(1<<i,x,d,poly),d,poly)<<i for i in range(b)) for x in states])
 scales=[F(10,9),F(7,6),F(11,9)];aa=np.array([float(F(1,2)/(K*L)) for L in scales]);M=4;m=np.arange(1,M+1)
 bm=4*729/29*((np.sin(np.pi*m*10/27)/(np.pi*m))**4+3*(np.sin(np.pi*m/9)/(np.pi*m))**4);cm=np.sqrt(2*bm/3);V=float(cm.sum());eps=F(1,128)
 logs=[];logs_dec=[]
 # Conditional theta-mixture exact finite enumeration, no theta quadrature claim.
 for theta in [F(1,7),F(4,9)]:
  def endpoint(v):
   off=v.numerator//v.denominator;r=v-off;u=(labels[(np.arange(K)+off)%K]+.5)/B
   return np.clip((float(r)-u+float(eps)/2)/float(eps),0,1)-float(r)
  H=np.array([endpoint(theta+K*L/2)-endpoint(theta-K*L/2) for L in scales]);xx=(np.arange(K)+float(theta))/K
  low=np.array([np.sum((cm*np.sin(np.pi*m*float(L)))[:,None]*np.cos(2*np.pi*m[:,None]*xx[None,:]),axis=0) for L in scales])
  logs.append(np.log1p(aa[:,None]*(H+low)))
  logs_dec.append(np.log1p(aa[:,None,None]*(H[:,:,None]+low[:,None,:])))
 def moment(tau,dec=False):
  arr=logs_dec if dec else logs
  return sum(np.exp(np.tensordot(tau,L,axes=(0,0))).mean() for L in arr)/len(arr)
 e0=np.array([1.,0.,0.]);A=moment(e0).real;At=moment(e0,True).real;assert abs(A-At)<2e-15
 rows=[];v=16
 for t in [[0.,0.,0.],[1.,0.,0.],[0.,1.,-1.],[3.,-2.,1.],[.1,.2,.3]]:
  eta=np.array(t)/math.sqrt(v);tau=e0+1j*eta;mu=moment(tau);mt=moment(tau,True);cap,par=bound(K,b,aa,eta,V);f=mu/A;ft=mt/At
  assert abs(mu-mt)<=cap+2e-14
  for n in [1,100,10000]:
   delta=abs(f**n-ft**n);ncap=min(2.,n*cap/min(A,At));assert delta<=ncap+5e-12
   rows.append(dict(t=t,v=v,n=n,mu=cpair(mu),decoupled_mu=cpair(mt),anchor=A,decoupled_anchor=At,one_coordinate_moment_gap=abs(mu-mt),anchor_aware_moment_bound=cap,characteristic_gap=delta,characteristic_bound=ncap,bound_parameters=par,passed=True))
 # The log-corollary comparison uses principal Log only inside a disk centered at 1.
 complex_rows=[]
 for tau in [np.array([.3+.7j,-.2+.1j,.4-.2j]),e0+1j*np.array([1.,-1.,0.])/4]:
  mu=moment(tau);mt=moment(tau,True);U=2*(b+1);W=U+V;r=W*max(aa);s=np.abs(tau);s1=float(s@aa);s2=float(s@(aa*aa));Z=(2+V)*s1/(1-(2+V)*max(aa));rho=math.expm1(Z)
  E=math.sqrt(K+1)/K*U*V*(s1*s1+s2)/(1-r)**2*math.exp(W*s1/(1-r));assert rho<1
  gap=abs(np.log(mu)-np.log(mt));cap=E/(1-rho);assert gap<=cap+1e-14
  complex_rows.append(dict(tau=[cpair(t) for t in tau],rho=rho,moment_gap=abs(mu-mt),moment_bound=E,principal_log_gap=gap,local_log_bound=cap,passed=True))
 large=[]
 with localcontext() as c:
  c.prec=80;n=D(10)**32;Kbig=D(2)**47-1;U=D(24);Vd=D(40)/9;W=U+Vd;p=D(99)/100;vbig=D(9)*(D(99)/100)**2*n/(D(29)*Kbig*Kbig);Ls=[D(10)/9,D(7)/6,D(11)/9];as_=[p/(Kbig*L) for L in Ls];r=W*max(as_)
  for t in [[1,0,0],[0,1,-1],[3,-2,1]]:
   eta=[D(abs(z))/vbig.sqrt() for z in t];s1=sum(z*a for z,a in zip(eta,as_));s2=sum(z*a*a for z,a in zip(eta,as_));cap=n*(Kbig+1).sqrt()/Kbig*U*Vd*(W*s1/(1-r)).exp()*((1+as_[0]*W)*(s1*s1+s2)/(1-r)**2+2*as_[0]*s1/(1-r))
   large.append(dict(n=str(int(n)),v=str(vbig),t=t,A_Haar='1 exactly',upper_bound_n_characteristic_gap=str(min(D(2),cap)),scope='Decimal evaluation only; v=9*p^2*n/(29*K^2), the accepted baseline normalization; not a probability certificate.'))
 out=dict(small_field=dict(d=d,modulus=poly,generator=2,b=b,K=K,M=M,phase_mixture=['1/7','4/9'],epsilon=str(eps),v=v),anchor_normalized_cases=rows,complex_log_cases=complex_rows,large_bound_evaluations=large,script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),scope='Finite sums are floating numerical evidence. Bounds are proved separately. No C_n log gain or event lower bound.')
 (HERE/'RESULTS.json').write_text(json.dumps(out,indent=2)+'\n')
 print(json.dumps({'cases':len(rows),'max_moment_gap':max(r['one_coordinate_moment_gap'] for r in rows),'max_CF_gap':max(r['characteristic_gap'] for r in rows),'complex_log_cases':complex_rows,'large':large},indent=2))
if __name__=='__main__':main()
