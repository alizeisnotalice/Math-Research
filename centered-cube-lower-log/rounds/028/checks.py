from pathlib import Path
from fractions import Fraction as F
from math import ceil,sqrt,sin,pi,log
import numpy as np,json
P=Path(__file__).resolve().parent
alpha=F(729,29);a=F(10,27);s=F(1,9)
b1lower=4*alpha*((7*a/11)**4+3*(7*s/11)**4);g0lower=b1lower*F(2,9)**2/6
threshold_root=ceil(F(380,21)/g0lower);K_sufficient=max(16,threshold_root**2)
CB=F(10000,243);AF=F(154,27);Jconstant=4*(F(81,100)*(2+2*CB)+F(6561,10000)*AF**2);assert Jconstant<360
b1=4*float(alpha)*((sin(pi*float(a))/pi)**4+3*(sin(pi*float(s))/pi)**4);g0=b1/6*sin(pi/9)**2
rows=[]
for M in [1,16,22,128,1024]:
 m=np.arange(1,M+1);b=4*float(alpha)*((np.sin(pi*m*float(a))/(pi*m))**4+3*(np.sin(pi*m*float(s))/(pi*m))**4)
 derivative_upper=pi*pi/3*np.sum(m*m*b);assert derivative_upper<float(CB)
 rows.append(dict(M=M,spectral_derivative_squared_upper=float(derivative_upper)))
char=[]
for eps,m in [(.001,1000000),(.001,100000000),(.0001,100000000)]:
 t=np.linspace(-1,1,20001);logphi=np.log1p(-2*np.sin(eps*t/2)**2);lhs=m*logphi;rhs=-m*eps**2*t*t/4
 assert np.max(lhs-rhs)<1e-12
 char.append(dict(eps=eps,m=m,max_log_envelope_violation=float(np.max(lhs-rhs)),integral=float(np.trapezoid(np.exp(lhs),t)),gaussian_integral_bound=2*sqrt(pi)/(eps*sqrt(m))))
scale=[]
for n in [10**8,10**16,10**24,10**40]:
 p=.5;target=16*log(n);d=1
 while 9*p*p*n/(29*((2**(d+1)-1)**2))>=target:d+=1
 K=2**d-1;v=9*p*p*n/(29*K*K);w=4;delta=1/(w*K*v);B=2**ceil(np.log2(v));M=ceil(sqrt(v));cx=40/9*((1-1/B)*sqrt(K+1)/K*np.sum(1/np.arange(1,M+1,dtype=float)**2)+pi*(1-1/B)**2/(2*K)*np.sum(1/np.arange(1,M+1,dtype=float)))
 scale.append(dict(n=str(n),K=K,v=v,B=B,M=M,sidon_sufficient_gate=(K+1>3*B**4),global_cross_gate_numeric=bool(cx<g0),cross_bound=cx,scaling_relative_error=abs(sqrt(K/delta)-3*p*sqrt(w/29)*sqrt(n))/(sqrt(K/delta))))
out=dict(b1=b1,g0=g0,b1_lower=str(b1lower),g0_lower=str(g0lower),uniform_cross_sufficient_K=K_sufficient,score_constant=str(Jconstant),score_constant_upper=360,spectral_checks=rows,characteristic_checks=char,scaling_checks=scale)
(P/'checks.json').write_text(json.dumps(out,indent=2));print(json.dumps(out,indent=2))
