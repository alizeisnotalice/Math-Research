"""Actual finite-field joint-vs-decoupled log-score MGF exploration."""
from pathlib import Path
from fractions import Fraction as Q
from decimal import Decimal as D,localcontext
import numpy as np,math,json,hashlib
import probe
H=Path(__file__).resolve().parent
def bound(n,K,b,aa,tt,V):
 U=2*(b+1);aa=np.asarray(aa);t=np.abs(tt);s1=float(np.dot(t,aa));s2=float(np.dot(t,aa*aa));delta=float(max(aa)*(U+V));z=(2+V)*s1/(1-(2+V)*max(aa))
 mass=U*V*(s1*s1+s2)*math.exp((U+V)*s1/(1-delta))/(1-delta)**2
 return n*math.exp(z)*math.sqrt(K+1)/K*mass

def main():
 target=H/'LOG_MGF.json';assert not target.exists();d=7;poly=0x83;K=127;b=3;B=8;g=2;states=[];s=1
 for _ in range(K):states.append(s);s=probe.mul(s,g,d,poly)
 assert s==1 and len(set(states))==K
 labels=np.array([sum(probe.trace(probe.mul(1<<i,x,d,poly),d,poly)<<i for i in range(b))for x in states]);ls=[Q(10,9),Q(7,6),Q(11,9)];p=Q(1,2);M=4;m=np.arange(1,M+1);bm=4*729/29*((np.sin(np.pi*m*10/27)/(np.pi*m))**4+3*(np.sin(np.pi*m/9)/(np.pi*m))**4);cm=np.sqrt(2*bm/3);V=float(cm.sum());aa=np.array([float(p/(K*L))for L in ls]);cases=[]
 for theta in [Q(1,7),Q(4,9)]:
  def end(v):
   off=v.numerator//v.denominator;r=v-off;u=(labels[(np.arange(K)+off)%K]+.5)/B;eps=1/128
   return np.clip((float(r)-u+eps/2)/eps,0,1)-float(r)
  hh=np.array([end(theta+K*L/2)-end(theta-K*L/2)for L in ls]);x=(np.arange(K)+float(theta))/K
  vv=np.array([np.sum((cm*np.sin(np.pi*m*float(L)))[:,None]*np.cos(2*np.pi*m[:,None]*x[None,:]),axis=0)for L in ls])
  for tt in [[1.,0.,0.],[.5,.5,0.],[1.,-1.,1.],[.1,.2,.3]]:
   joint=float(np.exp(np.array(tt)@np.log1p(aa[:,None]*(hh+vv))).mean());logarr=np.zeros((K,K))
   for j in range(3):logarr+=tt[j]*np.log1p(aa[j]*(hh[j,:,None]+vv[j,None,:]))
   dec=float(np.exp(logarr).mean());diff=abs(math.log(joint)-math.log(dec));cap=bound(1,K,b,aa,np.array(tt),V);assert diff<=cap+1e-12
   cases.append(dict(theta=str(theta),tau=tt,joint_MGF=joint,conditional_decoupled_MGF=dec,absolute_log_gap=diff,candidate_log_bound=cap,bound_passed=True))
 large=[]
 with localcontext()as c:
  c.prec=80;n=D(10)**32;kk=D(2)**47-1;pp=D(99)/100;Vd=D(40)/9;U=D(24)
  for J in [1,3,100,1000]:
   scales=[D(10)/9+(D(1)/9)*(D(i)/max(1,J-1))for i in range(J)];a=[pp/(kk*L)for L in scales];t=D(1)/J;s1=sum(a)/J;s2=sum(x*x for x in a)/J;delta=max(a)*(U+Vd);Z=(2+Vd)*s1/(1-(2+Vd)*max(a));cap=n*(kk+1).sqrt()/kk*U*Vd*(s1*s1+s2)*((U+Vd)*s1/(1-delta)+Z).exp()/(1-delta)**2
   large.append(dict(n=str(int(n)),J=J,tau_each=str(t),L1_tau='1',upper_bound_on_n_log_MGF_error=str(cap),A_norm_radius=str(delta),scope='Decimal evaluation of a proposed analytic error bound, not a simulation of n coordinates or a probability certificate.'))
 out=dict(small_orbit=dict(d=d,poly=poly,g=g,b=b,K=K,M=M,p=str(p),V_actual=V),actual_small_field_cases=cases,matched_parameter_bound_evaluation=large,script_SHA256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),scope='Small finite sums are numerical evidence; large rows evaluate bounds only. No C_n log gain.')
 target.write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
if __name__=='__main__':main()
