"""Actual scalar Fourier-minorant calibration, with no Gaussian replacement theorem.
Python3.14 has preinstalled SciPy; all target constants derive from accepted parameters.
"""
from pathlib import Path
from fractions import Fraction as F
import math,json,hashlib,sys
import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq,minimize_scalar
from scipy.special import ndtr
HERE=Path(__file__).resolve().parent


def parameters():
 K=2**47-1;q=K+1;n=10**32;p=F(99,100);L=F(7,6);U=F(24);V=F(40,9);W=U+V;a0=p/(K*L);v=9*p*p*n/(29*K*K);r=W*a0
 lam=math.sqrt(q)/K;vf=float(v);af=float(a0);rf=float(r);Wf=float(W);UV=float(U*V)
 beta=Wf*af/(math.sqrt(vf)*(1-rf));C= n*lam*UV*af*af
 k1=C/math.sqrt(vf)*((1+af*Wf)/(1-rf)**2+2/(1-rf));k2=C*(1+af*Wf)/(vf*(1-rf)**2)
 return dict(K=K,q=q,n=n,p=p,L=L,U=U,V=V,W=W,a0=a0,v=v,r=r,beta=beta,k1=k1,k2=k2)

def main():
 par=parameters();K=par['K'];q=par['q'];v=float(par['v']);beta=par['beta'];k1=par['k1'];k2=par['k2'];normal_const=1/math.sqrt(2*math.pi)
 def bound(t):
  if t==0:return 0.
  return min(2.,math.exp(min(100.,beta*t))*(k1*t+k2*t*t))
 tstar=brentq(lambda t:beta*t+math.log(k1*t+k2*t*t)-math.log(2),1e-5,1e6)
 def error(s,d):
  z=(s/d)**2;u0=s*tstar
  fun=lambda u:2*normal_const*math.exp(-u*u/2)*((1-z)+z*u*u)*bound(u/s)
  finite,err=quad(fun,0,u0,epsabs=1e-12,epsrel=1e-11,limit=150)
  tail=4*(ndtr(-u0)+z*u0*normal_const*math.exp(-u0*u0/2))
  return finite+tail,err,tail
 def gaussian(s,d,center):
  den=1+s*s;mass=s/math.sqrt(den)*math.exp(-center*center/(2*den))
  return mass*(1-(s*s/den+center*center*s**4/den**2)/(d*d))
 rows=[];checks=[]
 for width in [1,2,4]:
  d=width/(2*math.sqrt(v));center=(-2+width/2)/math.sqrt(v)
  def obj(z):
   s=z*d;return gaussian(s,d,center)-error(s,d)[0]
  grid=np.linspace(.025,1,100);scores=np.array([obj(z)for z in grid]);best=int(np.argmax(scores));lo=float(grid[max(0,best-1)]);hi=float(grid[min(len(grid)-1,best+1)])
  opt=minimize_scalar(lambda z:-obj(z),bounds=(lo,hi),method='bounded',options={'xatol':1e-12});candidates=[(float(opt.x),obj(opt.x)),(1.,obj(1.)),(float(grid[best]),float(scores[best]))];ratio,lower=max(candidates,key=lambda x:x[1]);s=ratio*d;err,quaderr,tail=error(s,d);G=gaussian(s,d,center)
  for z in [1/math.sqrt(3),ratio,1.]:
   ss=z*d;direct=ss*normal_const*quad(lambda u:(1-(ss*u/d)**2)*math.exp(-u*u/2-(center+ss*u)**2/2),-12,12,epsabs=1e-12)[0];assert abs(direct-gaussian(ss,d,center))<1e-11
   weightnorm=quad(lambda u:2*normal_const*math.exp(-u*u/2)*(1-z*z+z*z*u*u),0,12,epsabs=1e-12)[0];assert abs(weightnorm-1)<1e-11
  rows.append({'original_logscore_width':width,'normalized_halfwidth_d':d,'normalized_center':center,'h_definition':'v/2-2','s_over_d_numerical_optimizer':ratio,'Gaussian_E_f':G,'Fourier_error_integral_min2':err,'Gaussian_reference_lower_bridge':lower,'Gaussian_reference_is_positive':bool(lower>0),'quadrature_error_estimate':quaderr,'exact_min2_tail_evaluation_numeric':tail,'normal_interval_probability_reference':float(ndtr(center+d)-ndtr(center-d)),'scope':'Gaussian calibration only. Artificial-dec E f is not evaluated here; this is not a true-source probability bound.'})
  # Strict uniform nonpositivity certificate for every 0<s<=d in this Gaussian calibration.
  U,V,W,a0,r=map(par.get,['U','V','W','a0','r']);bracket=(1+a0*W)/(1-r)**2+2/(1-r)
  lamlo=F(math.isqrt(q),K);lamhi=F(math.isqrt(q)+1,K)
  clo=2*par['n']*lamlo*U*V*a0*a0*bracket/width;cui=2*par['n']*lamhi*U*V*a0*a0*bracket/width
  assert cui<=1
  floor_sqrt_v=math.isqrt(par['v'].numerator//par['v'].denominator);assert floor_sqrt_v**2<=par['v']<(floor_sqrt_v+1)**2
  dup=F(width,2*floor_sqrt_v);assert dup*dup<2
  error_lower=F(79,100)*F(6,7)*clo
  G_upper=2*dup*(1+dup*dup/2)/(3*F(17,10))
  assert G_upper<error_lower
  checks.append({'width':width,'kappa1_over_d_rational_lower':str(clo),'kappa1_over_d_rational_upper':str(cui),'uniform_error_integral_rational_lower':str(error_lower),'uniform_Gaussian_E_f_rational_upper':str(G_upper),'uniform_Gaussian_bridge_rational_upper':str(G_upper-error_lower),'uniform_bridge_upper_numeric':float(G_upper-error_lower),'all_s_gaussian_reference_strictly_negative_certified':True,'basis':'sqrt(2/pi)>79/100, 1-exp(-2)>6/7, sqrt3>17/10, sqrt(1+x)<=1+x/2; exact Fraction arithmetic.'})
 out={'parameters':{k:str(v)if isinstance(v,F)else v for k,v in par.items()},'tstar_min2_transition_numeric':tstar,'cases':rows,'uniform_Gaussian_nonpositivity_certificates':checks,'runtime_python':sys.version,'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'scope':'Strict Fourier bridge is proved separately. Gaussian comparison is calibration only; full-frequency error uses the exact min2 envelope, not bounded-T asymptotics.'}
 (HERE/'RESULTS.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({'parameters':out['parameters'],'cases':rows,'uniform_certified_gaussian_bridge_upper':[x['uniform_bridge_upper_numeric']for x in checks]},indent=2))
if __name__=='__main__':main()
