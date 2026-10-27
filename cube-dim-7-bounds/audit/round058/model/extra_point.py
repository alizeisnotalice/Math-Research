from fractions import Fraction as F
from pathlib import Path
import json,time,os,resource,hashlib
st=time.monotonic();cpu=time.process_time();here=Path(__file__).resolve().parent;rr=(here/'extra_registration.json').read_bytes();r=json.loads(rr);i,j,k=r['state'];assert (i,j,k)==(127,127,0)
def logint(x):
 z=(x-1)/(x+1);assert 0<=z<1;p=z;s=F(0)
 for k in range(48):s+=2*p/(2*k+1);p*=z*z
 return s,s+2*p/(97*(1-z*z))
d=F(1,5)+F(2*i,15*127);a=F(1,16)+F(3*j,8*127);b=F(9,16)+F(27*k,64*127);q=1-b*d;Q=1-a*d;R=q+d;l,h=logint(1/(1-d/(2*R)));ml=R/q*l*l;mu=R/q*h*h;eta=a*d/(2*Q);u,v=logint((Q+d)/R);cl=u*q/Q-2*mu;cu=v*q/Q-2*ml
assert mu<eta and cu<0
out={'state':[i,j,k],'mu_less_eta':True,'c_negative':True,'eta_minus_mu_upper':str(eta-mu),'c_upper':str(cu),'registration_sha256':hashlib.sha256(rr).hexdigest(),'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'pid':os.getpid(),'cpu_seconds':time.process_time()-cpu,'wall_seconds':time.monotonic()-st,'peak_rss_bytes':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
(here/'extra_point.json').write_text(json.dumps(out,indent=2)+'\n');print('mu<eta, c<0',float(eta-mu),float(cu))
