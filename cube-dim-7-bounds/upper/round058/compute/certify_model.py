from fractions import Fraction as F
import json
from pathlib import Path
P=Path(__file__).resolve().parent
def log_interval(x,N=48):
 y=(x-1)/(x+1);assert 0<=y<1
 s=sum((y**(2*k+1)/F(2*k+1) for k in range(N)),F(0))*2
 return s,s+2*y**(2*N+1)/(F(2*N+1)*(1-y*y))
top=json.loads((P/'summary.json').read_text())['top64'][0];d,a,b=(F(top[x]) for x in ('d','a','b'));qmin=1-b*d;qmax=1-a*d;R=qmin+d;phi0,phi1=log_interval(1/(1-d/(2*R)));mu0=R/qmin*phi0**2;mu1=R/qmin*phi1**2;eta=a*d/(2*qmax);ln0,ln1=log_interval((qmax+d)/R);c0=ln0*qmin/qmax-2*mu1;c1=ln1*qmin/qmax-2*mu0;assert 0<mu0<=mu1<eta and c0>0
out={'scope':'exact rational bounds only for conditional ideal numerical model, not an operator theorem','input':{k:top[k] for k in ('d','a','b','i','j','k')},'terms':48,'mu_lower':str(mu0),'mu_upper':str(mu1),'eta':str(eta),'c_lower':str(c0),'c_upper':str(c1),'verified_strict':['mu_upper < eta','c_lower > 0'],'decimal_display':{'c_lower':float(c0),'c_upper':float(c1),'eta_minus_mu_upper':float(eta-mu1)}};(P/'model_certificate.json').write_text(json.dumps(out,indent=2)+'\n');print(out['decimal_display'])
