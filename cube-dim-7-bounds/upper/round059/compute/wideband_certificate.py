from fractions import Fraction as Q
from pathlib import Path
import json,hashlib
P=Path(__file__).resolve().parent

def log_interval(x,N=48):
 y=(x-1)/(x+1); assert 0<=y<1
 s=2*sum((y**(2*k+1)/Q(2*k+1) for k in range(N)),Q(0))
 return s,s+2*y**(2*N+1)/(Q(2*N+1)*(1-y*y))
d=Q(391,1905);a=Q(283,2032);b=Q(63,64);qmin=1-b*d;qmax=1-a*d;R=qmin+d;vmin=a*d
rlo=R+Q(1,2**20);rhi=qmax+d-Q(1,2**20);chi=Q(147,10000);sig=Q(3,100);eret=Q(1,4096);ea=Q(1,65536);e=Q(1,2**20);P0=1-e-Q(1280,2**48)
Ulo=P0*(1-e)/(1+e);Uhi=(1+e)/(P0*(1-e));M=2*Uhi/qmin
rowb=(1-eret)*2*Ulo/qmax-2*e
loglo,loghi=log_interval(rhi/rlo);ml=(1-ea)*loglo;mh=(1-ea)*loghi
glo=(1-e)*ml*rowb/M;ghi=(1-e)*mh*rowb/M
pl,ph=log_interval(1/(1-d/(2*R)));mul=R/qmin*pl**2/(1-e);muh=R/qmin*ph**2/(1-e);cap=vmin/(2*qmax)
assert rowb>0 and rhi>rlo
verdict={'gamma_minus_3_97_gt_1_10':glo-Q(3,97)>Q(1,10),'mu_lt_chi':muh<chi,'chi_lt_cap':chi<cap}
vals={'d':d,'a':a,'b':b,'qmin':qmin,'qmax':qmax,'rlo':rlo,'rhi':rhi,'P0':P0,'Ulo':Ulo,'Uhi':Uhi,'M':M,'rowb':rowb,'gamma_lower':glo,'gamma_upper':ghi,'mu_lower':mul,'mu_upper':muh,'chi':chi,'caprate':cap,'gamma_margin_lower':glo-Q(3,97)-Q(1,10),'mu_margin_lower':chi-muh,'cap_margin':cap-chi}
out={'scope':'conditional scalar model only: no operator theorem or C_n bound','method':'ln x = 2 sum_(k=0..47) y^(2k+1)/(2k+1) + remainder; y=(x-1)/(x+1), 0 <= remainder <= 2*y^97/(97*(1-y^2))','verdict':verdict,'exact_rational':{k:str(v) for k,v in vals.items()},'decimal_diagnostics':{k:float(v) for k,v in vals.items()},'input_eps':'epsRet=1/4096 epsA=1/65536 epsJ=epsN=epsS=epsT=epsD=2^-20; P0=1-2^-20-1280*2^-48','code_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
path=P/'wideband_certificate.json';path.write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({'verdict':verdict,'decimal_diagnostics':{k:out['decimal_diagnostics'][k] for k in ('gamma_lower','gamma_margin_lower','mu_upper','chi','caprate','mu_margin_lower','cap_margin')},'code_sha256':out['code_sha256'],'certificate_sha256':hashlib.sha256(path.read_bytes()).hexdigest()},indent=2))
