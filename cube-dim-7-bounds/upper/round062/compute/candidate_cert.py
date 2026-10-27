from fractions import Fraction as F
from pathlib import Path
import json,hashlib
P=Path(__file__).resolve().parent
D=2**200
lo=lambda x:F((x.numerator*D)//x.denominator,D)
up=lambda x:-lo(-x)
def ln_bounds(x,N=48):
 y=(x-1)/(x+1);assert 0<=y<1
 s=2*sum((y**(2*k+1)/F(2*k+1) for k in range(N)),F(0));return lo(s),up(s+2*y**(2*N+1)/(97*(1-y*y)))
d=F(7,25);J=128;q0=F(18,25);q1=F(24,25);gap=F(1,2**30);e=F(1,2**20);eret=F(1,4096);ea=F(1,65536);P0=1-e-F(1280,2**48);ul=P0*(1-e)/(1+e);uh=(1+e)/(P0*(1-e));sigma=F(3,50);total=F(0);rows=[]
for j in range(J):
 l=q0+(q1-q0)*F(j,J);u=q0+(q1-q0)*F(j+1,J);rl=l+d+gap;rh=u+d-gap; ll,lh=ln_bounds(rh/rl);pl,ph=ln_bounds(1/(1-d/(2*(l+d))));mu=up((l+d)/l*ph**2/(1-e));ceilmu=(mu.numerator*10**6+mu.denominator-1)//mu.denominator;chi=F(ceilmu,10**6)+F(1,10000);eta=(1-u)/(2*u);margin=sigma/2-chi*(1-sigma);b=(1-eret)*2*ul/u-2*e;M=2*uh/l;term=lo((1-e)*(1-ea)*ll*b/M);total+=term
 vals={'j':j,'qlo':l,'qhi':u,'rlo':rl,'rhi':rh,'log_lower':ll,'log_upper':lh,'phi_upper':ph,'mu_upper':mu,'chi':chi,'eta':eta,'eta_margin':eta-chi,'sigma_margin':margin,'rowb':b,'gamma_lower':term};rows.append({k:str(v) if k!='j' else v for k,v in vals.items()})
out={'scope':'conditional scalar model; no operator theorem','method':'48 term atanh expansion remainder <=2*y^97/(97*(1-y*y)); all transcendent bounds and final positive terms rounded outward to dyadic denominator 2^200','constant':{'d':str(d),'J':J,'sigma':str(sigma),'Ulo':str(ul),'Uhi':str(uh)},'rows':rows,'gamma_total_lower':str(total),'gamma_minus_3_47_minus_3_20_lower':str(total-F(3,47)-F(3,20)),'verdict':{'all_eta_gt_chi':all(F(x['eta_margin'])>0 for x in rows),'all_sigma_margin_positive':all(F(x['sigma_margin'])>0 for x in rows),'gamma_minus_3_47_gt_3_20':total-F(3,47)>F(3,20)},'code_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()};f=P/'candidate_certificate.json';f.write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({'verdict':out['verdict'],'gamma_lower':float(total),'gamma_slack':float(total-F(3,47)-F(3,20)),'min_eta_slack':float(min(F(x['eta_margin']) for x in rows)),'min_sigma_slack':float(min(F(x['sigma_margin']) for x in rows)),'code_hash':out['code_sha256'],'certificate_hash':hashlib.sha256(f.read_bytes()).hexdigest()}))
