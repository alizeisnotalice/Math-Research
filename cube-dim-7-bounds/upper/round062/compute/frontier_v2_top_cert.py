from fractions import Fraction as F
from pathlib import Path
import json,hashlib,time
p=Path(__file__).resolve().parent;h=lambda b:hashlib.sha256(b).hexdigest();summary=p/'frontier_v2_summary.json';top=json.loads(summary.read_text())['top64'][0];d=F(top['d']);qhi=F(top['qhi']);J=top['J'];qlo=1-d;D=2**160;floor=lambda x:F((x.numerator*D)//x.denominator,D);ceil=lambda x:-floor(-x)
def bounds(x):
 y=(x-1)/(x+1);assert 0<=y<1;s=2*sum((y**(2*k+1)/F(2*k+1) for k in range(48)),F(0));return floor(s),ceil(s+2*y**97/(97*(1-y*y)))
rows=[];terms=F(0);mus=[];minslack=None
for j in range(J):
 l=qlo+(qhi-qlo)*F(j,J);u=qlo+(qhi-qlo)*F(j+1,J);R=l+d;ll,_=bounds((u+d)/R);_,phi=bounds(1/(1-d/(2*R)));mu=ceil(R/l*phi**2);eta=(1-u)/(2*u);v=floor(ll*l/u);terms+=v;mus.append(mu);slack=eta-mu;minslack=slack if minslack is None else min(minslack,slack);rows.append({'j':j,'qlo':str(l),'qhi':str(u),'mu_upper':str(mu),'eta':str(eta),'term_lower':str(v),'eta_margin':str(slack)})
c=terms-2*max(mus);out={'scope':'post-run rigorous arithmetic only for finite ideal model, not operator theorem or global optimality','input_summary_sha256':h(summary.read_bytes()),'source_sha256':h(Path(__file__).read_bytes()),'top_parameters':{k:top[k] for k in ('d','qhi','J','i','k')},'method':'48-term atanh log lower/upper with positive geometric remainder and dyadic outward rounding 2^-160','rows':rows,'c_model_lower':str(c),'min_eta_minus_mu_lower':str(minslack),'eligible_certified':minslack>0,'positive_certified':c>0,'executed_unix':time.time()};f=p/'frontier_v2_top_certificate.json';f.write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({'c_lower':float(c),'minimum_margin':float(minslack),'eligible':minslack>0,'positive':c>0,'certificate_sha256':h(f.read_bytes()),'source_sha256':out['source_sha256']}))
