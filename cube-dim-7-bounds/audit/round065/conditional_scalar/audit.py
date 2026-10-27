from fractions import Fraction as F
from pathlib import Path
import json,hashlib,time,os,resource,sys
sys.set_int_max_str_digits(0);st=time.monotonic();cpu=time.process_time();here=Path(__file__).resolve().parent;src=Path('/Users/zhengzhihao/Desktop/T/cube-r030-compute/cube-dim-7-bounds/upper/round065/compute');rr=(here/'registration.json').read_bytes();reg=json.loads(rr);sha=lambda b:hashlib.sha256(b).hexdigest()
for f,h in reg['source_sha256'].items():assert sha((src/f).read_bytes())==h
producer=json.loads((src/'candidate_certificate.json').read_bytes());assert producer['code_sha256']==reg['source_sha256']['candidate.py']
# log x = 2 sum y^(2k+1)/(2k+1) + tail <= 2 y^(2N+1)/((2N+1)(1-y²)).
def bounds(x):
 y=(x-1)/(x+1);assert 0<=y<1;term=y;value=F(0)
 for k in range(40):value+=2*term/(2*k+1);term*=y*y
 return value,value+2*term/(81*(1-y*y))
e=F(1,1024);s=F(3,2);ql=F(2,3);qh=F(909,1000);trim=F(1,2**30);rlo=1+trim;rhi=F(2727,2000)-trim;assert ql<qh and s*ql+trim==rlo and s*qh-trim==rhi and rlo<rhi
P0=1-e-F(1,2**50);ul=P0*(1-e)/(1+e);uh=(1+e)/(P0*(1-e));M=2*(1+e)*uh;mass=F(3097,10000);b=2*(1-e)**2*ul-2*e;assert b>0
lrlo,lrhi=bounds(rhi/rlo);lnlo,lnhi=bounds(F(6,5));qbar=F(18233,100000);gamma=(1-e)*mass*b/M;mu=s*qbar**2/(1-e);variance=2*s*qbar**3/(3*(1-e));eta=(1-qh)/(2*qh);chi=F(1,20);sigma=F(91,1000);tau=sigma/2-chi*(1-sigma);raw=F(61,200)-F(91,909)
checks={'mass_gt_m':(1-e)*lrlo>mass,'gamma_gt_61_200':gamma>F(61,200),'Qbar_gt_log_6_5':qbar>lnhi,'mu_lt_999_20000_lt_chi':mu<F(999,20000)<chi,'variance_lt_1_150':variance<F(1,150),'eta_minus_chi_identity':eta-chi==F(1,18180),'tau_identity':tau==F(1,20000),'raw_identity':raw==F(37249,181800) and raw-F(1,5)==F(889,181800)}
assert all(checks.values()) and all(producer['verdict'][k] for k in checks)
E=producer['exact'];eq={'m':mass,'gamma':gamma,'gamma_margin':gamma-F(61,200),'Qbar':qbar,'mu_upper':mu,'mu_margin':F(999,20000)-mu,'var_upper':variance,'eta_minus_chi':eta-chi,'tau_gap':tau,'raw':raw,'raw_minus_1_5':raw-F(1,5)}
assert all(F(E[k])==v for k,v in eq.items());assert (1-e)*lrlo<=F(E['mass_lower'])<=(1-e)*lrhi and lnlo<=F(E['log_6_5_upper'])<=lnhi
vals={'mass_lower_40':(1-e)*lrlo,'mass_margin':(1-e)*lrlo-mass,'gamma':gamma,'gamma_margin':gamma-F(61,200),'Qbar':qbar,'log_6_5_upper_40':lnhi,'mu_upper':mu,'mu_margin':F(999,20000)-mu,'variance_upper':variance,'variance_margin':F(1,150)-variance,'eta_minus_chi':eta-chi,'tau_gap':tau,'raw':raw,'raw_minus_1_5':raw-F(1,5)}
out={'scope':'conditional scalar only; geometry, operator theorem and continuous optimization not audited','checks':checks,'exact':{k:str(v) for k,v in vals.items()},'approx':{k:float(v) for k,v in vals.items()},'method':'40-term atanh positive series and geometric tail','registration_sha256':sha(rr),'script_sha256':sha(Path(__file__).read_bytes()),'pid':os.getpid(),'cpu_seconds':time.process_time()-cpu,'wall_seconds':time.monotonic()-st,'peak_rss_bytes':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
(here/'audit.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({'checks':checks,'approx':out['approx'],'pid':out['pid'],'cpu':out['cpu_seconds'],'wall':out['wall_seconds']},indent=2))
