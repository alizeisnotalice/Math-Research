from fractions import Fraction as F
from pathlib import Path
import json,hashlib,time,os,resource,sys
sys.set_int_max_str_digits(0);st=time.monotonic();cpu=time.process_time();here=Path(__file__).resolve().parent;src=Path('/Users/zhengzhihao/Desktop/T/cube-r030-compute/cube-dim-7-bounds/upper/round067/compute');rr=(here/'registration.json').read_bytes();reg=json.loads(rr);sha=lambda b:hashlib.sha256(b).hexdigest()
for f,h in reg['source_sha256'].items():assert sha((src/f).read_bytes())==h
old=json.loads((src/'candidate_certificate.json').read_bytes());assert old['code_sha256']==reg['source_sha256']['candidate.py']
def logint(x):
 z=(x-1)/(x+1);assert 0<=z<1;p=z;total=F(0)
 for k in range(40):total+=2*p/(2*k+1);p*=z*z
 return total,total+2*p/(81*(1-z*z))
s=F(3,2);ql=F(2,3);qh=F(917,1000);sigma=F(83,1000);chi=F(9,200);trim=F(1,2**30);e=F(1,2**20);ea=F(1,65536);er=F(1,4096);h=F(1,2**80);P0=1-e-F(1,2**60);assert 0<h<e and 0<P0<1
rl=s*ql+trim;rh=s*qh-trim;assert rl==1+trim and rh==F(2751,2000)-trim and rl<rh
ulo=P0*(1-e)/(1+e);uhi=(1+e)/(P0*(1-e));M=2*(1+e)*uhi;b=2*(1-er)*(1-e)*ulo-2*e;m=F(637,2000);gamma=(1-e)*m*b/M;assert b>0
loglo,loghi=logint(F(6,5));masslo,masshi=logint(rh/rl);mu0up=F(1,2)-F(5,2)*loglo;bridge=F(2211,50000)/(1-e);qbar=F(18233,100000);variance=qbar**3/(1-e);eta=(1-qh)/(2*qh);tau=sigma/2-chi*(1-sigma);raw=F(317,1000)-F(83,917);pingap=chi-F(223,5000)-F(1,5000)
checks={'mass_gt_m':(1-ea)*masslo>m,'gamma_gt_clean':gamma>F(317,1000),'mu0_lt_221_5000':mu0up<F(221,5000),'bridge_lt_223_5000':bridge<F(223,5000),'pin_gap_identity':pingap==F(1,5000),'eta_identity':eta-chi==F(47,183400),'tau_identity':tau==F(47,200000),'raw_identity':raw==F(207689,917000) and raw-F(9,40)==F(341,229250),'variance_lt_1_150':variance<F(1,150),'qbar_above_log_upper':qbar>loghi}
assert all(checks.values()) and all(old['verdict'][k] for k in old['verdict'])
E=old['exact'];vals={'m':m,'gamma':gamma,'gamma_margin':gamma-F(317,1000),'bridge_upper':bridge,'pin_gap':pingap,'eta_minus_chi':eta-chi,'tau':tau,'raw':raw,'raw_minus_9_40':raw-F(9,40),'variance_upper':variance}
assert all(F(E[k])==v for k,v in vals.items());assert (1-ea)*masslo<=F(E['logmass_lower'])<=(1-ea)*masshi and F(1,2)-F(5,2)*loghi<=F(E['mu0_upper'])<=mu0up
vals.update({'logmass_lower_40':(1-ea)*masslo,'mass_margin':(1-ea)*masslo-m,'mu0_upper_40':mu0up,'mu0_margin':F(221,5000)-mu0up,'profile_margin':F(223,5000)-bridge,'variance_margin':F(1,150)-variance})
out={'scope':'conditional scalar only; target profile, geometry and operator theorem not audited','checks':checks,'exact':{k:str(v) for k,v in vals.items()},'approx':{k:float(v) for k,v in vals.items()},'method':'40-term atanh series with rational geometric tail; mu0 upper uses lower log bound','registration_sha256':sha(rr),'script_sha256':sha(Path(__file__).read_bytes()),'pid':os.getpid(),'cpu_seconds':time.process_time()-cpu,'wall_seconds':time.monotonic()-st,'peak_rss_bytes':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
(here/'audit.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({'checks':checks,'approx':out['approx'],'pid':out['pid'],'cpu':out['cpu_seconds'],'wall':out['wall_seconds']},indent=2))
