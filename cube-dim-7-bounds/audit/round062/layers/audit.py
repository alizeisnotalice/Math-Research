from fractions import Fraction as F
from pathlib import Path
import json,time,hashlib,resource,os
st=time.monotonic();cpu=time.process_time();here=Path(__file__).resolve().parent;src=Path('/Users/zhengzhihao/Desktop/T/cube-r030-compute/cube-dim-7-bounds/upper/round062/compute');reg=json.loads((here/'registration.json').read_bytes());sha=lambda b:hashlib.sha256(b).hexdigest()
for name,h in reg['sha256'].items():assert sha((src/name).read_bytes())==h
fixed=json.loads((src/'candidate_certificate.json').read_bytes());gate=json.loads((src/'candidate_gate.json').read_bytes());top=json.loads((src/'frontier_v2_top_certificate.json').read_bytes());summary=json.loads((src/'frontier_v2_summary.json').read_bytes())
assert fixed['code_sha256']==sha((src/'candidate_cert.py').read_bytes()) and gate['input_sha256']==sha((src/'candidate_certificate.json').read_bytes()) and gate['source_sha256']==sha((src/'candidate_gate.py').read_bytes())
assert top['source_sha256']==sha((src/'frontier_v2_top_cert.py').read_bytes()) and top['input_summary_sha256']==sha((src/'frontier_v2_summary.json').read_bytes())
def floor(x,D):return F(x.numerator*D//x.denominator,D)
def ceil(x,D):return -floor(-x,D)
for D in [2**160,2**200]:
 for x in [F(-1,3),-F(1,2*D),F(0),F(1,2*D),F(1,3)]:assert floor(x,D)<=x<=ceil(x,D) and ceil(x,D)-floor(x,D)<=F(1,D)
def loginterval(x,N):
 y=(x-1)/(x+1);assert 0<=y<1; z=y;s=F(0)
 for k in range(N):s+=2*z/(2*k+1);z*=y*y
 return s,s+2*z/((2*N+1)*(1-y*y))
def bounds(x,D):
 # Reconstruct producer 48 terms and separately contain them within a different 40-term interval.
 s,t=loginterval(x,48);w,z=loginterval(x,40);assert w<=s<=t<=z
 l,u=floor(s,D),ceil(t,D);assert l<=s<=t<=u
 return l,u
d=F(7,25);n=128;q0=F(18,25);q1=F(24,25);gap=F(1,2**30);e=F(1,2**20);eret=F(1,4096);ea=F(1,65536);P0=1-e-F(1280,2**48);ul=P0*(1-e)/(1+e);uh=(1+e)/(P0*(1-e));sig=F(3,50);D=2**200
assert fixed['constant']=={'d':str(d),'J':128,'sigma':str(sig),'Ulo':str(ul),'Uhi':str(uh)} and len(fixed['rows'])==n
sumg=F(0);mineta=None;minsig=None
for j,row in enumerate(fixed['rows']):
 l=q0+(q1-q0)*F(j,n);u=q0+(q1-q0)*F(j+1,n);rl=l+d+gap;rh=u+d-gap;assert l<u and rl<rh
 ll,lh=bounds(rh/rl,D);_,ph=bounds(1/(1-d/(2*(l+d))),D);mu=ceil((l+d)/l*ph**2/(1-e),D)
 chi=F(-(-(mu*10**6).numerator//(mu*10**6).denominator),10**6)+F(1,10000);eta=(1-u)/(2*u);sm=sig/2-chi*(1-sig);rowb=(1-eret)*2*ul/u-2*e;M=2*uh/l;assert rowb>0 and eta>chi and sm>0
 gamma=floor((1-e)*(1-ea)*ll*rowb/M,D);assert gamma>0
 checks={'qlo':l,'qhi':u,'rlo':rl,'rhi':rh,'log_lower':ll,'log_upper':lh,'phi_upper':ph,'mu_upper':mu,'chi':chi,'eta':eta,'eta_margin':eta-chi,'sigma_margin':sm,'rowb':rowb,'gamma_lower':gamma}
 assert row['j']==j and all(F(row[k])==v for k,v in checks.items());sumg+=gamma;mineta=eta-chi if mineta is None else min(mineta,eta-chi);minsig=sm if minsig is None else min(minsig,sm)
assert F(fixed['gamma_total_lower'])==sumg and F(fixed['gamma_minus_3_47_minus_3_20_lower'])==sumg-F(3,47)-F(3,20)
assert sumg>F(429,2000) and F(gate['gamma_lower'])==sumg and F(gate['gamma_threshold'])==F(429,2000) and F(gate['gamma_margin_exact'])==sumg-F(429,2000)
raw=F(429,2000)-F(3,47);assert raw==F(14163,94000) and raw-F(3,20)==F(63,94000) and F(gate['raw_coefficient'])==raw and F(gate['raw_minus_3_20'])==F(63,94000)
# The ideal frontier certificate has no finite-error terms or trims.
t=top['top_parameters'];assert all(t[k]==summary['top64'][0][k] for k in ('d','qhi','J','i','k'));dd=F(t['d']);qh=F(t['qhi']);J=t['J'];assert J==512 and len(top['rows'])==J;ql=1-dd;assert 0<ql<qh<1;D=2**160
terms=F(0);mus=[];mine=None
for j,row in enumerate(top['rows']):
 l=ql+(qh-ql)*F(j,J);u=ql+(qh-ql)*F(j+1,J);R=l+dd;assert l<u
 ll,_=bounds((u+dd)/R,D);_,phi=bounds(1/(1-dd/(2*R)),D);mu=ceil(R/l*phi**2,D);eta=(1-u)/(2*u);term=floor(ll*l/u,D);assert eta>mu and term>=0
 checks={'qlo':l,'qhi':u,'mu_upper':mu,'eta':eta,'term_lower':term,'eta_margin':eta-mu};assert row['j']==j and all(F(row[k])==v for k,v in checks.items());terms+=term;mus.append(mu);mine=eta-mu if mine is None else min(mine,eta-mu)
c=terms-2*max(mus);assert c>0 and F(top['c_model_lower'])==c and F(top['min_eta_minus_mu_lower'])==mine
out={'scope':'128 finite-error layers and 512 ideal-model layers, exact arithmetic only; no operator theorem or global optimality','checks':{'128_layers_match':True,'128_eta_sigma_guards':True,'gamma_gt_429_over_2000':True,'raw_identity':True,'512_layers_match':True,'512_ideal_mu_eta_and_positive':True,'negative_and_positive_rounding':True,'40term_encloses_48term_each_log':True},'fixed_gamma_lower':str(sumg),'gamma_margin':str(sumg-F(429,2000)),'raw':str(raw),'raw_minus_3_20':str(raw-F(3,20)),'min_eta_margin':str(mineta),'min_sigma_margin':str(minsig),'ideal_c_lower':str(c),'ideal_min_margin':str(mine),'decimal':{'fixed_gamma':float(sumg),'gamma_margin':float(sumg-F(429,2000)),'min_eta_margin':float(mineta),'min_sigma_margin':float(minsig),'ideal_c_lower':float(c),'ideal_min_margin':float(mine)},'registration_sha256':sha((here/'registration.json').read_bytes()),'script_sha256':sha(Path(__file__).read_bytes()),'pid':os.getpid(),'cpu_seconds':time.process_time()-cpu,'wall_seconds':time.monotonic()-st,'peak_rss_bytes':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
(here/'audit.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({'checks':out['checks'],'decimal':out['decimal'],'pid':out['pid'],'cpu':out['cpu_seconds'],'wall':out['wall_seconds']},indent=2))
