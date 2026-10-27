from fractions import Fraction as F
from pathlib import Path
import json,hashlib,os,time,resource
st=time.monotonic();cpu=time.process_time();here=Path(__file__).resolve().parent;reg=json.loads((here/'registration.json').read_bytes());sha=lambda b:hashlib.sha256(b).hexdigest()
for path,h in reg['source_sha256'].items():assert sha(Path(path).read_bytes())==h
source=Path('/Users/zhengzhihao/Desktop/T/cube-r030-compute/cube-dim-7-bounds/upper/round059/compute/wideband_certificate.json');orig=json.loads(source.read_bytes())
# For y=(x−1)/(x+1), log x=2∑ y^(2k+1)/(2k+1), with 0≤tail≤2y^(97)/(97(1−y²)).
def log_bounds(x):
 y=(x-1)/(x+1);assert 0<=y<1;z=y;S=F(0)
 for k in range(48):S+=2*z/(2*k+1);z*=y*y
 return S,S+2*z/(97*(1-y*y))
d=F(391,1905);a=F(283,2032);b=F(63,64);trim=F(1,2**20);chi=F(147,10000);sigma=F(3,100);er=F(1,4096);ea=F(1,65536);ej=en=es=et=ed=F(1,2**20);P0=1-ej-F(1280,2**48)
qlo=1-b*d;qhi=1-a*d;R=qlo+d;rl=R+trim;rh=qhi+d-trim
assert 0<a<b<1 and 0<d<1 and 0<rl<rh and 0<P0<1 and 0<chi<sigma<1
Ulo=P0*(1-en)/(1+ej);Uhi=(1+en)/(P0*(1-ej));M=2*Uhi/qlo;row=(1-er)*2*Ulo/qhi-es-et;assert row>0
ll,lh=log_bounds(rh/rl);mlo=(1-ea)*ll;mhi=(1-ea)*lh;gl=(1-ed)*mlo*row/M;gh=(1-ed)*mhi*row/M
pl,ph=log_bounds(1/(1-d/(2*R)));mul=R/qlo*pl**2/(1-ej);muh=R/qlo*ph**2/(1-ej);cap=a*d/(2*qhi)
assert gl-F(3,97)>F(1,10) and muh<chi<cap
expected={'d':d,'a':a,'b':b,'qmin':qlo,'qmax':qhi,'rlo':rl,'rhi':rh,'P0':P0,'Ulo':Ulo,'Uhi':Uhi,'M':M,'rowb':row,'gamma_lower':gl,'gamma_upper':gh,'mu_lower':mul,'mu_upper':muh,'chi':chi,'caprate':cap,'gamma_margin_lower':gl-F(3,97)-F(1,10),'mu_margin_lower':chi-muh,'cap_margin':cap-chi}
assert all(F(orig['exact_rational'][k])==v for k,v in expected.items())
out={'scope':'conditional scalar only, not operator or C_n proof','checks':{'valid_parameters':True,'positive_row':True,'gamma_minus_3_over_97_gt_1_over_10':True,'mu_upper_lt_chi_lt_caprate':True,'all_producer_exact_fields_match':True},'exact':{k:str(v) for k,v in expected.items()},'decimal':{k:float(v) for k,v in expected.items() if k in ('gamma_lower','gamma_margin_lower','mu_upper','chi','caprate','mu_margin_lower','cap_margin')},'registration_sha256':sha((here/'registration.json').read_bytes()),'script_sha256':sha(Path(__file__).read_bytes()),'pid':os.getpid(),'cpu_seconds':time.process_time()-cpu,'wall_seconds':time.monotonic()-st,'peak_rss_bytes':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
(here/'audit.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({'checks':out['checks'],'decimal':out['decimal'],'pid':out['pid'],'cpu':out['cpu_seconds'],'wall':out['wall_seconds']},indent=2))
