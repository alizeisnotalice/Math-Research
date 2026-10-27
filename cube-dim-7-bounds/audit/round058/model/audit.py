from fractions import Fraction as F
from decimal import Decimal
from pathlib import Path
import json,hashlib,time,resource,os
st=time.monotonic();cpu=time.process_time();here=Path(__file__).resolve().parent;reg=json.loads((here/'registration.json').read_bytes());sha=lambda b:hashlib.sha256(b).hexdigest()
for path,h in reg['source_sha256'].items():assert sha(Path(path).read_bytes())==h
p=Path('/Users/zhengzhihao/Desktop/T/cube-r030-compute/cube-dim-7-bounds/upper/round058/compute');source_reg=json.loads((p/'registration.json').read_bytes());summary=json.loads((p/'summary.json').read_bytes());assert summary['registration_sha256']==sha((p/'registration.json').read_bytes()) and source_reg['code_sha256']==sha((p/'r058.py').read_bytes())
assert len(summary['receipts'])==128 and {x['file'] for x in summary['receipts']}=={f'chunk_{i:03}.json' for i in range(128)}
alltop=[];checked=eligible=positive=0;counts=[]
for z in summary['receipts']:
 f=p/z['file'];raw=f.read_bytes();assert sha(raw)==z['sha256'];x=json.loads(raw);i=x['index'];assert z['file']==f'chunk_{i:03}.json' and x['status']==z['status']=='completed' and x['error'] is None and x['code_sha256']==source_reg['code_sha256']
 assert x['checked']==z['checked']==128**2 and x['eligible']==z['eligible'] and x['positive']==z['positive'] and 0<=x['positive']<=x['eligible']<=x['checked'];assert len(x['top32'])==32
 for v in x['top32']:
  j,k=v['j'],v['k'];assert v['i']==i and 0<=j<128 and 0<=k<128
  assert F(v['d'])==F(1,5)+F(2*i,15*127) and F(v['a'])==F(1,16)+F(3*j,8*127) and F(v['b'])==F(9,16)+F(27*k,64*127)
 alltop+=x['top32'];checked+=x['checked'];eligible+=x['eligible'];positive+=x['positive'];counts.append({'i':i,'eligible':x['eligible'],'positive':x['positive']})
assert checked==128**3
alltop.sort(key=lambda v:Decimal(v['c_model']),reverse=True);assert alltop[:64]==summary['top64']
occupancy={i:sum(v['i']==i for v in summary['top64']) for i in range(128)}
# Independently evaluate five fixed extrema using exact 48-term atanh intervals.
def logbounds(x):
 z=(x-1)/(x+1);assert 0<=z<1;p=z;s=F(0)
 for n in range(48):s+=2*p/(2*n+1);p*=z*z
 return s,s+2*p/(97*(1-z*z))
tests=[]
for i,j,k in reg['five_points']:
 d=F(1,5)+F(2*i,15*127);a=F(1,16)+F(3*j,8*127);b=F(9,16)+F(27*k,64*127);q=1-b*d;Q=1-a*d;R=q+d;assert 0<a<b<1 and 0<q<Q
 l,u=logbounds(1/(1-d/(2*R)));ml=R/q*l*l;mu=R/q*u*u;eta=a*d/(2*Q);g0,g1=logbounds((Q+d)/R);cl=g0*q/Q-2*mu;cu=g1*q/Q-2*ml
 sign=lambda lo,hi:'positive' if lo>0 else 'negative' if hi<0 else 'uncertain'
 tests.append({'i':i,'j':j,'k':k,'mu_lt_eta':sign(eta-mu,eta-ml),'c_sign':sign(cl,cu),'mu_upper_approx':float(mu),'eta_approx':float(eta),'c_lower_approx':float(cl),'c_upper_approx':float(cu)})
out={'scope':'receipt and retained-top merge; 5 bounded point checks; global ranking conditional on producer per-chunk 80-digit results','chunks':128,'checked':checked,'eligible':eligible,'positive':positive,'retained_candidates':len(alltop),'top64_exact_match':True,'max_top64_chunk_occupancy':max(occupancy.values()),'five_points':tests,'chunk_counts':counts,'registration_sha256':sha((here/'registration.json').read_bytes()),'source_summary_sha256':sha((p/'summary.json').read_bytes()),'script_sha256':sha(Path(__file__).read_bytes()),'pid':os.getpid(),'wall_seconds':time.monotonic()-st,'cpu_seconds':time.process_time()-cpu,'peak_rss_bytes':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
(here/'audit.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({k:v for k,v in out.items() if k not in ('chunk_counts',)},indent=2))
