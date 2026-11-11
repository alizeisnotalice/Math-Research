import os
for s in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):os.environ[s]='1'
import concurrent.futures as cf,hashlib,json,math,time
from pathlib import Path
import numpy as np
from scipy.optimize import brentq,minimize_scalar
from scipy.special import gammainc,ndtr
from scipy.stats import binom
S=2000000
MU=brentq(lambda u:u-1-math.exp(-u),1,2)
def reference(n):
 z=brentq(lambda x:ndtr(-x)-math.sqrt(MU)*x/math.sqrt(n),1e-10,4);q=math.sqrt(MU)*z/math.sqrt(n);ks=np.arange(n+1)
 def quant(c):
  ws=binom.pmf(ks,n,1/(1+math.exp(c)))
  return brentq(lambda a:float(np.dot(ws,gammainc(n,np.maximum(a-c*ks,0))))-q,1e-12,n*(1+c)+20*math.sqrt(n),xtol=1e-10)
 opt=minimize_scalar(lambda c:-quant(c),bounds=(.001,6),method='bounded',options={'xatol':1e-7})
 return {'N':n,'c_ref':float(opt.x),'q_ref_deterministic':float(-opt.fun),'CDF_probability':q,'solver_success':bool(opt.success),'scope':'bounded scalar-search iid reference, not a certified global optimum'}
def job(spec):
 n,m,t,c1,c2,r,ref=spec;seed=int.from_bytes(hashlib.sha256(repr(spec).encode()).digest()[:8],'big');rng=np.random.default_rng(seed);c0=ref['c_ref'];p0=1/(1+math.exp(c0))
 def group(size,c):
  pc=1/(1+math.exp(c));lo=min(pc,p0);hi=max(pc,p0);k=rng.binomial(size,lo,S);extra=rng.binomial(size-k,(hi-lo)/(1-lo));kc=k+extra if pc>=p0 else k;kr=k if pc>=p0 else k+extra;g=rng.gamma(size,1,S);return g+c*kc,g+c0*kr
 ca,ra=group(m,c1);cb,rb=group(n-m,c2);v1=n*t/m;v2=n*(1-t)/(n-m);candidate=v1*ca+v2*cb;uniform=ra+rb;q=ref['CDF_probability'];eps=math.sqrt(math.log(4e8)/(2*S))+2/S;qs=(max(0,q-eps),q,min(1,q+eps));a=np.quantile(candidate,qs,overwrite_input=True);b=np.quantile(uniform,qs,overwrite_input=True)
 return {'N':n,'m':m,'t':t,'lambda1':v1,'lambda2':v2,'c1':c1,'c2':c2,'replicate':r,'seed':seed,'pid':os.getpid(),'samples':S,'candidate_quantile':float(a[1]),'reference_quantile':float(b[1]),'gain_vs_reference':float(a[1]-b[1]),'DKW_gain_interval':[float(a[0]-b[2]),float(a[2]-b[0])],'per_job_failure_bound':1e-8,'reference':ref}
if __name__=='__main__':
 out=Path('/tmp/normal_probe94');out.mkdir(exist_ok=True);refs={n:reference(n) for n in (8,16,32,64,128,256)};cs=(.25,.75,1.,MU,1.75,2.5,4.,6.);specs=[(n,m,t,c,d,r,refs[n]) for n in refs for m in sorted(set((1,n//8,n//4,n//2))) for t in sorted(set((m/n,.01,.025,.05,.1,.25,.5,.75,.9,.975,.99))) for c in cs for d in cs for r in (0,1)];start=time.monotonic();pids=set();num=0
 (out/'references.json').write_text(json.dumps(refs,indent=2)+'\n')
 with (out/'rows.jsonl').open('w') as f,cf.ProcessPoolExecutor(max_workers=12) as pool:
  for row in pool.map(job,specs,chunksize=1):
   f.write(json.dumps(row)+'\n');f.flush();num+=1;pids.add(row['pid'])
   if num%24==0:print(json.dumps({'completed':num,'total':len(specs),'workers':len(pids),'seconds':time.monotonic()-start}),flush=True)
 (out/'summary.json').write_text(json.dumps({'pass_batch':True,'cases':num,'samples_per_case':S,'actual_workers':len(pids),'pids':sorted(pids),'seconds':time.monotonic()-start,'family_failure_bound_union':num*1e-8,'scope':'normalized two-weight/two-action independent quantile search vs bounded scalar iid reference; no physical-halo certificate or asymptotic proof inferred'},indent=2)+'\n')
