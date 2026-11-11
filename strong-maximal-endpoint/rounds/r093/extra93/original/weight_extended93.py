import os
for name in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):os.environ[name]='1'
import concurrent.futures as cf,hashlib,json,math,time
from pathlib import Path
import numpy as np
from scipy.optimize import brentq
from scipy.special import ndtr
MU=brentq(lambda u:u-1-math.exp(-u),1,2)
SAMPLES=2000000

def job(spec):
 n,m,w,c,rep=spec;seed=int.from_bytes(hashlib.sha256(repr(spec).encode()).digest()[:8],'big');rng=np.random.default_rng(seed);p=1/(1+math.exp(c))
 base=rng.gamma(n-m,1,SAMPLES)+c*rng.binomial(n-m,p,SAMPLES);part=rng.gamma(m,1,SAMPLES)+c*rng.binomial(m,p,SAMPLES)
 full=base+part;weighted=base+w*part;minimum=float(np.min(full-weighted));assert minimum>=-1e-10
 L=math.log(n);sigma=math.sqrt(MU);z=brentq(lambda x:ndtr(-x)-sigma*x/math.sqrt(n),1e-12,max(2,math.sqrt(L)+1));qs=sorted(set((.02,.05,.1,sigma*z/math.sqrt(n))))
 eps=math.sqrt(math.log(4e8)/(2*SAMPLES))+2/SAMPLES
 ts=sorted(set(t for q in qs for t in (max(0,q-eps),q,min(1,q+eps))));fd=dict(zip(ts,map(float,np.quantile(full,ts,overwrite_input=True))));wd=dict(zip(ts,map(float,np.quantile(weighted,ts,overwrite_input=True))))
 rows=[{'CDF_probability':q,'quantile_gap':fd[q]-wd[q],'gap_per_weight_deficit':(fd[q]-wd[q])/((1-w)*m),'DKW_gap_interval':[max(0,fd[max(0,q-eps)]-wd[min(1,q+eps)]),max(0,fd[min(1,q+eps)]-wd[max(0,q-eps)])]} for q in qs]
 return {'N':n,'m':m,'w':w,'c':c,'replicate':rep,'seed':seed,'pid':os.getpid(),'samples':SAMPLES,'CDF_error_with_interpolation_allowance':eps,'per_job_failure_bound':1e-8,'minimum_samplewise_gap':minimum,'rows':rows}

if __name__=='__main__':
 out=Path('/tmp/strong_weight_probe93b');out.mkdir(exist_ok=True);specs=[(n,m,w,c,r) for n in (4,8,16,32,64,128) for m in sorted(set((1,n//4,n//2))) for w in (.25,.5,.75,.9,.95,.98,.995) for c in tuple(j/8 for j in range(1,49))+(MU,) for r in range(4) if not (r<2 and w in (.5,.75,.9,.98) and c in (.5,1.,MU,2.,4.))];start=time.monotonic();pids=set();count=0
 with (out/'rows.jsonl').open('w') as log,cf.ProcessPoolExecutor(max_workers=12) as pool:
  for row in pool.map(job,specs,chunksize=1):
   log.write(json.dumps(row)+'\n');log.flush();count+=1;pids.add(row['pid'])
   if count%12==0:print(json.dumps({'completed':count,'total':len(specs),'workers':len(pids),'seconds':time.monotonic()-start}),flush=True)
 summary={'pass_batch':True,'jobs':count,'samples_per_job':SAMPLES,'total_sample_vectors':count*SAMPLES,'actual_workers':len(pids),'pids':sorted(pids),'seconds':time.monotonic()-start,'numpy':np.__version__,'family_failure_bound_union':count*1e-8,'scope':'Common-random-number independent weighted Gamma--Binomial quantile diagnostics; not a physical halo certificate or proof of uniform derivative bounds.'};(out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary),flush=True)
