"""Numerical exploration of specified continuous count-profile grid witnesses.
NOT a certified halo volume or general asymptotic bound. Single-thread BLAS.
"""
import os
for s in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS']:os.environ[s]='1'
import numpy as np,json,time,sys,math
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor
from scipy.special import gammaln,logsumexp
from scipy.stats import binom
HERE=Path(__file__).resolve().parent

def task(args):
 n,R,b,T=args;started=time.monotonic();p=2/(b+1);q=2/(min(T,R)+1);s=1 if T<=R else (R+1)/(T+1)
 pm=[binom.pmf(np.arange(i+1),i,p) for i in range(n+1)];qm=[binom.pmf(np.arange(i+1),i,q) for i in range(n+1)];lv=np.full(n+1,-np.inf);margin=np.full(n+1,np.inf);logfact=gammaln(np.arange(n+1)+1)
 for c in range(n+1):
  scale=s**c
  if scale<=.5:continue
  for j in range(n-c+1):
   a=n-j-c;dist=np.convolve(pm[j],qm[c]);tail=np.cumsum(dist[::-1])[::-1];av=np.zeros(n+1);av[:a+1]=scale
   if j+c:av[a+1:]=scale*tail[1:]
   keep=av>.5+1e-10
   weight=logfact[n]-logfact[a]-logfact[j]-logfact[c]+j*math.log(b-1)+c*math.log(T-b)
   lv[keep]=np.logaddexp(lv[keep],weight);margin[keep]=np.minimum(margin[keep],av[keep]-.5)
 kk=np.arange(n+1);terms=logfact[n]-logfact[kk]-logfact[n-kk]+kk*math.log(1/R)+(n-kk)*math.log1p(-1/R);source=np.array([logsumexp(terms[k:]) for k in range(n+1)])+n*math.log(R)
 score=lv-source;best=int(np.argmax(score[1:])+1)
 return dict(n=n,R=R,b=b,T=T,best_k=best,numerical_log_witness_ratio=float(score[best]),ratio_rate=float(score[best]/n),minimum_selected_average_margin=float(margin[best]),scope='Numerical grid inner approximation for explicitly chosen rectangles; no exact certificate and no full-halo claim',all_k_log_ratio=[float(x) if np.isfinite(x) else None for x in score],elapsed=time.monotonic()-started)
if __name__=='__main__':
 started=time.monotonic();jobs=[(n,R,b,T) for n in [8,16,32,64,128,256] for R in [4,8,16,64,256] for b in [math.sqrt(R),(R+1)/2] for T in [R,1.5*R]]
 with ProcessPoolExecutor(max_workers=12) as pool:res=list(pool.map(task,jobs))
 out=dict(status='complete',workers=12,cases=len(res),elapsed=time.monotonic()-started,seed='deterministic',average_threshold=.5,selection_margin=1e-10,certified=False,results=res)
 (HERE/'count_scan.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({k:out[k] for k in ['status','workers','cases','elapsed']}));print(json.dumps([max((r for r in res if r['n']==n),key=lambda r:r['ratio_rate']) for n in [8,16,32,64,128,256]],default=str)[:3000])
