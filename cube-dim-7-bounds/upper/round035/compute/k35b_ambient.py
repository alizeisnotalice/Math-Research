#!/usr/bin/env python3
"""14 ambient row/column full-kernel means; no pointwise or theorem claim."""
import os,sys,time,json,hashlib,resource,math,concurrent.futures
from pathlib import Path
os.environ.update({k:'1' for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS')})
import numpy as np
import mpmath as mp
mp.mp.dps=90
LS=(64,128,256,512,1024,4096,16384);ORDERS=(64,128,256);TASKS=[(L,k) for L in LS for k in ('row','column')]
NODES={k:np.polynomial.legendre.leggauss(k) for k in ORDERS}
def run(task):
 L0,kind=task;pid=os.getpid();start=time.monotonic();cpu=time.process_time();m=math.ceil(L0/math.log(2));n=mp.power(2,m);L=m*mp.log(2);x=L-3*mp.log(L);assert x>0;sqrtx=mp.sqrt(x);logt=L-mp.log(x);T=sqrtx*logt-1;assert T>0
 if kind=='row':fixedr=mp.exp(1/sqrtx)
 else:fixedq=mp.exp(-1/sqrtx)
 def G_over_n(v):
  if kind=='row':r=fixedr;q=mp.exp(-v/sqrtx)
  else:r=mp.exp(v/sqrtx);q=fixedq
  # n/(4rρ)=x/(4rq); use log1p(-z/n) to avoid cancellation for huge n.
  ndef=x*(r-q)**2/(4*r*q);fixeddef=x*(r-q)/(2*q)
  logM=mp.log1p(-ndef/n)
  product=mp.exp((n-2)*logM)
  return mp.mpf('0.5')*product*(mp.exp(logM)/n+(1-1/n)*(1+q/r)*(1-fixeddef/n))
 def integral(k):
  nodes,weights=NODES[k];zmax=T/(1+T);tot=mp.mpf('0')
  for node,weight in zip(nodes,weights):
   z=zmax*(mp.mpf(float(node))+1)/2;v=z/(1-z)
   tot+=mp.mpf(float(weight))*G_over_n(v)/(1-z)**2
  return tot*zmax/2
 estimates={str(k):mp.nstr(integral(k),30) for k in ORDERS};e=[mp.mpf(estimates[str(k)]) for k in ORDERS]
 return {'L0':L0,'kind':kind,'m':m,'n':'2^'+str(m),'L_actual':mp.nstr(L,25),'x':mp.nstr(x,25),'fixed_scaled_radius':1,'normalization':'sqrt(x)/n times ∫ G dlogrho for row, sqrt(x)/n times ∫ G dlogr for column','estimates':estimates,'scaled_diff_64_128':mp.nstr(abs(e[0]-e[1])/(1+abs(e[1])),15),'scaled_diff_128_256':mp.nstr(abs(e[1]-e[2])/(1+abs(e[2])),15),'stable_128_256_1e-10':bool(abs(e[1]-e[2])/(1+abs(e[2]))<mp.mpf('1e-10')),'status':'numeric only, not an interval enclosure or pointwise good-set certificate','pid':pid,'wall_seconds':time.monotonic()-start,'cpu_seconds':time.process_time()-cpu,'peak_rss_bytes':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
if __name__=='__main__':
 path=Path(__file__);start=time.monotonic();print('START parent_pid',os.getpid(),'workers=14 cases=14',flush=True)
 with concurrent.futures.ProcessPoolExecutor(max_workers=14) as pool:rows=list(pool.map(run,TASKS))
 out={'batch':'K35B','script_sha256':hashlib.sha256(path.read_bytes()).hexdigest(),'start_pid':os.getpid(),'workers_requested':14,'wall_seconds':time.monotonic()-start,'method':'64/128/256 Gauss-Legendre on exact finite domain via z=v/(1+v); 90dps arithmetic and float64 nodes/weights','results':rows};path.with_name('k35b_results.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps([{'L0':r['L0'],'kind':r['kind'],'value':r['estimates']['256'],'diff':r['scaled_diff_128_256'],'pid':r['pid'],'cpu':r['cpu_seconds'],'wall':r['wall_seconds']} for r in rows],indent=2))
