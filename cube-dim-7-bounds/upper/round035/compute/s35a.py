#!/usr/bin/env python3
"""S35A: expected exact face coefficients for typical shell; diagnostic only."""
import sys,os,time,json,hashlib,resource,concurrent.futures,math
from pathlib import Path
os.environ.update({k:'1' for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS')})
import mpmath as mp
import numpy as np
mp.mp.dps=90
LS=(16,32,64,128,256,1024);PS=('0.5','0.75','1')
NODES={k:np.polynomial.legendre.leggauss(k) for k in (64,128)}

def run(task):
 start=time.monotonic();cpu=time.process_time();L0,pstr=task;p=mp.mpf(pstr);m=math.ceil(L0/math.log(2));n=mp.power(2,m);L=m*mp.log(2);y=mp.sqrt(L);x=L-3*mp.log(L);assert x>0
 invn=1/n;r=1+mp.mpf('0.21')/(2*y);eta=1-mp.mpf('0.29')/y;d=1+mp.mpf('0.21')/y-eta;dr=r-eta
 # Uniform[-r,r] at a fixed negative face: n−1 iid coordinates and one pinned −r.
 # S/n = E[(|w|-eta)+]+O(1/n), K/n = Pr(|w|>eta)+O(1/n).
 ES=(1-invn)*dr**2/(2*r)+invn*dr; EK=(1-invn)*dr/r+invn
 markovS=ES/(4*d**2);markovK=EK/(8*d);markov_bad=min(mp.mpf(1),markovS+markovK)
 zl=mp.mpf('0.10')/L**p;zh=mp.mpf('0.29')/L**p;assert 0<zl<zh<1
 def integrand(u):
  # rho=t*u, n/(2rho)=x/(2u); all products use log1p so no 1−exp(−L) cancellation.
  f=x/(2*u);plus=(r-u)**2/(2*r);ec=(1-u)+(r-1)**2/(4*r)+(r-u)**2/(4*r)
  v=f*plus;vc=f*ec;v0=f*(r-u);vc0=f*(1+r-2*u)
  ellprod=mp.exp((n-2)*mp.log1p(-v*invn));cprod=mp.exp((n-2)*mp.log1p(-vc*invn))
  EN=1+u/r;EP=(r+u)/(2*r)
  full=mp.mpf('0.5')*(invn*ellprod*(1-v*invn)+(1-invn)*EN*(1-v0*invn)*ellprod)
  inter=mp.mpf('0.5')*(1-invn)*EP*(1-vc0*invn)*cprod
  return full/u,inter/u,(full-inter)/u
 def quad(k):
  ns,ws=NODES[k];out=[mp.mpf(0)]*3
  for node,weight in zip(ns,ws):
   u=1-(zl+zh)/2+(zh-zl)/2*mp.mpf(float(node)); vals=integrand(u)
   for j in range(3):out[j]+=mp.mpf(float(weight))*vals[j]
  return [z*(zh-zl)/2*y for z in out]
 a=quad(64);b=quad(128);rel=max(abs(a[i]-b[i])/(1+abs(b[i])) for i in range(3))
 # tau upper from R034 analytic n/L, normalized by n/y; avoid B asymptotic surrogate.
 return {'L0':L0,'m':m,'n_representation':'2^'+str(m),'actual_L':mp.nstr(L,20),'p':pstr,'a0':'.10','b0':'.29','r_minus_one':mp.nstr(r-1,20),'markov_ES_over_n':mp.nstr(ES,20),'markov_EK_over_n':mp.nstr(EK,20),'markov_union_upper':mp.nstr(markov_bad,20),'unconditional_face_mass_not_good_set':True,'shell_full_128_over_n_times_y':mp.nstr(b[0],24),'shell_corner_intersection_128_over_n_times_y':mp.nstr(b[1],24),'shell_net_128_over_n_times_y':mp.nstr(b[2],24),'two_rule_max_scaled_difference':mp.nstr(rel,12),'stable_under_1e-11':bool(rel<mp.mpf('1e-11')),'tau_R034_upper_over_n_times_y':mp.nstr(1/y,20),'quadrature_nodes_per_case':192,'pid':os.getpid(),'wall_seconds':time.monotonic()-start,'cpu_seconds':time.process_time()-cpu,'peak_rss_bytes':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,'scope':'unconditional face expectation, not uniform good-set bound or operator certificate'}
if __name__=='__main__':
 start=time.monotonic();tasks=[(L,p) for L in LS for p in PS]
 with concurrent.futures.ProcessPoolExecutor(max_workers=6) as ex:r=list(ex.map(run,tasks))
 out={'batch':'S35A','status':'numerical diagnostic only','model':'negative face -r, n−1 independent uniform coordinates; exact coordinate expectations of R033 local captured and original intersection','script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'workers':6,'precision_decimal_digits':90,'quadrature':'Gauss-Legendre 64 and 128, numpy float64 nodes and weights; 192 total per case','wall_seconds':time.monotonic()-start,'results':r}
 Path(__file__).with_name('s35a_results.json').write_text(json.dumps(out,indent=2)+'\n'); print(json.dumps([{'L0':z['L0'],'p':z['p'],'net':z['shell_net_128_over_n_times_y'],'markov_bad':z['markov_union_upper'],'tau_upper':z['tau_R034_upper_over_n_times_y'],'stable':z['stable_under_1e-11'],'pid':z['pid']} for z in r],indent=2))
