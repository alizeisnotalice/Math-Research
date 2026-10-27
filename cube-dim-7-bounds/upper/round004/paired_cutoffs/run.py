#!/usr/bin/env python3
import os
for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS'):os.environ[k]='1'
import numpy as np, scipy, scipy.special as sp, mpmath as mp
import json,hashlib,time,resource,datetime,sys,subprocess,threading
from pathlib import Path
from fractions import Fraction
from concurrent.futures import ProcessPoolExecutor,as_completed
D=Path(__file__).resolve().parent
E=np.array([1/256,1/1024,1/4096,1/16384]); NB=1024; CB=1024
BINS=np.array([0,1e-8,1e-7,1e-6,1e-5,1e-4,1e-3,1e-2,1,np.inf]); UB=np.array([0,1e-6,1e-4,1e-3,1e-2,.05,.1,.25,1.0000000001])
V=np.array([[1]+[1 if (i>>j)&1 else -1 for j in range(6)] for i in range(64)],float)
KD=np.exp(sp.gammaln(32.5)-sp.gammaln(32)-3.5*np.log(np.pi))
now=lambda:datetime.datetime.now(datetime.timezone.utc).isoformat()
def slow(y,x):
 # Re-evaluate same binary64 gamma draw at high precision, preserving all coordinates.
 with mp.workdps(80):
  yy=list(map(lambda z:mp.mpf(float(z)),y)); s=sum(yy);p=[z/s for z in yy];v=[[int(z) for z in row] for row in V];A=mp.matrix(7)
  for k in range(64):
   for i in range(7):
    for j in range(7):A[i,j]+=v[k][i]*v[k][j]/p[k]
  X=mp.matrix([mp.mpf(float(z)) for z in x]);q=(X.T*mp.lu_solve(A,X))[0]; det=mp.det(A)
  if q<=0 or det<=0:raise ValueError('high precision invalid Q or determinant')
  t=[]
  for e in E:
   z=q/(q+mp.mpf(float(e))**2);t.append(mp.gamma(mp.mpf('32.5'))/(mp.pi**mp.mpf('3.5')*mp.gamma(32))*det**(-mp.mpf('.5'))*q**-3*mp.betainc(3,mp.mpf('32.5'),0,z,regularized=True))
  return np.array([float(v) for v in t]),np.array([float(t[i+1]-t[i]) for i in range(3)])
def task(pid,rep,x):
 h=f'p{pid}s{rep}'; start=time.perf_counter(); cpu=time.process_time(); rng=np.random.default_rng(np.random.SeedSequence([2026102704,pid,rep])); sums=np.zeros(7); sq=np.zeros(7); subs=[]; vtot=np.zeros(7); bins_n=np.zeros((4,9),dtype=np.int64); bins_s=np.zeros((4,9)); un=np.zeros(8,dtype=np.int64); us=np.zeros((8,7)); fallback=[]; minp=1.;minq=float('inf'); minrat=1.
 try:
  for c in range(NB):
   y=rng.gamma(.5,1,(CB,64))
   if not np.isfinite(y).all() or (y<=0).any():raise ValueError(f'nonpositive gamma in chunk {c}; no redraw')
   p=y/y.sum(axis=1)[:,None]; minp=min(minp,float(p.min())); w=1/p; scale=w.max(axis=1); a=np.matmul(V.T[None,:,:]*((w/scale[:,None])[:,None,:]),V)
   try:L=np.linalg.cholesky(a); bad=np.zeros(CB,dtype=bool)
   except np.linalg.LinAlgError:
    L=np.zeros_like(a);bad=np.zeros(CB,dtype=bool)
    for k in range(CB):
     try:L[k]=np.linalg.cholesky(a[k])
     except np.linalg.LinAlgError:bad[k]=True;L[k]=np.eye(7)
   diag=np.diagonal(L,axis1=1,axis2=2); ratio=diag.min(axis=1)/diag.max(axis=1); minrat=min(minrat,float(ratio.min())); bad|=(ratio<1e-5)
   u=np.linalg.solve(L,np.broadcast_to(x[:,None],(CB,7,1))).squeeze(-1); Q=(u*u).sum(axis=1)/scale; det=np.prod(diag,axis=1)**-1*scale**-3.5
   bad|=~np.isfinite(Q)|~np.isfinite(det)|(Q<=0)|(det<=0);good=~bad
   T=np.zeros((CB,7)); z=Q[good,None]/(Q[good,None]+E[None,:]**2); ib=sp.betainc(3,32.5,z); comp=sp.betaincc(3,32.5,z); coef=KD*det[good]*Q[good]**-3; vals=coef[:,None]*ib; dif=np.where(ib[:,:-1]>.5,comp[:,:-1]-comp[:,1:],ib[:,1:]-ib[:,:-1])*coef[:,None];T[good]=np.column_stack((vals,dif))
   for k in np.flatnonzero(bad):
    vals,dif=slow(y[k],x);T[k]=np.r_[vals,dif];fallback.append({'chunk':c,'offset':int(k),'gamma_sha256':hashlib.sha256(y[k].tobytes()).hexdigest(),'reason':'Cholesky/diagonal ratio or nonfinite'} )
   if not np.isfinite(T).all() or (T[:,:4]<=0).any() or (T[:,4:]<0).any():raise ValueError(f'invalid kernel or band in chunk {c}')
   minq=min(minq,float(np.min(Q[good])) if good.any() else minq); ss=T.sum(axis=0);sums+=ss;sq+=(T*T).sum(axis=0);vtot+=ss
   if (c+1)%64==0:subs.append((vtot/(64*CB)).tolist());vtot[:]=0
   for j in range(4):
    ix=np.searchsorted(BINS,T[:,j],side='right')-1;bins_n[j]+=np.bincount(ix,minlength=9);bins_s[j]+=np.bincount(ix,weights=T[:,j],minlength=9)
   ix=np.searchsorted(UB,p[:,55]+p[:,63],side='right')-1;un+=np.bincount(ix,minlength=8)
   for j in range(7):us[:,j]+=np.bincount(ix,weights=T[:,j],minlength=8)
   if c in (0,63,255,511,767,1023):
    with (D/'progress.jsonl').open('a') as f:f.write(json.dumps(dict(handle=h,status='progress',chunks_completed=c+1,utc=now()))+'\n')
  return dict(handle=h,status='completed',seed=[2026102704,pid,rep],n=NB*CB,chunks=NB,sum=sums.tolist(),sum_squares=sq.tolist(),block_means=subs,kernel_bin_counts=bins_n.tolist(),kernel_bin_sums=bins_s.tolist(),u_bin_counts=un.tolist(),u_bin_sums=us.tolist(),min_p=minp,min_q_good=minq,min_cholesky_diag_ratio=minrat,fallbacks=fallback,wall_seconds=time.perf_counter()-start,cpu_seconds=time.process_time()-cpu,peak_rss_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
 except BaseException as e:
  return dict(handle=h,status='failed',seed=[2026102704,pid,rep],completed_chunks=c,error=repr(e),fallbacks=fallback,wall_seconds=time.perf_counter()-start,cpu_seconds=time.process_time()-cpu,peak_rss_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
def monitor(stop):
 while not stop.is_set():
  try:
   ps=subprocess.run(['ps','-axo','pid,ppid,pcpu,rss,command'],capture_output=True,text=True).stdout
   swap=subprocess.run(['sysctl','vm.swapusage'],capture_output=True,text=True).stdout.strip()
   with (D/'resource-samples.jsonl').open('a') as f:f.write(json.dumps(dict(utc=now(),ps=[s for s in ps.splitlines() if 'paired_cutoffs/run.py' in s or 'multiprocessing.spawn' in s],swap=swap))+'\n')
  except Exception as e:
   with (D/'resource-samples.jsonl').open('a') as f:f.write(json.dumps(dict(utc=now(),error=repr(e)))+'\n')
  stop.wait(5)
def main():
 pp=json.loads((D/'points.json').read_text())['points'];pts=[np.array([float(Fraction(s)) for s in pp[i]['exact_rational']]) for i in range(2)]+[np.array([float.fromhex(s) for s in pp[2]['exact_floathex']])];assert len(set(map(abs,pts[1])))==len(set(map(abs,pts[2])))==7
 start=now();t=time.perf_counter();stop=threading.Event();thread=threading.Thread(target=monitor,args=(stop,),daemon=True);thread.start(); results=[]
 try:
  with ProcessPoolExecutor(max_workers=15) as ex:
   fs={ex.submit(task,i,j,pts[i]):f'p{i}s{j}' for i in range(3) for j in range(5)}
   for f in as_completed(fs):
    try:r=f.result()
    except BaseException as e:r=dict(handle=fs[f],status='failed',error=repr(e))
    results.append(r)
    with (D/'terminal.jsonl').open('a') as o:o.write(json.dumps(r,allow_nan=False)+'\n');o.flush();os.fsync(o.fileno())
 finally:stop.set();thread.join()
 result=dict(status='completed' if len(results)==15 and all(x['status']=='completed' for x in results) else 'failed',start_utc=start,finish_utc=now(),wall_seconds=time.perf_counter()-t,sum_task_cpu_seconds=sum(r.get('cpu_seconds',0) for r in results),max_process_rss_bytes=max(r.get('peak_rss_bytes',0) for r in results),epsilons=E.tolist(),kernel_bins=BINS[:-1].tolist()+['inf'],u_bins=UB.tolist(),python=sys.version,numpy=np.__version__,scipy=scipy.__version__,results=sorted(results,key=lambda r:r['handle']))
 (D/'results.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n');print(result['status'],result['wall_seconds'],result['sum_task_cpu_seconds'],flush=True)
 if result['status']!='completed':sys.exit(1)
if __name__=='__main__':main()
