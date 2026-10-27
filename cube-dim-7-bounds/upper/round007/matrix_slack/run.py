#!/usr/bin/env python3
import os
for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS'):os.environ[k]='1'
import numpy as np, scipy, scipy.special as sp,mpmath as mp
import json,hashlib,time,resource,datetime,sys,subprocess,threading,multiprocessing as multiprocessing
from concurrent.futures import ProcessPoolExecutor,as_completed
from pathlib import Path
D=Path(__file__).resolve().parent; C=1024; NCH=1024;CAP=4096
V=np.array([[1]+[(-1)**((i>>k)&1) for k in range(6)] for i in range(64)],float)
code=sorted(sum(((r&u).bit_count()%2)<<(r-1) for r in range(1,7)) for u in range(8));unused=set(range(64));blocks=[]
while unused:
 q=min(unused);b=sorted(q^c for c in code);assert set(b)<=unused;blocks.append(b);unused-=set(b)
blocks=np.array(blocks); assert blocks.shape==(8,8) and np.all(np.matmul(np.swapaxes(V[blocks],1,2),V[blocks])==8*np.eye(7)[None,:,:])
pairs=blocks.reshape(4,16);Vp=V[pairs]; Dv=float(2*np.exp(sp.gammaln(32.5)-sp.gammaln(32)));upper=Dv/64
EDGES=[(0,1),(1,2),(2,3),(0,5),(5,1),(5,4),(4,2)]; PAIRS=[(i,j) for i in range(6) for j in range(i+1,6)];tol=2e-12*upper
now=lambda:datetime.datetime.now(datetime.timezone.utc).isoformat()
def slow(y):
 with mp.workdps(128):
  Y=[mp.mpf(float(z)) for z in y];total=sum(Y);w=[total/q for q in Y];vm=[[int(v) for v in row] for row in V];
  def matrix(ids):
   M=mp.matrix(7)
   for k in ids:
    for i in range(7):
     for j in range(7):M[i,j]+=w[k]*vm[k][i]*vm[k][j]
   return M
  bb=[matrix([int(x) for x in row]) for row in blocks];A=mp.matrix(7)
  for b in bb:A+=b
  r=[mp.det(b)**(mp.mpf(1)/7) for b in bb];rr=[mp.det(bb[2*i]+bb[2*i+1]) for i in range(4)];dd=2*mp.gamma(mp.mpf('32.5'))/mp.gamma(32)
  z=[dd*mp.det(A)**(-mp.mpf(1)/14),dd*sum(r)**(-mp.mpf(1)/2),dd/2*mp.fprod([(r[2*i]+r[2*i+1])**(-mp.mpf(1)/8) for i in range(4)]),dd/mp.sqrt(8)*mp.fprod([t**(-mp.mpf(1)/16) for t in r]),dd/2*mp.fprod([v**(-mp.mpf(1)/56) for v in rr]),dd*sum(v**(mp.mpf(1)/7) for v in rr)**(-mp.mpf(1)/2)]
  if any(z[i]>z[j] for i,j in EDGES) or any(x>dd/64 or x<=0 for x in z):raise ArithmeticError('high precision hierarchy failure')
  return np.array([float(v) for v in z])
def chol_matrices(a):
 try:return np.linalg.cholesky(a),np.zeros(a.shape[0],dtype=bool)
 except np.linalg.LinAlgError:
  L=np.zeros_like(a);bad=np.zeros(a.shape[0],bool)
  for k in range(a.shape[0]):
   try:L[k]=np.linalg.cholesky(a[k])
   except np.linalg.LinAlgError:bad[k]=True;L[k]=np.eye(7)
  return L,bad
def calc(y):
 p=y/y.sum(axis=1)[:,None];w=1/p;scale=w.max(axis=1);A=np.matmul(V.T[None,:,:]*(w/scale[:,None])[:,None,:],V);L,bad=chol_matrices(A);dg=np.diagonal(L,axis1=1,axis2=2);ratio=dg.min(axis=1)/dg.max(axis=1);bad|=ratio<1e-5;loga=2*np.log(dg).sum(axis=1)+7*np.log(scale)
 wp=w[:,pairs];spair=wp.max(axis=2);aa=np.matmul(np.swapaxes(Vp,1,2)[None,:,:,:]*(wp/spair[:,:,None])[:,:,None,:],Vp[None,:,:,:]);flat=aa.reshape(-1,7,7);Lp,bp=chol_matrices(flat);dp=np.diagonal(Lp,axis1=1,axis2=2).reshape(-1,4,7);ratio2=dp.min(axis=2)/dp.max(axis=2);bad|=bp.reshape(-1,4).any(axis=1)|(ratio2<1e-5).any(axis=1);logpair=2*np.log(dp).sum(axis=2)+7*np.log(spair)
 wb=w[:,blocks];logblock=18*np.log(2)+np.log(wb).sum(axis=2)+np.log((1/wb).sum(axis=2));r=np.exp(logblock/7);lr=np.log(r[:,::2]+r[:,1::2]);
 z=np.column_stack((Dv*np.exp(-loga/14),Dv/np.sqrt(r.sum(axis=1)),Dv/2*np.exp(-lr.sum(axis=1)/8),Dv/np.sqrt(8)*np.exp(-np.log(r).sum(axis=1)/16),Dv/2*np.exp(-logpair.sum(axis=1)/56),Dv/np.sqrt(np.exp(logpair/7).sum(axis=1))))
 bad|=(~np.isfinite(z)).any(axis=1)|(z<=0).any(axis=1);raw=np.zeros((len(y),8),bool)
 for k,(i,j) in enumerate(EDGES):raw[:,k]=z[:,i]>z[:,j]
 raw[:,7]=(z>upper).any(axis=1);bad|=raw.any(axis=1)
 return z,bad,raw,p,ratio,ratio2,logblock

def task(sid,semaphore,stop):
 t=time.perf_counter();cpu=time.process_time();rng=np.random.default_rng(np.random.SeedSequence([2026102707,sid]));su=np.zeros(21);sq=np.zeros(21);mn=np.full(21,np.inf);mx=np.full(21,-np.inf);partial=np.zeros(21);sub=[];fallback=[];rawcount=np.zeros(8,dtype=int);tolcount=np.zeros(8,dtype=int);minp=1.;minratio=1.;c=-1;completed=0
 try:
  for c in range(NCH):
   if stop.is_set():raise RuntimeError('shared stop before next chunk')
   y=rng.gamma(.5,1,(C,64))
   if not np.isfinite(y).all() or (y<=0).any():raise ValueError(f'zero/nonfinite gamma chunk {c}, no redraw')
   z,bad,raw,p,ratio,ratio2,lb=calc(y);bad|=y.max(axis=1)/y.min(axis=1)>2**40;minp=min(minp,float(p.min()));minratio=min(minratio,float(ratio.min()),float(ratio2.min()));rawcount+=raw.sum(axis=0)
   for k,(i,j) in enumerate(EDGES):tolcount[k]+=int((z[:,i]-z[:,j]>tol).sum())
   tolcount[7]+=int(((z-upper)>tol).any(axis=1).sum())
   if c==0:
    for k in range(4):
     for b in range(8):
      idx=blocks[b];w=1/p[k,idx];g=np.matmul(V[idx].T*w,V[idx]);logexact=2*np.log(np.diag(np.linalg.cholesky(g))).sum();
      if abs(logexact-lb[k,b])>1e-8:bad[k]=True
   for k in np.flatnonzero(bad):
    if not semaphore.acquire(False):raise RuntimeError(f'shared fallback cap {CAP} exhausted at chunk {c}, offset {int(k)}')
    fallback.append(dict(chunk=c,offset=int(k),gamma_sha256=hashlib.sha256(y[k].tobytes()).hexdigest(),gamma_floathex=[float(v).hex() for v in y[k]],raw_order=raw[k].tolist(),status='started'));vals=slow(y[k]);z[k]=vals;fallback[-1]['status']='completed'
   if not np.isfinite(z).all() or (z<=0).any() or any((z[:,i]>z[:,j]).any() for i,j in EDGES) or (z>upper*(1+1e-14)).any():raise ArithmeticError(f'final numerical hierarchy invalid chunk {c}')
   x=np.column_stack([z]+[(z[:,j]-z[:,i]) for i,j in PAIRS]);su+=x.sum(axis=0);sq+=(x*x).sum(axis=0);mn=np.minimum(mn,x.min(axis=0));mx=np.maximum(mx,x.max(axis=0));partial+=x.sum(axis=0);completed=c+1
   if completed%128==0:sub.append((partial/(128*C)).tolist());partial[:]=0
   if c in (0,127,383,639,895,1023):
    with (D/'progress.jsonl').open('a') as f:f.write(json.dumps(dict(handle=f's{sid:02}',chunks_completed=completed,utc=now()))+'\n')
  status='completed';error=None
 except BaseException as e:stop.set();status='failed';error=repr(e)
 return dict(handle=f's{sid:02}',status=status,error=error,seed=[2026102707,sid],completed_chunks=completed,completed_shapes=completed*C,sum=su.tolist(),sum_squares=sq.tolist(),minimum=[float(v) if np.isfinite(v) else None for v in mn],maximum=[float(v) if np.isfinite(v) else None for v in mx],subblock_means=sub,raw_violations=rawcount.tolist(),tolerance_violations=tolcount.tolist(),fallbacks=fallback,min_p=minp,min_cholesky_diag_ratio=minratio,wall_seconds=time.perf_counter()-t,cpu_seconds=time.process_time()-cpu,peak_rss_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
def monitor(stop,pid):
 while not stop.is_set():
  try:
   rows=subprocess.run(['ps','-axo','pid,ppid,pcpu,rss,command'],capture_output=True,text=True).stdout.splitlines();mine=[]
   for line in rows:
    bits=line.split(None,4)
    if len(bits)==5 and (bits[0]==str(pid) or bits[1]==str(pid)):mine.append(line)
   sw=subprocess.run(['sysctl','vm.swapusage'],capture_output=True,text=True).stdout.strip()
   with (D/'resources.jsonl').open('a') as f:f.write(json.dumps(dict(utc=now(),parent_pid=pid,processes=mine,swap=sw))+'\n')
  except Exception as e:
   with (D/'resources.jsonl').open('a') as f:f.write(json.dumps(dict(utc=now(),error=repr(e)))+'\n')
  stop.wait(5)
def main():
 assert not (D/'terminal.jsonl').exists();start=now();t=time.perf_counter();monstop=threading.Event();thread=threading.Thread(target=monitor,args=(monstop,os.getpid()),daemon=True);thread.start();res=[]
 try:
  with multiprocessing.Manager() as manager:
   sem=manager.BoundedSemaphore(CAP);stop=manager.Event()
   with ProcessPoolExecutor(max_workers=15) as ex:
    fs={ex.submit(task,i,sem,stop):i for i in range(15)}
    for f in as_completed(fs):
     try:r=f.result()
     except BaseException as e:stop.set();r=dict(handle=f's{fs[f]:02}',status='failed',error=repr(e))
     res.append(r)
     with (D/'terminal.jsonl').open('a') as o:o.write(json.dumps(r,allow_nan=False)+'\n');o.flush();os.fsync(o.fileno())
 finally:monstop.set();thread.join()
 status='completed' if len(res)==15 and all(r['status']=='completed' for r in res) else 'failed';R=dict(status=status,start_utc=start,finish_utc=now(),wall_seconds=time.perf_counter()-t,sum_task_cpu_seconds=sum(x.get('cpu_seconds',0) for x in res),peak_worker_rss_bytes=max((x.get('peak_rss_bytes',0) for x in res),default=0),task_count=len(res),completed_shapes=sum(x.get('completed_shapes',0) for x in res),fallback_count=sum(len(x.get('fallbacks',[])) for x in res),fallback_cap=CAP,upper_envelope=upper,tolerance=tol,labels=['Z0','Z1','Z2','Z3','Zpairmat','ZpairSum']+[f'Z{j}-Z{i}' for i,j in PAIRS],hierarchy_edges=EDGES,cosets=blocks.tolist(),python=sys.version,numpy=np.__version__,scipy=scipy.__version__,processes=15,blas_threads=1,results=sorted(res,key=lambda x:x['handle']))
 (D/'results.json').write_text(json.dumps(R,indent=2,allow_nan=False)+'\n');print(json.dumps({k:v for k,v in R.items() if k!='results'},indent=2));sys.exit(0 if status=='completed' else 1)
if __name__=='__main__':main()
