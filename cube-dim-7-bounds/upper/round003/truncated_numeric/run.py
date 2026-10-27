#!/usr/bin/env python3
import os
for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS'): os.environ[k]='1'
import numpy as np, scipy, scipy.special as sp
from scipy.stats import qmc
from concurrent.futures import ProcessPoolExecutor, as_completed
from fractions import Fraction
import hashlib,json,time,resource,datetime,sys
from pathlib import Path
D=Path(__file__).resolve().parent
C=Path('/Users/zhengzhihao/Documents/Codex/2026-10-27/n-n-5-1-1-verify/work/cache/gaussian-round001/original-8192-shapes.npz')
V=Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube-dim-7/cube-dim-7-bounds/upper/round002/g_probe/verification.json')
E=np.array([1/1024,1/256,1/64,1/16]); sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest(); now=lambda:datetime.datetime.now(datetime.timezone.utc).isoformat()
def load():
 global inv,det,L,scale,kd
 f=np.load(C); inv=f['A_inverse'];det=f['det_A_minus_half'];L=f['scaled_A_cholesky'];scale=f['A_scale_wmax'];kd=float(f['radial_D']*f['potential_C'])
 assert inv.shape==(8192,7,7) and np.isfinite(inv).all() and np.isfinite(det).all() and (det>0).all()
def task(i,pts,warm=False):
 t=time.perf_counter(); c=time.process_time(); out=[]; rel=0.; qmin=float('inf'); tmax=0.
 for p in pts:
  Q=np.einsum('i,nij,j->n',p,inv,p,optimize=False)
  if not np.isfinite(Q).all() or (Q<=0).any(): raise ValueError('invalid Q')
  if warm and len(out)<4:
   ids=np.linspace(0,8191,32,dtype=int); u=np.linalg.solve(L[ids],np.broadcast_to(p[:,None],(32,7,1))).squeeze(-1); q2=(u*u).sum(axis=1)/scale[ids];rel=max(rel,float(np.max(abs(Q[ids]-q2)/q2)))
   if rel>1e-8: raise ValueError('Q mismatch')
  qmin=min(qmin,float(Q.min())); z=Q[None,:]/(Q[None,:]+E[:,None]**2); terms=kd*det[None,:]*Q[None,:]**-3*sp.betainc(3,32.5,z)
  if not np.isfinite(terms).all() or (terms<=0).any() or (np.diff(terms,axis=0)>0).any(): raise ValueError('invalid terms')
  tmax=max(tmax,float(terms.max()));out.append(terms.reshape(4,8,1024).mean(axis=2))
 return dict(id=i,status='completed',points=len(pts),point_epsilon_shape_terms=len(pts)*4*8192,wall_seconds=time.perf_counter()-t,cpu_seconds=time.process_time()-c,process_peak_rss_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,min_q=qmin,max_shape_term=tmax,warm_cholesky_max_relative_error=rel,blocks=np.array(out))
def main():
 start=now(); t=time.perf_counter(); assert sha(C)=='6dfa48bf8954a7d02270af5bda9c1d33d4fcaf0ba4f5e486589186b7f7c09552'
 y=qmc.Sobol(6,scramble=True,seed=2026102703).random_base2(12); old=qmc.Sobol(6,scramble=False).random_base2(6); assert len(np.unique(y,axis=0))==4096 and not any(np.array_equal(a,b) for a in y for b in old)
 audits=json.loads(V.read_text())['audits']; assert len(audits)==8 and all(a['original_exact_floathex_point_in_S3'] for a in audits)
 oldpts=np.array([[float(Fraction(x)) for x in a['original_floathex_as_exact_rationals']] for a in audits]); pts=np.vstack((np.column_stack((np.ones(4096),y)),oldpts)); np.save(D/'points.npy',pts)
 load(); ledger=[]; blocks=np.empty((4104,4,8)); path=D/'task-ledger.jsonl'; assert not path.exists()
 def record(r):
  b=r.pop('blocks');i=r['id'];blocks[16*i:16*i+len(b)]=b;ledger.append(r)
  with path.open('a') as f:f.write(json.dumps(r,allow_nan=False)+'\n');f.flush();os.fsync(f.fileno())
 tasks=[(i,pts[16*i:min(16*i+16,4104)]) for i in range(257)];record(task(*tasks[0],warm=True))
 with ProcessPoolExecutor(max_workers=15,initializer=load) as ex:
  fs={ex.submit(task,*a):a[0] for a in tasks[1:]}
  for f in as_completed(fs):
   try:record(f.result())
   except BaseException as e:
    with path.open('a') as o:o.write(json.dumps(dict(id=fs[f],status='failed',error=repr(e)))+'\n')
    raise
 assert sorted(r['id'] for r in ledger)==list(range(257)); means=blocks.mean(axis=2); mins=[]
 for j,e in enumerate(E):
  i=int(np.argmin(means[:4096,j])); b=blocks[i,j]; ranges=(blocks[:4096,j].max(axis=1)-blocks[:4096,j].min(axis=1))/means[:4096,j]
  mins.append(dict(epsilon=float(e),lowest_sobol_index=i,point_floathex=[float(x).hex() for x in pts[i]],empirical_mean=float(means[i,j]),block_means=b.tolist(),block_min_over_mean=float(b.min()/b.mean()),block_max_over_mean=float(b.max()/b.mean()),all_sobol_block_relative_range_median=float(np.median(ranges)),all_sobol_block_relative_range_max=float(ranges.max())))
 oldout=[dict(handle=audits[j]['handle'],point_floathex=[float(x).hex() for x in pts[4096+j]],empirical_means=means[4096+j].tolist(),block_means=blocks[4096+j].tolist()) for j in range(8)]
 out=dict(status='completed',start_utc=start,finish_utc=now(),wall_seconds=time.perf_counter()-t,sum_task_cpu_seconds=sum(x['cpu_seconds'] for x in ledger),sum_task_wall_seconds=sum(x['wall_seconds'] for x in ledger),peak_single_process_rss_bytes=max(x['process_peak_rss_bytes'] for x in ledger),task_count=257,task_point_epsilon_shape_terms=sum(x['point_epsilon_shape_terms'] for x in ledger),workers_after_warmup=15,warmup_points_included_in_total=16,min_q=min(x['min_q'] for x in ledger),max_shape_term=max(x['max_shape_term'] for x in ledger),warm_cholesky_max_relative_error=ledger[0]['warm_cholesky_max_relative_error'],cache_sha256=sha(C),point_file_sha256=sha(D/'points.npy'),python=sys.version,numpy=np.__version__,scipy=scipy.__version__,minima=mins,old_S3_points=oldout,note='finite-cache diagnostics only; no deterministic expectation error, g or C7 bound')
 (D/'results.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n'); print(json.dumps({k:v for k,v in out.items() if k not in ('minima','old_S3_points')},indent=2))
if __name__=='__main__':main()
