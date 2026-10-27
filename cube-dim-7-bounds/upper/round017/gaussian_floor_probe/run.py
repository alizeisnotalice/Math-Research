#!/usr/bin/env python3
import os
for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS'):os.environ[k]='1'
import numpy as np,scipy,scipy.special as sp
from pathlib import Path
from fractions import Fraction
import hashlib,json,itertools,time,datetime,resource,sys
D=Path(__file__).resolve().parent;C=Path('/Users/zhengzhihao/Documents/Codex/2026-10-27/n-n-5-1-1-verify/work/cache/gaussian-round001/original-8192-shapes.npz');EXPECTED='6dfa48bf8954a7d02270af5bda9c1d33d4fcaf0ba4f5e486589186b7f7c09552';now=lambda:datetime.datetime.now(datetime.timezone.utc).isoformat()
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def fr(d):return Fraction(int(d[0]))*Fraction(2)**int(d[1])
def main():
 start=now();t=time.perf_counter();cpu=time.process_time();assert sha(C)==EXPECTED
 old=json.loads((D/'r016-manifest.json').read_text());assert sha(D/'r016-result.json')==old['files']['result.json'];prior=json.loads((D/'r016-result.json').read_text());gd=prior['enclosures']['g_required_B_RB'];gref=float((fr(gd['lower'])+fr(gd['upper']))/2)
 f=np.load(C); B=f['A_inverse'];d=f['det_A_minus_half'];L=f['scaled_A_cholesky'];scale=f['A_scale_wmax'];kd=float(f['potential_C']*f['radial_D']);ind=float(np.exp(sp.gammaln(32.5)-sp.gammaln(32)-3.5*np.log(np.pi)));assert abs(kd/ind-1)<5e-14
 assert B.shape==(8192,7,7) and d.shape==(8192,) and L.shape==(8192,7,7) and scale.shape==(8192,) and np.isfinite(B).all() and np.isfinite(d).all() and np.isfinite(L).all() and np.isfinite(scale).all() and (d>0).all() and (scale>0).all()
 asym=np.max(np.abs(B-np.swapaxes(B,1,2)),axis=(1,2))/np.max(np.abs(B),axis=(1,2));assert float(asym.max())<1e-10
 eig=np.linalg.eigvalsh(B);assert np.isfinite(eig).all() and np.all(eig[:,0]>0);condition=eig[:,-1]/eig[:,0];assert np.all(eig[:,-1]<=1/2048*(1+1e-10))
 V=np.array(list(itertools.product((-1.,1.),repeat=7)));assert V.shape==(128,7)
 qmax=np.empty(8192);tr=np.trace(B,axis1=1,axis2=2);assert np.isfinite(tr).all() and np.all(tr>0);terminals=[]
 for block in range(8):
  sl=slice(1024*block,1024*(block+1));q=np.einsum('vi,nij,vj->nv',V,B[sl],V,optimize=False);assert np.isfinite(q).all() and (q>0).all();qmax[sl]=q.max(axis=1);assert np.all(qmax[sl]>=tr[sl]*(1-1e-12))
  terminals.append(dict(block=block,status='completed',start=block*1024,end=(block+1)*1024,shapes=1024,vertex_quadratic_evaluations=1024*128,trace_evaluations=1024))
 lv=kd*d/qmax**3;lt=kd*d/tr**3;bound=kd/(np.sqrt(2048)*343);assert np.isfinite(lv).all() and np.isfinite(lt).all() and (lv>0).all() and (lt>0).all() and np.all(lt>=lv*(1-1e-12)) and np.all(lt<=bound*(1+1e-10)) and np.all(lv<=bound*(1+1e-10))
 check=[]
 for i in np.argsort(condition)[-8:]:
  i=int(i);u=np.linalg.solve(L[i],V.T);q=np.sum(u*u,axis=0)/scale[i];eye=np.linalg.solve(L[i],np.eye(7));stabletr=np.sum(eye*eye)/scale[i];stabled=1/(np.prod(np.diag(L[i]))*scale[i]**3.5);qs=np.einsum('vi,ij,vj->v',V,B[i],V);dq=float(np.max(abs(q-qs)/q));dt=abs(stabletr-tr[i])/stabletr;dd=abs(stabled-d[i])/stabled;assert max(dq,dt,dd)<=1e-6
  check.append(dict(shape_index=i,condition_estimate=float(condition[i]),max_vertex_relative_discrepancy=dq,trace_relative_discrepancy=float(dt),detweight_relative_discrepancy=float(dd),vertices_recomputed=128))
 np.savez_compressed(D/'observables.npz',L_vertex=lv,L_trace=lt,max_vertex_Q=qmax,trace_B=tr)
 def stats(x):
  blocks=x.reshape(8,1024).mean(axis=1);return dict(mean=float(x.mean()),sample_standard_deviation=float(x.std(ddof=1)),nominal_shape_SE=float(x.std(ddof=1)/np.sqrt(8192)),block_means=blocks.tolist(),nominal_block_SE=float(blocks.std(ddof=1)/np.sqrt(8)),min=float(x.min()),max=float(x.max()),ratio_to_required_g=float(x.mean()/gref))
 out=dict(status='completed',cache_sha256=sha(C),new_shape_samples=0,shape_count=8192,vertices_per_shape=128,vertex_quadratic_evaluations=8192*128,trace_evaluations=8192,same_input_cross_checks=check,cross_check_vertex_evaluations=8*128,cross_check_trace_evaluations=8,kd=kd,kd_independent_relative_discrepancy=abs(kd/ind-1),uniform_bound=bound,g_required_reference=gref,vertex=stats(lv),trace=stats(lt),max_condition_estimate=float(condition.max()),max_relative_asymmetry=float(asym.max()),min_B_eigenvalue=float(eig[:,0].min()),max_B_eigenvalue=float(eig[:,-1].max()),blocks=terminals,interpretation='empirical finite-cache diagnostic, no certified expectation error, global g or Cn')
 (D/'results.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n');r=dict(status='completed',start_utc=start,finish_utc=now(),wall_seconds=time.perf_counter()-t,cpu_seconds=time.process_time()-cpu,peak_rss_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,python=sys.version,numpy=np.__version__,scipy=scipy.__version__,shape_count=8192,block_terminals=8,MC_new=0,cache_read_only=True,threads=1,cross_checks=8);(D/'receipt.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(r,indent=2));print(json.dumps({'vertex':out['vertex'],'trace':out['trace'],'uniform_bound':bound,'g_reference':gref,'cross_checks':check},indent=2))
if __name__=='__main__':main()
