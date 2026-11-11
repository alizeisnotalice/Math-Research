"""Deterministic global optimized-Laplace curvature diagnostics, not interval proofs."""
import os
for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):os.environ[key]='1'
import json,multiprocessing as mp,time,math
from pathlib import Path
import mpmath as m

def job(j):
 rows=[]
 for wp in (80,160):
  m.mp.dps=wp+max(0,math.ceil(-j*math.log10(2)/16))+80;u=m.power(2,m.mpf(j)/16)
  mu=m.findroot(lambda c:c-1-m.exp(-c),m.mpf('1.28'));c=mu if u<=1 else m.log(2*u+2)/u
  lo=m.mpf(0);hi=mu
  for it in range(100):
   f=m.expm1(u*c)/u-1-m.exp(-c);df=m.exp(u*c)+m.exp(-c)
   if abs(f)<m.power(10,-wp-30):break
   if f>0:hi=c
   else:lo=c
   nxt=c-f/df
   if not lo<nxt<hi:nxt=(lo+hi)/2
   step=abs(nxt-c);c=nxt
   if step<m.power(10,-wp-30):break
  res=abs(m.expm1(u*c)/u-1-m.exp(-c));p=1/(1+m.exp(c));A=1+u-p;y=c*(1-p);v=u/(1-p)
  num=u-p+2*p*c*A-c*c*p*A*(1-p);d2=num/(u*(1+u)**2);lower=1/(2*(1+u)**2)
  direct=m.log1p(m.exp(-(1+u)*c))-m.log1p(m.exp(-c))-m.log1p(u)

  if not (res<m.power(10,-wp) and d2>=lower and abs(direct+u*c)<m.power(10,-wp)*u):raise ArithmeticError((j,wp,m.nstr(u,10),m.nstr(res,10),m.nstr(d2/lower,10),m.nstr(abs(direct+u*c)/u,10),it))
  rows.append({'wp':wp,'guarded_digits':m.mp.dps,'u':m.nstr(u,65),'c':m.nstr(c,65),'ell':m.nstr(-u*c,65),'curvature':m.nstr(d2,65),'lower_curvature':m.nstr(lower,65),'curvature_ratio':m.nstr(d2/lower,65),'root_residual':m.nstr(res,65),'stationary_value_error_over_u':m.nstr(abs(direct+u*c)/u,65)})
 return {'dyadic_exponent_numerator':j,'dyadic_exponent_denominator':16,'pid':os.getpid(),'rows':rows}
if __name__=='__main__':
 out=Path('/tmp/laplace94');out.mkdir(exist_ok=True);start=time.monotonic();n=0;pids=set()
 with (out/'rows.jsonl').open('w') as f,mp.get_context('spawn').Pool(12) as pool:
  for r in pool.imap_unordered(job,range(-4096,129),chunksize=1):
   f.write(json.dumps(r)+'\n');f.flush();n+=1;pids.add(r['pid'])
   if n%48==0:print(json.dumps({'done':n,'total':4225,'workers':len(pids),'seconds':time.monotonic()-start}),flush=True)
 (out/'summary.json').write_text(json.dumps({'pass_batch':True,'cases':n,'actual_workers':len(pids),'pids':sorted(pids),'seconds':time.monotonic()-start,'scope':'Global optimized Laplace curvature diagnostics; scalar proof separate, no CDF or physical-halo ordering inferred.'},indent=2)+'\n')
