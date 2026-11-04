from fractions import Fraction as F
import multiprocessing as mp
from pathlib import Path
import json,os,time
from certify import case
BASE=Path(__file__).resolve().parent

def run(n,j):
 u=F(1,2**j);k=F(8*n,3);lam=F(11,10)*k**3/(8*n*n);a=(2*u)**n;v=2*u;delta=u/(1+2*u);w=delta+lam*u**6;rho=1-2*w;e=k*u**3;d=1-e**2;mout=(2*w)**n;M=a*rho**n;m0=(2*delta)**n;N=2**(7*j)
 row=case(n,k,2**j);row['j']=j;row['N']=str(N);row['N_rule']='ceil(u^-7)';row['definitions']={'cell':'T([-N,N]^n), T=[[1,e],[e,1]] direct sum I_(n-2)','good':'G=outer_cube minus the non-strict-witness set','selection':'A=G intersect {x1<c}, chosen so |A|=M; c>0 follows from G subset outer cube and M>mout/2','finite_input':'E_N=T([-N,N]^n) minus union over z in Z^n of (Tz+A)'}
 row['center_in_A_conditions']={'M_gt_mout_over_2':2*M>mout,'center_strict_witness':'checked below; G inclusion in outer cube and slice-existence are analytic inputs'}
 if row.get('budget_sign')!='positive':row['finite_status']='budget_or_guards_not_certified';return row
 radius=rho+w;neighbor=2*e**2*(2*w)**(n-2);defect_upper=(M+neighbor)/(2*radius)**n;strict_margin=u**n-defect_upper
 finite_gap=F(2*N-5,2*N)**n*M/(d-M)-m0/(1-m0)
 row['center']={'radius':str(radius),'defect_upper':str(defect_upper),'strict_margin':str(strict_margin),'strict_margin_positive':strict_margin>0,'interpretation':'actual center defect <= defect_upper; mean >= 1-defect_upper, not equality'}
 row['finite_halo_minus_Ln']={'exact_fraction':str(finite_gap),'strictly_positive':finite_gap>0,'finite_halo_formula':'1+((2N-5)/(2N))^n M/(d-M)'}
 row['finite_status']='certified' if 2*M>mout and strict_margin>0 and finite_gap>0 else 'finite_check_failed';return row

def worker(n):
 start=time.perf_counter();attempts=[]
 for j in (16,18):
  row=run(n,j);attempts.append(row)
  if row['finite_status']=='certified':break
 return {'n':n,'pid':os.getpid(),'elapsed':time.perf_counter()-start,'attempts':attempts}
def dispatch(n,barrier,queue):
 barrier.wait();queue.put(worker(n))

if __name__=='__main__':
 t=time.perf_counter();print(json.dumps({'master':os.getpid(),'workers_requested':15,'n':list(range(2,17)),'first_u':'2^-16','optional_retry':'2^-18'}),flush=True)
 ctx=mp.get_context('spawn');barrier=ctx.Barrier(15);queue=ctx.Queue();processes=[ctx.Process(target=dispatch,args=(n,barrier,queue)) for n in range(2,17)]
 for p in processes:p.start()
 res=sorted((queue.get() for _ in processes),key=lambda x:x['n'])
 for p in processes:p.join();assert p.exitcode==0
 data={'summary':{'cases':len(res),'worker_pids':sorted(set(x['pid'] for x in res)),'actual_worker_count':len(set(x['pid'] for x in res)),'elapsed_seconds':time.perf_counter()-t,'certified':sum(x['attempts'][-1]['finite_status']=='certified' for x in res),'failed':sum(x['attempts'][-1]['finite_status']!='certified' for x in res)},'cases':res}
 target=BASE/'center_finite.json';target.write_text(json.dumps(data,indent=2)+'\n');print(json.dumps(data['summary']),flush=True)
