#!/usr/bin/env python3
import os
for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS'):os.environ[k]='1'
from flint import fmpz_mat
import flint,json,itertools,math,time,datetime,resource,subprocess,threading,sys
from fractions import Fraction
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor,wait,FIRST_COMPLETED
D=Path(__file__).resolve().parent;I=json.loads((D/'inputs.json').read_text());V=I['body_vectors'];BODY=I['prefix_normals'];N=I['arrangement_normals']; CHUNK=1024;TOTAL=math.comb(36,6);now=lambda:datetime.datetime.now(datetime.timezone.utc).isoformat()
def chunk_task(taskid,start,combos,pilot=False):
 t=time.perf_counter();cpu=time.process_time();ranklow=0;mixed=0;good=0;rays={}; strata=[0]*7;cofactor_checked=0
 try:
  for offset,ids in enumerate(combos):
   mat=fmpz_mat([N[i] for i in ids]);rank=mat.rank();strata[sum(i<7 for i in ids)]+=1
   if pilot:
    cof=[(-1)**j*int(fmpz_mat([[N[i][k] for k in range(7) if k!=j] for i in ids]).det()) for j in range(7)];cofactor_checked+=1
    if (rank==6)!=any(cof):raise ArithmeticError("pilot cofactor rank disagreement")
   if rank<6:ranklow+=1;continue
   if rank!=6:raise ArithmeticError('invalid rank')
   kernel,dim=mat.nullspace()
   if dim!=1:raise ArithmeticError('rank/nullity mismatch')
   w=[int(kernel[i,0]) for i in range(7)];assert any(w) and all(sum(N[r][i]*w[i] for i in range(7))==0 for r in ids)
   if pilot:
    k=next(i for i,v in enumerate(cof) if v);assert all(w[i]*cof[k]==cof[i]*w[k] for i in range(7))
   if min(w)<0 and max(w)>0:mixed+=1;continue
   if max(w)<=0:w=[-v for v in w]
   if min(w)<0 or max(w)<=0:raise ArithmeticError('invalid orientation')
   g=math.gcd(*w);w=tuple(v//g for v in w);good+=1
   if w not in rays:rays[w]=list(ids)
  return dict(id=taskid,start=start,end=start+len(combos),status='completed',rank_below_6=ranklow,rank6_mixed= mixed,rank6_feasible=good,unique_local=len(rays),wall_strata=strata,pilot=pilot,cofactor_checked=cofactor_checked,wall_seconds=time.perf_counter()-t,cpu_seconds=time.process_time()-cpu,peak_rss_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,rays=[(list(w),v) for w,v in rays.items()])
 except BaseException as e:return dict(id=taskid,start=start,end=start+len(combos),status='failed',error=repr(e),offset=offset,rank_below_6=ranklow,rank6_mixed=mixed,rank6_feasible=good,wall_seconds=time.perf_counter()-t,cpu_seconds=time.process_time()-cpu,peak_rss_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,rays=[])
def monitor(stop,pid):
 while not stop.is_set():
  try:
   rows=subprocess.run(['ps','-axo','pid,ppid,pcpu,rss,command'],capture_output=True,text=True).stdout.splitlines();mine=[]
   for line in rows:
    bits=line.split(None,4)
    if len(bits)==5 and (bits[0]==str(pid) or bits[1]==str(pid)):mine.append(line)
   swap=subprocess.run(['sysctl','vm.swapusage'],capture_output=True,text=True).stdout.strip()
   with (D/'resources.jsonl').open('a') as f:f.write(json.dumps(dict(utc=now(),parent_pid=pid,processes=mine,swap=swap))+'\n')
  except Exception as e:
   with (D/'resources.jsonl').open('a') as f:f.write(json.dumps(dict(utc=now(),error=repr(e)))+'\n')
  stop.wait(5)
def main():
 assert not (D/'terminal.jsonl').exists();started=now();t=time.perf_counter();stop=threading.Event();th=threading.Thread(target=monitor,args=(stop,os.getpid()),daemon=True);th.start();rays={};ledger=[];source=itertools.combinations(range(36),6);nextid=8;nextstart=4096;middle_remaining=TOTAL-8192
 def nextjob():
  nonlocal nextid,nextstart,middle_remaining
  if not middle_remaining:return None
  batch=list(itertools.islice(source,min(CHUNK,middle_remaining)))
  assert batch
  job=(nextid,nextstart,batch,False);nextid+=1;nextstart+=len(batch);middle_remaining-=len(batch);return job
 def record(row):
  nonlocal rays
  points=row.pop('rays',[]);ledger.append(row)
  for w,ids in points:
   key=tuple(w);old=rays.get(key)
   if old is None or tuple(ids)<tuple(old):rays[key]=ids
  with (D/'terminal.jsonl').open('a') as f:f.write(json.dumps(row)+'\n');f.flush();os.fsync(f.fileno())
 failed=False
 try:
  for i in range(4):
   batch=list(itertools.islice(source,CHUNK));r=chunk_task(i,i*CHUNK,batch,True);record(r)
   if r['status']!='completed':raise RuntimeError('prefix pilot failed')
  suffix=itertools.islice(itertools.combinations(range(36),6),TOTAL-4096,None)
  for i in range(4):
   batch=list(itertools.islice(suffix,CHUNK));r=chunk_task(4+i,TOTAL-4096+i*CHUNK,batch,True);record(r)
   if r['status']!='completed':raise RuntimeError('suffix pilot failed')
  assert sum(x['cofactor_checked'] for x in ledger)==8192
  with ProcessPoolExecutor(max_workers=15) as ex:
   fs={}
   for _ in range(45):
    job=nextjob()
    if job:fs[ex.submit(chunk_task,*job)]=job[:2]
   while fs:
    done,_=wait(fs,return_when=FIRST_COMPLETED)
    for f in done:
     taskid,start=fs.pop(f)
     try:row=f.result()
     except BaseException as e:row=dict(id=taskid,start=start,status='failed',error=repr(e),rays=[])
     record(row)
     if row['status']!='completed':failed=True
     if not failed:
      job=nextjob()
      if job:fs[ex.submit(chunk_task,*job)]=job[:2]
  if failed:raise RuntimeError('one or more chunks failed')
  assert len(ledger)==1903 and nextstart==TOTAL-4096 and len({x['id'] for x in ledger})==1903
  spans=sorted((x['start'],x['end']) for x in ledger);assert spans[0][0]==0 and spans[-1][1]==TOTAL and all(spans[i][1]==spans[i+1][0] for i in range(len(spans)-1))
  assert [sum(x['wall_strata'][i] for x in ledger) for i in range(7)]==[475020,831285,498771,127890,14210,609,7]
  assert all(tuple(int(i==j) for i in range(7)) in rays for j in range(7))
  cat=[]
  for i,(w,ids) in enumerate(sorted(rays.items())):
   x=[sum(w[k] for k in range(i,7)) for i in range(7)];den=sum(x);assert den==sum((i+1)*w[i] for i in range(7)) and den>0
   active=[j for j,n in enumerate(N) if sum(n[k]*w[k] for k in range(7))==0];assert set(ids)<=set(active)
   signs=[(v>0)-(v<0) for n in BODY for v in [sum(n[k]*w[k] for k in range(7))]]
   minor=None
   for omit in range(7):
    cols=[k for k in range(7) if k!=omit];det=int(fmpz_mat([[N[j][k] for k in cols] for j in ids]).det())
    if det:minor=dict(columns=cols,determinant=det);break
   if minor is None:raise ArithmeticError('missing rank-six witness minor')
   cat.append(dict(id=i,w_primitive=list(w),xi_primitive=x,normalization_denominator=den,xi_normalized=[str(Fraction(v,den)) for v in x],active_normals=active,normal_incidence=[int(j in active) for j in range(36)],body_dot_signs=signs,witness_rows=ids,witness_minor=minor))
  (D/'rays.json').write_text(json.dumps(cat,indent=2)+'\n');status='completed'
 except BaseException as e:status='failed';print('ERROR',repr(e),file=sys.stderr)
 finally:stop.set();th.join()
 out=dict(status=status,start_utc=started,finish_utc=now(),wall_seconds=time.perf_counter()-t,sum_chunk_cpu_seconds=sum(x.get('cpu_seconds',0) for x in ledger),peak_worker_rss_bytes=max((x.get('peak_rss_bytes',0) for x in ledger),default=0),combination_total_expected=TOTAL,generated_combination_count=sum(x.get('end',x['start'])-x['start'] for x in ledger),terminal_chunks=len(ledger),processed_combination_count=sum(x.get('end',x['start'])-x['start'] for x in ledger if x['status']=='completed'),rank_below_6=sum(x.get('rank_below_6',0) for x in ledger),rank6_mixed=sum(x.get('rank6_mixed',0) for x in ledger),rank6_feasible=sum(x.get('rank6_feasible',0) for x in ledger),unique_rays=len(rays),python=sys.version,flint=flint.__version__,processes=15,threads_per_process=1,pilot_chunks_included=8, pilot_cofactor_checked=sum(x.get('cofactor_checked',0) for x in ledger),wall_strata=[sum(x.get('wall_strata',[0]*7)[i] for x in ledger) for i in range(7)])
 (D/'receipt.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2));sys.exit(0 if status=='completed' else 1)
if __name__=='__main__':main()
