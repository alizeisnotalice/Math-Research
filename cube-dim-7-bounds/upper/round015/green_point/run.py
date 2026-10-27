#!/usr/bin/env python3
import os
for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS'):os.environ[k]='1'
import green_kernel as kernel
from flint import arb,fmpq,ctx
from fractions import Fraction
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor,wait,FIRST_COMPLETED
import hashlib,json,time,datetime,resource,subprocess,threading,sys,flint
D=Path(__file__).resolve().parent;ctx.prec=256;old=json.loads((D/'r012-manifest.json').read_text())
for f,k in [('r012-result.json','result.json'),('cells.jsonl','cells.jsonl'),('r009-rays.json','r009-rays.json'),('r009-inputs.json','r009-inputs.json'),('r010-manifest.json','r010-manifest.json'),('testcases.json','testcases.json')]:assert hashlib.sha256((D/f).read_bytes()).hexdigest()==old['files'][k],f
oldres=json.loads((D/'r012-result.json').read_text());assert oldres['status']=='completed' and oldres['halfphases']==322560 and oldres['counts']['zero_ray_knots']==0 and oldres['point_floathex']==kernel.T['generic_case']['x_floathex'];x=[Fraction.from_float(float.fromhex(v)) for v in oldres['point_floathex']];assert len(set(abs(v) for v in x))==7 and all(v>0 for v in x)
K=['top_requests','unique_top','unique_subdd','dd_cache_hits','log_evaluations','hermite_evaluations','zero_ray_knots','repeated_top_nodes'];now=lambda:datetime.datetime.now(datetime.timezone.utc).isoformat()
def work(j):
 try:
  r=kernel.generic_chunk(j,256);rec=r['phase_records'];assert r['precision']==256 and [x['phase_id'] for x in rec]==list(range(64*j,64*(j+1)));counts={k:sum(x['counts'][k] for x in rec) for k in K};assert counts['top_requests']==64*212 and counts['zero_ray_knots']==0
  return dict(chunk=j,status='completed',start_phase=64*j,end_phase=64*(j+1),precision=256,interval=r['interval'],counts=counts,wall_seconds=r['wall_seconds'],cpu_seconds=r['cpu_seconds'],peak_rss_bytes=r['peak_rss_bytes'])
 except BaseException as e:return dict(chunk=j,status='failed',error=repr(e))
def monitor(stop,pid):
 while not stop.is_set():
  try:
   rows=subprocess.run(['ps','-axo','pid,ppid,pcpu,rss,command'],capture_output=True,text=True).stdout.splitlines();mine=[]
   for line in rows:
    s=line.split(None,4)
    if len(s)==5 and (s[0]==str(pid) or s[1]==str(pid)):mine.append(line)
   swap=subprocess.run(['sysctl','vm.swapusage'],capture_output=True,text=True).stdout.strip()
   with (D/'resources.jsonl').open('a') as f:f.write(json.dumps(dict(utc=now(),parent_pid=pid,processes=mine,swap=swap))+'\n')
  except Exception as e:
   with (D/'resources.jsonl').open('a') as f:f.write(json.dumps(dict(utc=now(),error=repr(e)))+'\n')
  stop.wait(5)
def main():
 assert not (D/'terminal.jsonl').exists();t=time.perf_counter();cpu=time.process_time();start=now();stop=threading.Event();th=threading.Thread(target=monitor,args=(stop,os.getpid()),daemon=True);th.start();rows=[];status='failed';error=None
 def record(x):
  rows.append(x)
  with (D/'terminal.jsonl').open('a') as o:o.write(json.dumps(x)+'\n');o.flush();os.fsync(o.fileno())
 try:
  cnt=kernel.counter();probe=kernel.dd_real((1,)*7,0,{},cnt);assert probe.contains(arb(fmpq(-1,6)))
  for j in range(4):
   r=work(j);record(r)
   if r['status']!='completed':raise RuntimeError('pilot failed')
  ps=sum((kernel.ball(r['interval']) for r in rows),arb(0));pw=kernel.width(ps);ptarget=arb(fmpq(1,10**12))*arb(fmpq(256,322560));passed=bool(pw<=ptarget)
  (D/'pilot.json').write_text(json.dumps(dict(status='passed' if passed else 'width_failed',chunks=[0,1,2,3],halfphases=256,same_node_logabs_probe='-1/6, exact integer knot=1',probe_counts=cnt,interval=kernel.ep(ps),width=kernel.ep(pw),target='1e-12*256/322560'),indent=2)+'\n')
  if not passed:raise RuntimeError('pilot width failed; no further chunks')
  pending=iter(range(4,5040));failed=False
  with ProcessPoolExecutor(max_workers=15) as ex:
   fs={}
   for _ in range(45):
    j=next(pending,None)
    if j is not None:fs[ex.submit(work,j)]=j
   while fs:
    done,_=wait(fs,return_when=FIRST_COMPLETED)
    for f in done:
     j=fs.pop(f)
     try:r=f.result()
     except BaseException as e:r=dict(chunk=j,status='failed',error=repr(e))
     record(r)
     if r['status']!='completed':failed=True
     if not failed:
      j=next(pending,None)
      if j is not None:fs[ex.submit(work,j)]=j
  if failed:raise RuntimeError('failed chunk; remaining undispatched, no retry')
  assert len(rows)==5040 and {x['chunk'] for x in rows}==set(range(5040)) and all(x['status']=='completed' for x in rows)
  total=sum((kernel.ball(x['interval']) for x in sorted(rows,key=lambda x:x['chunk'])),arb(0));wid=kernel.width(total);widthpass=bool(wid<=arb(fmpq(1,10**12)));positive=bool(total>0);status='completed' if widthpass and positive else 'width_failed' if not widthpass else 'positivity_unresolved'
  result=dict(status=status,full_G=kernel.ep(total),width=kernel.ep(wid),width_pass=widthpass,positive=positive,precision=256,point_floathex=oldres['point_floathex'],epsilon=0,chunks=5040,pilot_chunks=4,halfphases=322560,top_requests_before_cache=sum(x['counts']['top_requests'] for x in rows),counts={k:sum(x['counts'][k] for x in rows) for k in K},no_extra_normalization=True,not_a_C7_lower_bound=True)
  assert result['top_requests_before_cache']==68382720 and result['counts']['zero_ray_knots']==0
  (D/'result.json').write_text(json.dumps(result,indent=2)+'\n')
 except BaseException as e:error=repr(e)
 finally:stop.set();th.join()
 out=dict(status=status,error=error,start_utc=start,finish_utc=now(),wall_seconds=time.perf_counter()-t,parent_cpu_seconds=time.process_time()-cpu,sum_chunk_cpu_seconds=sum(x.get('cpu_seconds',0) for x in rows),peak_worker_rss_bytes=max((x.get('peak_rss_bytes',0) for x in rows),default=0),terminal_chunks=len(rows),new_G_chunks=sum(x.get('status')=='completed' for x in rows),precision_bits=256,workers=15,python=sys.version,flint=flint.__version__)
 (D/'receipt.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2));sys.exit(0 if status=='completed' else 1)
if __name__=='__main__':main()
