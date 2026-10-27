#!/usr/bin/env python3
import os
for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS'):os.environ[k]='1'
from flint import arb,acb,acb_series,fmpq,ctx
import flint,sys,json,math,time,datetime,resource,subprocess,threading
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor,as_completed
D=Path(__file__).resolve().parent
ctx.prec=512;ctx.cap=129
now=lambda:datetime.datetime.now(datetime.timezone.utc).isoformat()
def endpoint(a):
 lo=a.lower().man_exp();hi=a.upper().man_exp();return {'lower':[str(lo[0]),int(lo[1])], 'upper':[str(hi[0]),int(hi[1])]}
def panel(j):
 t=time.perf_counter();c=time.process_time()
 try:
  a=arb(fmpq(1,16));center=arb(fmpq(2*j+1,32));z=acb_series([acb(a,center),acb(0,1)])
  F=z.gamma()*(4+z).gamma()*(4+8*z/7).rgamma()*(arb(fmpq(1,2))+z/7).gamma()**8*((-4*z/7)*acb(arb(2).log())).exp()/arb(fmpq(1,2)).gamma()**8
  if j<4:
   zz=acb(a,center); direct=zz.gamma()*(4+zz).gamma()*(4+8*zz/7).rgamma()*(arb(fmpq(1,2))+zz/7).gamma()**8*((-4*zz/7)*acb(arb(2).log())).exp()/arb(fmpq(1,2)).gamma()**8
   if not (F[0]-direct).contains(0) or not ((zz+1).gamma()-zz*zz.gamma()).contains(0):raise ArithmeticError('pilot direct or gamma recurrence overlap failed')
  conj=acb_series([F[k].conjugate() for k in range(129)]); prod=F*conj
  if not all(prod[k].is_finite() and prod[k].imag.contains(0) for k in range(129)):raise ArithmeticError('nonfinite or imaginary coefficient not enclosing zero')
  h=arb(fmpq(1,32));res=arb(0);power=h
  for k in range(0,129,2):res+=2*power*prod[k].real/(k+1);power*=h*h
  if not res.is_finite() or not res>0:raise ArithmeticError('panel integral not finite strictly positive')
  return dict(panel=j,status='completed',integral=endpoint(res),wall_seconds=time.perf_counter()-t,cpu_seconds=time.process_time()-c,peak_rss_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
 except BaseException as e:return dict(panel=j,status='failed',error=repr(e),wall_seconds=time.perf_counter()-t,cpu_seconds=time.process_time()-c,peak_rss_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
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
 assert not (D/'panels.jsonl').exists();start=now();t=time.perf_counter(); stop=threading.Event();th=threading.Thread(target=monitor,args=(stop,os.getpid()),daemon=True);th.start();results=[]
 def record(r):
  results.append(r)
  with (D/'panels.jsonl').open('a') as o:o.write(json.dumps(r)+'\n');o.flush();os.fsync(o.fileno())
 try:
  for i in range(4):
   r=panel(i);record(r)
   if r['status']!='completed':raise RuntimeError('warmup failed')
  with ProcessPoolExecutor(max_workers=15) as ex:
   fs={ex.submit(panel,i):i for i in range(4,512)}
   for f in as_completed(fs):
    try:r=f.result()
    except BaseException as e:r=dict(panel=fs[f],status='failed',error=repr(e))
    record(r)
  status='completed' if len(results)==512 and all(x['status']=='completed' for x in results) else 'failed'
 except BaseException as e:status='failed';print('ERROR',repr(e),file=sys.stderr)
 finally:stop.set();th.join()
 receipt=dict(status=status,start_utc=start,finish_utc=now(),wall_seconds=time.perf_counter()-t,sum_panel_cpu_seconds=sum(x.get('cpu_seconds',0) for x in results),max_worker_rss_bytes=max([x.get('peak_rss_bytes',0) for x in results],default=0),panel_terminal_count=len(results),pilot_count_included=4,python=sys.version,flint=flint.__version__,precision_bits=512,series_coefficients=129,thread_limit=1,processes=15)
 (D/'run-receipt.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt,indent=2));sys.exit(0 if status=='completed' else 1)
if __name__=='__main__':main()
