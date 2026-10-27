#!/usr/bin/env python3
import os
for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS'):os.environ[k]='1'
from flint import arb,acb,acb_series,fmpq,ctx
import flint,sys,json,time,datetime,resource,subprocess,threading
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor,as_completed
D=Path(__file__).resolve().parent;ctx.prec=512;ctx.cap=129;now=lambda:datetime.datetime.now(datetime.timezone.utc).isoformat()
def ep(a):
 l=a.lower().man_exp();u=a.upper().man_exp();return {'lower':[str(l[0]),int(l[1])],'upper':[str(u[0]),int(u[1])]}
def G(z,j):
 e=arb(fmpq(1,7));ln2=arb(2).log();norm=arb(fmpq(1,2)).gamma()**8
 return z.gamma()*((-4*z/7)*acb(ln2)-acb(2*j*e*ln2)).exp()*(4+z-j*e).gamma()*(4+8*z/7-j*e).rgamma()*(arb(fmpq(1,2))+z/7-e).gamma()**j*(arb(fmpq(1,2))+z/7).gamma()**(8-j)/norm
def panel(key):
 j,k=key;t=time.perf_counter();c=time.process_time()
 try:
  a=arb(fmpq(9,16));center=arb(fmpq(2*k+1,4));z=acb_series([acb(a,center),acb(0,1)]);F=G(z,j);B=G(z,7-j);K=F*acb_series([B[v].conjugate() for v in range(129)])
  if not all(K[v].is_finite() for v in range(129)):raise ArithmeticError('nonfinite signed coefficients')
  if k<4:
   zz=acb(a,center);nzz=zz.conjugate();direct=G(zz,j)*G(nzz,7-j);direct_neg=G(nzz,j)*G(zz,7-j)
   if not (K[0]-direct).contains(0) or not (direct_neg-direct.conjugate()).contains(0) or not ((zz+1).gamma()-zz*zz.gamma()).contains(0):raise ArithmeticError('pilot constant/conjugacy/recurrence overlap failed')
  h=arb(fmpq(1,4));power=h;val=arb(0);terms=[]
  for n in range(0,129,2):
   term=2*power*K[n].real/(n+1);terms.append(ep(term));val+=term;power*=h*h
  if not val.is_finite():raise ArithmeticError('nonfinite signed integral')
  return dict(group=j,panel=k,status='completed',integral=ep(val),signed_integrated_even_terms=terms,constant_imag=ep(K[0].imag),negative_even_real_coefficients=sum(1 for n in range(0,129,2) if K[n].real<0),wall_seconds=time.perf_counter()-t,cpu_seconds=time.process_time()-c,peak_rss_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
 except BaseException as exc:return dict(group=j,panel=k,status='failed',error=repr(exc),wall_seconds=time.perf_counter()-t,cpu_seconds=time.process_time()-c,peak_rss_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
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
  stop.wait(2)
def main():
 assert not (D/'panels.jsonl').exists();start=now();t=time.perf_counter();stop=threading.Event();th=threading.Thread(target=monitor,args=(stop,os.getpid()),daemon=True);th.start();results=[]
 def record(r):
  results.append(r)
  with (D/'panels.jsonl').open('a') as f:f.write(json.dumps(r)+'\n');f.flush();os.fsync(f.fileno())
 try:
  for j in range(4):
   for k in range(4):
    r=panel((j,k));record(r)
    if r['status']!='completed':raise RuntimeError('pilot failed')
  with ProcessPoolExecutor(max_workers=15) as ex:
   fs={ex.submit(panel,(j,k)):(j,k) for j in range(4) for k in range(4,64)}
   for f in as_completed(fs):
    try:r=f.result()
    except BaseException as e:r=dict(group=fs[f][0],panel=fs[f][1],status='failed',error=repr(e))
    record(r)
  status='completed' if len(results)==256 and all(x['status']=='completed' for x in results) else 'failed'
 except BaseException as e:status='failed';print('ERROR',repr(e),file=sys.stderr)
 finally:stop.set();th.join()
 receipt=dict(status=status,start_utc=start,finish_utc=now(),wall_seconds=time.perf_counter()-t,sum_panel_cpu_seconds=sum(x.get('cpu_seconds',0) for x in results),max_worker_rss_bytes=max((x.get('peak_rss_bytes',0) for x in results),default=0),panel_terminal_count=len(results),pilot_count_included=16,python=sys.version,flint=flint.__version__,precision_bits=512,series_coefficients=129,thread_limit=1,processes=15)
 (D/'run-receipt.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt,indent=2));sys.exit(0 if status=='completed' else 1)
if __name__=='__main__':main()
