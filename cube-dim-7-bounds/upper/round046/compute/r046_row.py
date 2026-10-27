#!/usr/bin/env python3
"""R046: one-dimensional inverse-column-weight row diagnostics, no Cn claim."""
import os,sys,time,json,hashlib,resource,concurrent.futures,threading,subprocess
from pathlib import Path
from fractions import Fraction as Q
os.environ.update({k:'1' for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS')})
import mpmath as mp
D=Path(__file__).resolve().parent;DEN=(128,64,32,16,8);RHOS=(4,16,256);sha=lambda b:hashlib.sha256(b).hexdigest()
def reg():
 cases=[{'R':'1','g':str(Q(1,d)),'rho':str(Q(r)),'h_over_g':[str(Q(i,16)) for i in range(1,17)],'r_formula':'1+2*h*k/15 for k=0,...,15','combos':256} for d in DEN for r in RHOS]
 out={'batch':'R046','objective':'numerically cross-check J closed form vs direct piecewise integral, its face mean and second moment; not a C_n theorem','code_sha256':sha(Path(__file__).read_bytes()),'workers':15,'precision_dps':80,'cases':cases,'combo_total':3840,'per_unique_r_integrals':'one direct z integral, one face mean, one face second moment; all piecewise; identical r cached','stop':'256 registered combinations per task or stop case on first failing comparison/bound, saving exact parameters','monitor':'first 0.05s then 5s','registered_unix_time':time.time()};(D/'registration.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({'event':'registered','sha256':sha((D/'registration.json').read_bytes()),'code_sha256':out['code_sha256'],'tasks':15,'total_combos':3840}),flush=True)
def task(case):
 mp.mp.dps=80;start=time.monotonic();cpu=time.process_time();gq=Q(case['g']);rhoq=Q(case['rho']);g=mp.mpf(gq.numerator)/gq.denominator;rho=mp.mpf(rhoq.numerator)/rhoq.denominator;a=1-g;tol=mp.mpf('1e-65');cache={};rows=[];calls={'direct':0,'mean':0,'second':0};fail=None
 def E(m):
  if m>=g:return mp.mpf(0)
  q=g-m;return -2*mp.log1p(-q/2)-q
 e0=E(mp.mpf(0))
 def delta(w):
  w=abs(w);m=max(mp.mpf(0),w-a);return (E(m)+e0-m)/(2*rho)
 def direct(w):
  cuts=sorted(set([ -rho,rho,-rho+g,rho-g,max(-rho,min(rho,w-rho-a)),max(-rho,min(rho,w+rho+a)) ]))
  def f(z):
   ell=min(mp.mpf(1),(1+rho+a-abs(z))/2)
   return (1/ell if abs(z-w)<rho+a else 0)-1
  calls['direct']+=1
  return mp.fsum(mp.quad(f,[lo,hi]) for lo,hi in zip(cuts,cuts[1:]) if hi>lo)/(2*rho)
 def stats(rq):
  if rq in cache:return cache[rq]
  r=mp.mpf(rq.numerator)/rq.denominator;cuts=sorted(set([mp.mpf(0),min(r,a),min(r,mp.mpf(1)),r]))
  calls['mean']+=1;mean=mp.fsum(mp.quad(lambda w:delta(w),[lo,hi]) for lo,hi in zip(cuts,cuts[1:]) if hi>lo)/r
  calls['second']+=1;second=mp.fsum(mp.quad(lambda w:delta(w)**2,[lo,hi]) for lo,hi in zip(cuts,cuts[1:]) if hi>lo)/r
  Dv=direct(r);err=abs(Dv-delta(r));out=(mean,second,err,delta(r));cache[rq]=out;return out
 for i in range(1,17):
  hq=gq*Q(i,16);h=mp.mpf(hq.numerator)/hq.denominator
  for k in range(16):
   rq=Q(1)+2*hq*Q(k,15);mean,m2,err,point=stats(rq);lower=-6*h*g/rho;upper=31*g**3/rho**2
   flags={'piecewise_matches':bool(err<tol),'mean_upper':bool(mean<=tol),'mean_lower':bool(mean>=lower-tol),'second_bound':bool(m2<=upper+tol),'at_R_equal_one':bool(rq!=1 or abs(mean)<tol)}
   row={'h_over_g':str(Q(i,16)),'r':str(rq),'mean_minus_one':mp.nstr(mean,32),'second_moment':mp.nstr(m2,32),'point_delta':mp.nstr(point,32),'direct_abs_error':mp.nstr(err,8),'checks':flags};rows.append(row)
   if not all(flags.values()):fail={'parameters':{'R':'1','g':str(gq),'rho':str(rhoq),'h':str(hq),'r':str(rq)},'values':row};break
  if fail:break
 out={'R':'1','g':str(gq),'rho':str(rhoq),'status':'counterexample_or_numerical_mismatch' if fail else 'completed','failure':fail,'combos_completed':len(rows),'distinct_r':len(cache),'integral_calls':calls,'results':rows,'precision_dps':80,'pid':os.getpid(),'wall_seconds':time.monotonic()-start,'cpu_seconds':time.process_time()-cpu,'peak_rss_bytes':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,'scope':'numerical evidence only; no pointwise/global theorem or rigorous interval enclosure'}
 name='case_g'+str(gq.denominator)+'_rho'+str(rhoq.numerator)+'.json';(D/name).write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({'event':'case','g':str(gq),'rho':str(rhoq),'pid':out['pid'],'status':out['status'],'combos':len(rows),'integrals':calls,'wall':out['wall_seconds']}),flush=True);return {'file':name,'sha256':sha((D/name).read_bytes()),'status':out['status'],'pid':out['pid'],'wall_seconds':out['wall_seconds'],'cpu_seconds':out['cpu_seconds'],'peak_rss_bytes':out['peak_rss_bytes'],'integral_calls':calls}
def monitor(stop,parent):
 start=time.monotonic();delay=.05
 with (D/'monitor.jsonl').open('w') as f:
  while not stop.wait(delay):
   try:
    txt=subprocess.check_output(['ps','-axo','pid=,ppid=,%cpu=,rss=,time=,command='],text=True);ws=[]
    for line in txt.splitlines():
     z=line.strip().split(None,5)
     if len(z)==6 and int(z[1])==parent and 'spawn_main' in z[5]:ws.append({'pid':int(z[0]),'cpu_percent':float(z[2]),'rss_kib':int(z[3]),'cpu_time':z[4]})
    row={'wall_seconds':time.monotonic()-start,'worker_count':len(ws),'total_cpu_percent':sum(x['cpu_percent'] for x in ws),'total_rss_kib':sum(x['rss_kib'] for x in ws),'workers':ws}
   except Exception as ex:row={'wall_seconds':time.monotonic()-start,'monitor_error':repr(ex)}
   f.write(json.dumps(row,separators=(',',':'))+'\n');f.flush();delay=5
if __name__=='__main__':
 if sys.argv[1:]==['--register']:reg();sys.exit()
 assert sys.argv[1:]==['--run'];r=json.loads((D/'registration.json').read_bytes());assert r['code_sha256']==sha(Path(__file__).read_bytes());start=time.monotonic();stop=threading.Event();th=threading.Thread(target=monitor,args=(stop,os.getpid()),daemon=True);th.start();print(json.dumps({'event':'start','pid':os.getpid(),'workers':15,'cases':15,'registration_sha256':sha((D/'registration.json').read_bytes()),'monitor':str(D/'monitor.jsonl')}),flush=True)
 try:
  with concurrent.futures.ProcessPoolExecutor(max_workers=15) as ex:results=list(ex.map(task,r['cases']))
 finally:stop.set();th.join(timeout=10)
 out={'batch':'R046','registration_sha256':sha((D/'registration.json').read_bytes()),'wall_seconds':time.monotonic()-start,'results':results};(D/'summary.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({'event':'complete','cases':len(results),'wall_seconds':out['wall_seconds']}),flush=True)
