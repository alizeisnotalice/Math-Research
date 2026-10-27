import os,sys,time,json,hashlib,heapq,resource,concurrent.futures,threading,subprocess
from pathlib import Path
from fractions import Fraction as F
os.environ.update({k:'1' for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS')})
import mpmath as mp
P=Path(__file__).resolve().parent
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
def vals(i,j,k):
 return F(1,5)+F(2,15)*F(i,127),F(1,16)+F(3,8)*F(j,127),F(9,16)+F(27,64)*F(k,127)
def asmp(v):return mp.mpf(v.numerator)/v.denominator
def register():
 x={'batch':'R058','identity':'conditional ideal sourceband model; no operator bound','code_sha256':sha(Path(__file__)),'d':'1/5+(2/15)*i/127','a':'1/16+(3/8)*j/127','b':'9/16+(27/64)*k/127','i_j_k':'0..127 inclusive','tasks':128,'tuples_per_task':16384,'total_tuples':2097152,'precision_decimal':80,'workers':15,'qualifying':'a<b, mu<eta and c_model>0','objective':'log((qmax+d)/(qmin+d))*qmin/qmax - 2*mu','stop':'finite all tuples or first exception in affected chunk; persist receipts even on exception','registration_time_unix':time.time(),'output':'per chunk top32, count, source hash and timing; merged top64; no per-tuple storage'}
 (P/'registration.json').write_text(json.dumps(x,indent=2)+'\n');print('registration sha256',sha(P/'registration.json'),flush=True)
def task(i):
 mp.mp.dps=80; start=time.monotonic(); cpu=time.process_time(); heap=[];checked=eligible=positive=0;error=None
 try:
  for k in range(128):
   df,_,bf=vals(i,0,k);d=asmp(df);b=asmp(bf);qmin=1-b*d;Rmin=qmin+d;phi=-mp.log1p(-d/(2*Rmin));mu=(Rmin/qmin)*phi**2
   for j in range(128):
    _,af,_=vals(i,j,k);a=asmp(af);qmax=1-a*d;eta=a*d/(2*qmax);obj=mp.log((qmax+d)/Rmin)*qmin/qmax-2*mu;checked+=1
    if mu<eta:
     eligible+=1
     if obj>0:
      positive+=1;rec={'i':i,'j':j,'k':k,'d':str(df),'a':str(af),'b':str(bf),'c_model':mp.nstr(obj,75),'mu':mp.nstr(mu,75),'eta':mp.nstr(eta,75),'eta_minus_mu':mp.nstr(eta-mu,75)};entry=(obj,j,k,rec)
      if len(heap)<32:heapq.heappush(heap,entry)
      elif entry[:3]>heap[0][:3]:heapq.heapreplace(heap,entry)
 except Exception as e:error={'i':i,'j':j,'k':k,'exception':repr(e)}
 out={'index':i,'status':'error' if error else 'completed','error':error,'checked':checked,'eligible':eligible,'positive':positive,'top32':[x[3] for x in sorted(heap,reverse=True)],'pid':os.getpid(),'cpu_seconds':time.process_time()-cpu,'wall_seconds':time.monotonic()-start,'peak_rss_bytes':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,'code_sha256':sha(Path(__file__))}
 path=P/f'chunk_{i:03}.json';path.write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({'chunk':i,'status':out['status'],'pid':out['pid'],'checked':checked,'wall':out['wall_seconds']}),flush=True);return {'file':path.name,'sha256':sha(path),'status':out['status'],'checked':checked,'eligible':eligible,'positive':positive,'cpu_seconds':out['cpu_seconds']}
def monitor(stop,parent):
 start=time.monotonic();delay=.5
 with (P/'monitor.jsonl').open('w') as f:
  while not stop.wait(delay):
   try:
    rows=subprocess.check_output(['ps','-axo','pid=,ppid=,%cpu=,rss=,time=,command='],text=True); ws=[]
    for line in rows.splitlines():
     z=line.strip().split(None,5)
     if len(z)==6 and int(z[1])==parent and 'spawn_main' in z[5]:ws.append({'pid':int(z[0]),'cpu_percent':float(z[2]),'rss_kib':int(z[3]),'cpu_time':z[4]})
    row={'wall_seconds':time.monotonic()-start,'worker_count':len(ws),'total_cpu_percent':sum(x['cpu_percent'] for x in ws),'total_rss_kib':sum(x['rss_kib'] for x in ws),'workers':ws}
   except Exception as e:row={'monitor_error':repr(e)}
   f.write(json.dumps(row)+'\n');f.flush();delay=5
if __name__=='__main__':
 if sys.argv[1:]==['--register']:register();sys.exit()
 assert sys.argv[1:]==['--run'];reg=json.loads((P/'registration.json').read_text());assert reg['code_sha256']==sha(Path(__file__));start=time.monotonic();stop=threading.Event();th=threading.Thread(target=monitor,args=(stop,os.getpid()),daemon=True);th.start();print(json.dumps({'event':'start','pid':os.getpid(),'registration_sha256':sha(P/'registration.json'),'workers':15,'tasks':128}),flush=True)
 try:
  with concurrent.futures.ProcessPoolExecutor(max_workers=15) as e: rs=list(e.map(task,range(128)))
 finally:stop.set();th.join(timeout=10)
 mp.mp.dps=80;top=[]
 for z in rs:
  for x in json.loads((P/z['file']).read_text())['top32']:top.append(x)
 top.sort(key=lambda x:mp.mpf(x['c_model']),reverse=True)
 out={'batch':'R058','registration_sha256':sha(P/'registration.json'),'wall_seconds':time.monotonic()-start,'receipts':rs,'top64':top[:64]};(P/'summary.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({'event':'complete','wall_seconds':out['wall_seconds'],'tasks':len(rs),'checked':sum(z['checked'] for z in rs),'failed':sum(z['status']!='completed' for z in rs),'best':top[0] if top else None}),flush=True)
