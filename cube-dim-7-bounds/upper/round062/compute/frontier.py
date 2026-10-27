import os,json,time,hashlib,heapq,resource,concurrent.futures,threading,subprocess,sys
from pathlib import Path
from fractions import Fraction as F
os.environ.update({k:'1' for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS')})
import mpmath as mp
P=Path(__file__).resolve().parent;sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest();Js=(16,32,64,128)
def m(x):return mp.mpf(x.numerator)/x.denominator
def register():
 reg={'batch':'R062 frontier pilot','scope':'finite layered ideal model only; no C_n theorem','code_sha256':sha(Path(__file__)),'d':'9/50+(23/150)*i/31, i=0..31','qhi':'9/10+(9/100)*k/31, k=0..31','q_global_min':'1-d','J':list(Js),'tuple_count':4096,'band_evaluations':32*32*sum(Js),'tasks':128,'workers':15,'decimal_digits':80,'objective':'sum_j log((qhi_j+d)/(qlo_j+d))*qlo_j/qhi_j - 2 max_j mu_j','mu':'(qlo_j+d)/qlo_j * log(1-d/[2(qlo_j+d)])^2','eligible':'all mu_j < (1-qhi_j)/(2*qhi_j), qhi>qlo','stop':'finite all tuples, each task halts on exception with receipt','registered_unix':time.time()};(P/'frontier_registration.json').write_text(json.dumps(reg,indent=2)+'\n');print(sha(P/'frontier_registration.json'),flush=True)
def task(x):
 mp.mp.dps=80;i,J=x;d=F(9,50)+F(23,150)*F(i,31); start=time.monotonic();cpu=time.process_time();checked=bands=eligible=0;heap=[];error=None
 try:
  for k in range(32):
   hi=F(9,10)+F(9,100)*F(k,31);q0=1-d;assert hi>q0;s=mp.mpf(0);maxmu=mp.mpf(0);gap=mp.inf
   for j in range(J):
    ql=q0+(hi-q0)*F(j,J);qu=q0+(hi-q0)*F(j+1,J);D,L,U=map(m,(d,ql,qu));R=L+D;phi=-mp.log1p(-D/(2*R));mu=R/L*phi**2;eta=(1-U)/(2*U);s+=mp.log((U+D)/R)*L/U;maxmu=max(maxmu,mu);gap=min(gap,eta-mu);bands+=1
   c=s-2*maxmu;checked+=1
   if gap>0:
    eligible+=1;rec={'i':i,'k':k,'d':str(d),'qhi':str(hi),'J':J,'c_model':mp.nstr(c,70),'min_eta_minus_mu':mp.nstr(gap,70),'max_mu':mp.nstr(maxmu,70)};entry=(c,k,rec)
    if len(heap)<32:heapq.heappush(heap,entry)
    elif entry[:2]>heap[0][:2]:heapq.heapreplace(heap,entry)
 except Exception as e:error={'error':repr(e),'i':i,'J':J,'k':k,'j':j}
 result={'i':i,'J':J,'status':'error' if error else 'completed','error':error,'tuples':checked,'bands':bands,'eligible':eligible,'top32':[x[2] for x in sorted(heap,reverse=True)],'pid':os.getpid(),'cpu_seconds':time.process_time()-cpu,'wall_seconds':time.monotonic()-start,'peak_rss_bytes':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss};p=P/f'frontier_{i:02}_{J}.json';p.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({'i':i,'J':J,'pid':result['pid'],'status':result['status'],'tuples':checked,'bands':bands}),flush=True);return {'file':p.name,'sha256':sha(p),'status':result['status'],'tuples':checked,'bands':bands,'eligible':eligible,'cpu_seconds':result['cpu_seconds']}
def monitor(stop,parent):
 start=time.monotonic();delay=.5
 with (P/'frontier_monitor.jsonl').open('w') as f:
  while not stop.wait(delay):
   try:
    raw=subprocess.check_output(['ps','-axo','pid=,ppid=,%cpu=,rss=,time=,command='],text=True);ws=[]
    for line in raw.splitlines():
     z=line.strip().split(None,5)
     if len(z)==6 and int(z[1])==parent and 'spawn_main' in z[5]:ws.append({'pid':int(z[0]),'cpu_percent':float(z[2]),'rss_kib':int(z[3]),'cpu_time':z[4]})
    row={'wall_seconds':time.monotonic()-start,'workers':ws,'worker_count':len(ws),'total_cpu_percent':sum(x['cpu_percent'] for x in ws),'total_rss_kib':sum(x['rss_kib'] for x in ws)}
   except Exception as e:row={'error':repr(e)}
   f.write(json.dumps(row)+'\n');f.flush();delay=5
if __name__=='__main__':
 if sys.argv[1:]==['--register']:register();sys.exit()
 assert sys.argv[1:]==['--run'];r=json.loads((P/'frontier_registration.json').read_text());assert r['code_sha256']==sha(Path(__file__));tasks=[(i,j) for i in range(32) for j in Js];stop=threading.Event();th=threading.Thread(target=monitor,args=(stop,os.getpid()),daemon=True);th.start();start=time.monotonic();print(json.dumps({'event':'start','pid':os.getpid(),'tasks':len(tasks),'workers':15,'registration_sha256':sha(P/'frontier_registration.json')}),flush=True)
 try:
  with concurrent.futures.ProcessPoolExecutor(max_workers=15) as ex:rs=list(ex.map(task,tasks))
 finally:stop.set();th.join(timeout=10)
 tops=[]
 for rr in rs:tops.extend(json.loads((P/rr['file']).read_text())['top32'])
 mp.mp.dps=80;tops.sort(key=lambda x:mp.mpf(x['c_model']),reverse=True);out={'registration_sha256':sha(P/'frontier_registration.json'),'receipts':rs,'top64':tops[:64],'wall_seconds':time.monotonic()-start};(P/'frontier_summary.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({'event':'complete','wall_seconds':out['wall_seconds'],'tuples':sum(z['tuples'] for z in rs),'bands':sum(z['bands'] for z in rs),'errors':sum(z['status']!='completed' for z in rs),'eligible':sum(z['eligible'] for z in rs),'top':tops[0] if tops else None}),flush=True)
