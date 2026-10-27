#!/usr/bin/env python3
"""R052 exact finite candidate enumeration; conditional scalar model only."""
import os,sys,time,json,hashlib,resource,concurrent.futures,threading,subprocess,heapq
from fractions import Fraction as Q
from pathlib import Path
os.environ.update({k:'1' for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS')})
if hasattr(sys,'set_int_max_str_digits'):sys.set_int_max_str_digits(0)
D=Path(__file__).resolve().parent;DEN=2**20;T=(1,2,4,8,16,32,64,128);SOURCE=Path('/Users/zhengzhihao/Desktop/T/cube-r030-soft/cube-dim-7-bounds/upper/round045/layered_source/PROOF.md');sha=lambda b:hashlib.sha256(b).hexdigest()
class Score:
 __slots__=('n','d')
 def __init__(self,n,d):self.n=n;self.d=d
 def __lt__(self,o):return self.n*o.d<o.n*self.d
 def __eq__(self,o):return self.n*o.d==o.n*self.d

def integer(A,k,t):
 # a=A/8, d=k/2^20, s=(4A+t)/16, rho=t/16. M=4 verified globally.
 if 4*A*k>=DEN:return None,'d_guard'
 pn=(4*A+t)*k*k;pd=16*DEN*DEN-pn
 if 2*pn>=16*DEN*DEN:return None,'sigma_guard'
 Fnum=512*3579*(A-24)-8200*A
 if Fnum<=0:return None,'positive_deficit_guard'
 num0=64449*k*Fnum;den0=131072*(4*DEN+3*k)*4198400*A
 return (num0*pd-pn*den0,den0*pd,A,k,t),None

def fraction(A,k,t):
 a=Q(A,8);d=Q(k,DEN);s=2*a+Q(t,16);P=Q(3579,4096);eps=Q(1,1024)
 assert a>3 and s>2*a and 0<d<=Q(1,16) and d<Q(1,32)/a and s*d*d<Q(1,2)
 ulo=P/(2*(1+eps));uhi=1/(P*(1-d/2)*(1-eps));mass=(1-Q(1,64))*d/(2*(1+3*d/4));M=max(Q(4),2*uhi/(1-3*d/4));term=(1-3/a)*ulo-2*eps;lim=mass*term/M;delta=(1-eps)*lim;sig=s*d*d;c=delta-sig/(1-sig)
 return {'a':str(a),'rho':str(Q(t,16)),'s':str(s),'d':str(d),'u_lo':str(ulo),'u_hi':str(uhi),'mass':str(mass),'M':str(M),'bracket':str(term),'delta_lim':str(lim),'delta':str(delta),'sigma':str(sig),'penalty':str(sig/(1-sig)),'c_raw':str(c),'c_raw_decimal':format(float(c),'.17g')}

def check():
 # M=4 throughout the registered superset 0<d<=1/256, verified at its worst endpoint.
 p=Q(3579,4096);eps=Q(1,1024);dmax=Q(4096,DEN);assert p*(1-dmax/2)*(1-eps)*(1-3*dmax/4)>Q(1,2)
 tests=[]
 for A in (32,64,128,256):
  for k in (1,63,511,1024,4096):
   for t in (1,16,128):
    v,reason=integer(A,k,t)
    if reason:continue
    q=fraction(A,k,t);assert q['M']=='4' and Q(q['c_raw'])==Q(v[0],v[1]);tests.append([A,k,t])
 out={'status':'PASS','M_equals_four_globally':'Plo*(1-dmax/2)*(1-eps)*(1-3dmax/4)>1/2 by exact fractions','dmax':str(dmax),'integer_vs_Fraction_cases':len(tests),'test_points':tests,'source_sha256':sha(Path(__file__).read_bytes())};(D/'equivalence.json').write_text(json.dumps(out,indent=2)+'\n');print('equivalence PASS',len(tests),'valid rational points',flush=True)

def register():
 check();proofhash=sha(SOURCE.read_bytes());assert proofhash=='22db93b275416b90e9f15a0118a88cd412fee0c2d0b20c4cce1ed32769f259ac'
 out={'batch':'R052','model_status':'conditional scalar candidate only; no operator-bound claim before R051 audit and full cost comparison','R045_source_sha256':proofhash,'R051_source_status':'proof in progress; no frozen theorem hash yet','code_sha256':sha(Path(__file__).read_bytes()),'equivalence_sha256':sha((D/'equivalence.json').read_bytes()),'A_inclusive':[32,256],'a':'A/8 = 4+k/8, k=0..224','rho_numerators_over16':list(T),'d':'k/2^20, k=1..4096','tuple_count':7372800,'guard_policy':'reject d>=1/(32a), sigma>=1/2, nonpositive bracket; count all guards; other declared guards hold on registered ranges','chunk':'one A and 512 consecutive k, all eight rho; 1800 chunks','workers':15,'thread_env':1,'local_top':32,'global_top':64,'stop':'all 7372800 tuples exactly once; no time cap, no unbounded grid','monitor':'0.5s initial, then 5s','registered_unix_time':time.time()};(D/'registration.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({'event':'registered','registration_sha256':sha((D/'registration.json').read_bytes()),'code_sha256':out['code_sha256'],'tuples':out['tuple_count']}),flush=True)

def task(chunk):
 A,k0,k1=chunk;start=time.monotonic();cpu=time.process_time();top=[];count=0;guard={'d_guard':0,'sigma_guard':0,'positive_deficit_guard':0};valid=0
 for k in range(k0,k1+1):
  for t in T:
   count+=1;v,reason=integer(A,k,t)
   if reason:guard[reason]+=1;continue
   valid+=1;item=(Score(v[0],v[1]),-A,-t,-k,v)
   if len(top)<32:heapq.heappush(top,item)
   elif top[0]<item:heapq.heapreplace(top,item)
 return {'A':A,'k0':k0,'k1':k1,'pid':os.getpid(),'count':count,'valid':valid,'guards':guard,'top':[x[-1] for x in top],'wall_seconds':time.monotonic()-start,'cpu_seconds':time.process_time()-cpu,'peak_rss_bytes':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}

def push(heap,v,limit):
 num,den,A,k,t=v;item=(Score(num,den),-A,-t,-k,v)
 if len(heap)<limit:heapq.heappush(heap,item)
 elif heap[0]<item:heapq.heapreplace(heap,item)

def monitor(stop,parent):
 start=time.monotonic();delay=.5
 with (D/'monitor.jsonl').open('w') as f:
  while not stop.wait(delay):
   try:
    raw=subprocess.check_output(['ps','-axo','pid=,ppid=,%cpu=,rss=,time=,command='],text=True);workers=[]
    for line in raw.splitlines():
     v=line.strip().split(None,5)
     if len(v)==6 and int(v[1])==parent and 'spawn_main' in v[5]:workers.append({'pid':int(v[0]),'cpu_percent':float(v[2]),'rss_kib':int(v[3]),'cpu_time':v[4]})
    row={'wall_seconds':time.monotonic()-start,'worker_count':len(workers),'cpu_percent':sum(x['cpu_percent'] for x in workers),'rss_kib':sum(x['rss_kib'] for x in workers),'workers':workers}
   except Exception as e:row={'wall_seconds':time.monotonic()-start,'error':repr(e)}
   f.write(json.dumps(row,separators=(',',':'))+'\n');f.flush();delay=5
if __name__=='__main__':
 if sys.argv[1:]==['--register']:register();sys.exit()
 assert sys.argv[1:]==['--run'];r=json.loads((D/'registration.json').read_bytes());assert r['code_sha256']==sha(Path(__file__).read_bytes()) and r['equivalence_sha256']==sha((D/'equivalence.json').read_bytes()) and r['R045_source_sha256']==sha(SOURCE.read_bytes());chunks=[(A,k,k+511) for A in range(32,257) for k in range(1,4097,512)];assert len(chunks)==1800
 start=time.monotonic();stop=threading.Event();th=threading.Thread(target=monitor,args=(stop,os.getpid()),daemon=True);th.start();print(json.dumps({'event':'start','pid':os.getpid(),'workers':15,'chunks':len(chunks),'registration_sha256':sha((D/'registration.json').read_bytes()),'monitor':str(D/'monitor.jsonl')}),flush=True)
 totals={'tuples':0,'valid':0,'d_guard':0,'sigma_guard':0,'positive_deficit_guard':0};globaltop=[];workers={};done=0
 try:
  with concurrent.futures.ProcessPoolExecutor(max_workers=15) as ex:
   futures=[ex.submit(task,c) for c in chunks]
   for fu in concurrent.futures.as_completed(futures):
    z=fu.result();done+=1;totals['tuples']+=z['count'];totals['valid']+=z['valid']
    for k,v in z['guards'].items():totals[k]+=v
    w=workers.setdefault(z['pid'],{'chunks':0,'tuples':0,'cpu_seconds':0.,'peak_rss_bytes':0,'top':[]});w['chunks']+=1;w['tuples']+=z['count'];w['cpu_seconds']+=z['cpu_seconds'];w['peak_rss_bytes']=max(w['peak_rss_bytes'],z['peak_rss_bytes'])
    for v in z['top']:push(globaltop,v,64);push(w['top'],v,32)
    if z['A']==32 and z['k0']==1:
     vv=max(z['top'],key=lambda t:Q(t[0],t[1]));print(json.dumps({'event':'warm','A':32,'k_range':[1,512],'tuples':z['count'],'valid':z['valid'],'best':fraction(vv[2],vv[3],vv[4]),'pid':z['pid']}),flush=True)
    if done%50==0:print(json.dumps({'event':'progress','chunks':done,'tuples':totals['tuples'],'wall':time.monotonic()-start}),flush=True)
 finally:stop.set();th.join(timeout=10)
 assert done==1800 and totals['tuples']==7372800 and totals['valid']+sum(totals[x] for x in ('d_guard','sigma_guard','positive_deficit_guard'))==totals['tuples']
 ranked=sorted([e[-1] for e in globaltop],key=lambda v:(-Q(v[0],v[1]),v[2],v[4],v[3]));outworkers={str(pid):{**{k:v for k,v in w.items() if k!='top'},'top':[fraction(x[-1][2],x[-1][3],x[-1][4]) for x in sorted(w['top'],reverse=True)]} for pid,w in workers.items()}
 out={'batch':'R052','model_status':r['model_status'],'registration_sha256':sha((D/'registration.json').read_bytes()),'totals':totals,'chunks_completed':done,'wall_seconds':time.monotonic()-start,'workers':outworkers,'top64':[fraction(v[2],v[3],v[4]) for v in ranked]};(D/'results.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({'event':'complete','tuples':totals,'wall_seconds':out['wall_seconds'],'best':out['top64'][0]}),flush=True)
