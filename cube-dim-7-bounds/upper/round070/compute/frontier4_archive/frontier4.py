import os,sys,time,json,math,heapq,hashlib,resource,concurrent.futures,threading,subprocess
from pathlib import Path
from fractions import Fraction as F
os.environ.update({k:'1' for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS')})
P=Path(__file__).resolve().parent;sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest();HDS=(256,384,512,768);eJ=2**-20;eB=2**-19;eN=eB;eS=eJ;eT=eB;eD=eJ;ea=1/65536;er=eJ;trim=2**-19;qmin=2/3
Lval=lambda i:F(4096)+F(4096*(2*i+1),40)
qval=lambda i:F(22,25)+F(1,20)*F(2*i+1,128)
sval=lambda i:F(7,100)+F(3,50)*F(2*i+1,128)
cval=lambda i:F(1,25)+F(1,100)*F(2*i+1,64)
def register():
 out={'batch':'R070 conditional exp-prefix full-cost frontier','scope':'conditional scalar proposal, remaining finite guards and operator theory unverified','source_sha256':sha(Path(__file__)),'L':'4096+4096*(2*i+1)/40,i=0..19','qmax':'22/25+(1/20)*(2*j+1)/128,j=0..63','sigma':'7/100+(3/50)*(2*k+1)/128,k=0..63','chi':'1/25+(1/100)*(2*m+1)/64,m=0..31','hden':list(HDS),'tuple_count':20*64*64*32*4,'tasks':1280,'workers':15,'x':'(1-sigma)*L','rho_dimension':'n=exp(L), L natural log','cost':'(1-eD)*mL*(2*p*(1-eRad)*Ulo-eS-eT)/(2*(1+eRad)*Uhi)+ln(1-sigma) - (exp(-L)+exp(-sigma L)+2 sqrt(L) exp(-sigma L/2))','retain':'gap^2/(.00531*x+gap^2), gap=chi*x-.183-2211*x/(50000*(1-eJ))','cap':'ln(qmax/qmin)*rhi/h * exp(-((1-exp(-L))*eta-chi)*x)','tau':'rhi*sqrt(L)/h*exp(-(sigma/2-chi*(1-sigma))*L)','screening':'double float, rejected margins not rigorous; high precision and positive S64 prefix recheck top candidates required','stop':'finite counts; each chunk records its counts and top32; on exception affected chunk stops','registered_unix':time.time()};(P/'frontier4_registration.json').write_text(json.dumps(out,indent=2)+'\n');print(sha(P/'frontier4_registration.json'),flush=True)
def task(z):
 il,iq=z;L=float(Lval(il));q=float(qval(iq));rhi=1.5*q-trim;rlo=1+trim;eta=(1-q)/(2*q);Csrc=math.log(q/qmin);m=(1-ea)*math.log(rhi/rlo);start=time.monotonic();cpu=time.process_time();counts={'domain':0,'payment':0,'eligible':0};heap=[];error=None
 try:
  for isg in range(64):
   sigma=float(sval(isg));x=(1-sigma)*L;cost0=math.log1p(-sigma)-(math.exp(-L)+math.exp(-sigma*L)+2*math.sqrt(L)*math.exp(-sigma*L/2))
   for ic in range(32):
    chi=float(cval(ic));gap=chi*x-.183-2211*x/(50000*(1-eJ));cs=(eta-chi)*x;ct=(sigma/2-chi*(1-sigma))*L
    for D in HDS:
     if min(gap,cs,ct)<=0 or 2/(D*L)>=trim:counts['domain']+=1;continue
     cap=Csrc*rhi*D*L*math.exp(-cs);tau=rhi*math.sqrt(L)*D*L*math.exp(-ct)
     if cap>eS or tau>eT:counts['payment']+=1;continue
     p=gap*gap/(.00531*x+gap*gap);P0=1-eB-6/D;ul=P0*(1-eN)/(1+eJ);uh=(1+eN)/(P0*(1-eJ));b=2*p*(1-er)*ul-eS-eT;gamma=(1-eD)*m*b/(2*(1+er)*uh);cost=gamma+cost0;counts['eligible']+=1
     row={'il':il,'iq':iq,'isigma':isg,'ichi':ic,'hden':D,'L':str(Lval(il)),'qmax':str(qval(iq)),'sigma':str(sval(isg)),'chi':str(cval(ic)),'cost_float':cost,'gamma_float':gamma,'retention_float':p,'cap_float':cap,'tau_float':tau};z=(cost,isg,ic,D,row)
     if len(heap)<32:heapq.heappush(heap,z)
     elif z[:4]>heap[0][:4]:heapq.heapreplace(heap,z)
 except Exception as ex:error={'error':repr(ex),'isigma':isg,'ichi':ic,'hden':D}
 out={'il':il,'iq':iq,'status':'error' if error else 'completed','failure':error,'counts':counts,'top32':[z[4] for z in sorted(heap,reverse=True)],'source_sha256':sha(Path(__file__)),'pid':os.getpid(),'cpu_seconds':time.process_time()-cpu,'wall_seconds':time.monotonic()-start,'peak_rss_bytes':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss};f=P/f'frontier4_{il:02}_{iq:02}.json';f.write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({'il':il,'iq':iq,'pid':out['pid'],'status':out['status'],'counts':counts}),flush=True);return {'file':f.name,'sha256':sha(f),'status':out['status'],'counts':counts,'cpu_seconds':out['cpu_seconds']}
def monitor(stop,parent):
 start=time.monotonic();delay=.5
 with (P/'frontier4_monitor.jsonl').open('w') as f:
  while not stop.wait(delay):
   try:
    text=subprocess.check_output(['ps','-axo','pid=,ppid=,%cpu=,rss=,time=,command='],text=True);ws=[]
    for line in text.splitlines():
     z=line.strip().split(None,5)
     if len(z)==6 and int(z[1])==parent and 'spawn_main' in z[5]:ws.append({'pid':int(z[0]),'cpu_percent':float(z[2]),'rss_kib':int(z[3]),'cpu_time':z[4]})
    row={'wall_seconds':time.monotonic()-start,'worker_count':len(ws),'total_cpu_percent':sum(x['cpu_percent'] for x in ws),'total_rss_kib':sum(x['rss_kib'] for x in ws),'workers':ws}
   except Exception as e:row={'error':repr(e)}
   f.write(json.dumps(row)+'\n');f.flush();delay=5
if __name__=='__main__':
 if sys.argv[1:]==['--register']:register();sys.exit()
 assert sys.argv[1:]==['--run'];r=json.loads((P/'frontier4_registration.json').read_text());assert r['source_sha256']==sha(Path(__file__));start=time.monotonic();stop=threading.Event();th=threading.Thread(target=monitor,args=(stop,os.getpid()),daemon=True);th.start();print(json.dumps({'event':'start','pid':os.getpid(),'workers':15,'tasks':1280,'registration_sha256':sha(P/'frontier4_registration.json')}),flush=True)
 try:
  with concurrent.futures.ProcessPoolExecutor(max_workers=15) as ex:rs=list(ex.map(task,[(i,j) for i in range(20) for j in range(64)]))
 finally:stop.set();th.join(timeout=10)
 tops=[]
 for x in rs:tops.extend(json.loads((P/x['file']).read_text())['top32'])
 tops.sort(key=lambda z:z['cost_float'],reverse=True);out={'registration_sha256':sha(P/'frontier4_registration.json'),'wall_seconds':time.monotonic()-start,'receipts':rs,'top64_float_proposals':tops[:64]};(P/'frontier4_summary.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({'event':'complete','wall_seconds':out['wall_seconds'],'counts':{k:sum(x['counts'][k] for x in rs) for k in ('domain','payment','eligible')},'failures':sum(x['status']!='completed' for x in rs),'best':tops[0] if tops else None}),flush=True)
