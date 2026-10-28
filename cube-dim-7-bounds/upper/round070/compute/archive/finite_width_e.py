import os,sys,time,json,hashlib,resource,concurrent.futures,threading,subprocess
from fractions import Fraction as F
from pathlib import Path
os.environ.update({k:'1' for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS')})
import mpmath as mp
P=Path(__file__).resolve().parent;H=lambda p:hashlib.sha256(p.read_bytes()).hexdigest();S=F(3,2);Ls=(4096,8192);Hs=(256,512);TH=(F(-1,2),F(0),F(1,2),F(1));mp.mp.dps=90
nodes={}
def v(x):return mp.mpf(x.numerator)/x.denominator
def rule(N):
 if N not in nodes:
  X,W=mp.gauss_quadrature(N,'legendre');nodes[N]=(list(X),list(W))
 return nodes[N]
def calc(i,L,Hden,theta):
 q0=F(9,10)+F(1,50)*F(2*i+1,512);r=S*q0;h=F(1,Hden*L);q=q0+theta*h/S;eps=h/(2*S*q0);R=S*q-h;a=(1-eps)*q;g=R-a;assert R>0 and eps*q<=h/S and g>0
 n=mp.exp(L);x=mp.mpf(9)*L/10;rho=n*v(q)/x;rr=v(r);R0=v(R);a0=v(a);g0=v(g);Q0=-mp.log1p(-g0/(2*R0));half=v(h);cuts=sorted(set([-rr,-R0,-a0,a0,R0,rr]));cuts=[z for z in cuts if -rr<=z<=rr];
 def evalf(w):
  delta=max(mp.mpf(0),abs(w)-a0);z=g0-delta;Q1=-mp.log1p(-z/(2*R0)) if z>0 else mp.mpf(0);f2=R0*(Q0**3+Q1**3)/3;f4=R0*(Q0**5+Q1**5)/5
  def E(m):
   b=g0-m;return -2*R0*mp.log1p(-b/(2*R0))-b if b>0 else mp.mpf(0)
  Jcorr=(E(mp.mpf(0))+E(delta)-delta)/(2*rho);return f2,f4,Jcorr
 def integ(N):
  xs,ws=rule(N);s1=mp.mpf(0);minvar=mp.inf;maxcorr=mp.mpf(0);pieces=[]
  for l,u in zip(cuts,cuts[1:]):
   mid=(l+u)/2;wid=(u-l)/2
   for z,wgt in zip(xs,ws):
    f2,f4,delta=evalf(mid+wid*z);weight=wid*wgt;s1+=weight*f2;pieces.append((weight,f2));inverse_correction=delta/(1+delta);cv=f4-f4*inverse_correction-f2*f2/(rho*(1+delta)**2);minvar=min(minvar,cv);maxcorr=max(maxcorr,abs(inverse_correction))
  mean=s1/(2*rr);variance=sum((weight*(f2-mean)**2 for weight,f2 in pieces),mp.mpf(0))/(2*rr);return mean,variance,minvar,maxcorr
 m32,var32,_,_=integ(32);m64,var,mincv,maxcorr=integ(64);err=max(abs(m64-m32),abs(var-var32));assert mincv>=0 and var>=0
 mean=(n/rho)*m64;variance=n/(rho*rho)*var;cantelli={str(t):mp.nstr(variance/(variance+(v(t)*mean)**2),20) for t in (F(1,100),F(1,20),F(1,4))};baseline=None
 if theta==1:
  pin=1-g0/(2*R0);Q=-mp.log(pin);baseline=abs(m64-2*R0*(2-pin*(Q*Q+2*Q+2)))
 return {'L':L,'h_den':Hden,'theta':str(theta),'q0':str(q0),'q':str(q),'r':str(r),'rho':'exp(L)*q/((9/10)*L)','target_mean_nH2':mp.nstr(mean,35),'target_variance_sum':mp.nstr(variance,15),'min_scaled_conditional_source_variance':mp.nstr(mincv,25),'quad32_64_maxdiff':mp.nstr(err,12),'max_J_inverse_correction':mp.nstr(maxcorr,20),'baseline_scaled_error':None if baseline is None else mp.nstr(baseline,12),'cantelli_upper':cantelli}
def register():
 out={'batch':'R070 exp(L) finite-width variance and Cantelli diagnostic','source_sha256':H(Path(__file__)),'scope':'new finite-width conditional model; numerical 32/64 high-precision Gaussian quadrature; no Cn theorem','q0':'9/10+(1/50)*(2*i+1)/512, i=0..255','s':'3/2','L':list(Ls),'h':'1/(Hden*L), Hden 256 or 512','qsource':'q0+theta*h/s, theta=-1/2,0,1/2,1','r':'s*q0','epsilon':'h/(2*s*q0)','rho':'exp(L)*qsource/((9/10)*L)','target_w':'uniform [-r,r], integrating kink-split H2 and square','observables':'mean and variance of H2; normalized source conditional variance from H4/J-(H2/J)^2; Cantelli at relative gaps 1%,5%,25%','states':4096,'tasks':256,'workers':14,'precision':'90 decimal digits; carry rho*(J-1)=B/2 and keep 1/J correction separately; final sum may round subprecision correction; centered variance quadrature','stop':'finite tasks; any exception saved in chunk; no theoretical certificate','registered_unix':time.time()};(P/'width_e_registration.json').write_text(json.dumps(out,indent=2)+'\n');print(H(P/'width_e_registration.json'),flush=True)
def task(i):
 start=time.monotonic();cpu=time.process_time();rows=[];err=None
 try:
  for L in Ls:
   for k in Hs:
    for z in TH:rows.append(calc(i,L,k,z))
 except Exception as e:err={'L':L,'h_den':k,'theta':str(z),'error':repr(e)}
 out={'i':i,'status':'exception_or_candidate_counterexample' if err else 'completed','failure':err,'rows':rows,'source_sha256':H(Path(__file__)),'pid':os.getpid(),'wall_seconds':time.monotonic()-start,'cpu_seconds':time.process_time()-cpu,'peak_rss_bytes':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss};f=P/f'width_e_chunk_{i:03}.json';f.write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({'i':i,'status':out['status'],'pid':out['pid'],'states':len(rows)}),flush=True);return {'file':f.name,'sha256':H(f),'status':out['status'],'states':len(rows),'cpu_seconds':out['cpu_seconds']}
def monitor(stop,parent):
 start=time.monotonic();delay=.5
 with (P/'width_e_monitor.jsonl').open('w') as f:
  while not stop.wait(delay):
   try:
    rows=subprocess.check_output(['ps','-axo','pid=,ppid=,%cpu=,rss=,time=,command='],text=True);ws=[]
    for line in rows.splitlines():
     z=line.strip().split(None,5)
     if len(z)==6 and int(z[1])==parent and 'spawn_main' in z[5]:ws.append({'pid':int(z[0]),'cpu_percent':float(z[2]),'rss_kib':int(z[3]),'cpu_time':z[4]})
    row={'wall_seconds':time.monotonic()-start,'worker_count':len(ws),'total_cpu_percent':sum(x['cpu_percent'] for x in ws),'total_rss_kib':sum(x['rss_kib'] for x in ws),'workers':ws}
   except Exception as e:row={'error':repr(e)}
   f.write(json.dumps(row)+'\n');f.flush();delay=5
if __name__=='__main__':
 if sys.argv[1:]==['--register']:register();sys.exit()
 assert sys.argv[1:]==['--run'];r=json.loads((P/'width_e_registration.json').read_text());assert r['source_sha256']==H(Path(__file__));start=time.monotonic();stop=threading.Event();th=threading.Thread(target=monitor,args=(stop,os.getpid()),daemon=True);th.start();print(json.dumps({'event':'start','pid':os.getpid(),'workers':14,'tasks':256,'registration_sha256':H(P/'width_e_registration.json')}),flush=True)
 try:
  with concurrent.futures.ProcessPoolExecutor(max_workers=14) as ex:rs=list(ex.map(task,range(256)))
 finally:stop.set();th.join(timeout=10)
 out={'registration_sha256':H(P/'width_e_registration.json'),'wall_seconds':time.monotonic()-start,'receipts':rs};(P/'width_e_summary.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({'event':'complete','wall_seconds':out['wall_seconds'],'states':sum(x['states'] for x in rs),'failures':sum(x['status']!='completed' for x in rs)}),flush=True)
