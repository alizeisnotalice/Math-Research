import os,sys,time,json,hashlib,resource,concurrent.futures,threading,subprocess
from pathlib import Path
from fractions import Fraction as F
os.environ.update({k:'1' for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS')})
import mpmath as mp
P=Path(__file__).resolve().parent;H=lambda p:hashlib.sha256(p.read_bytes()).hexdigest();S=F(3,2);HF=[F(1,2**x) for x in (24,48,96,128)];RF=[F(4),F(16),F(2**16),F(2**64),F(2**128),F(2**256)];EF=[F(0),F(1,4),F(1,2)];BITS=(24,64,128,192)
def v(x):return mp.mpf(x.numerator)/x.denominator
def calc(q,h,eps,t,w,side=None):
 rho=t*q;R=S*q-h;g=R-(1-eps)*q;a=(1-eps)*q;w0=abs(w);m=max(F(0),w0-a);r0=v(R);g0=v(g);q0=v(q);ep=v(eps);rp=v(S)
 def e(x,xp):
  z=g0-x
  if z<=0:return mp.mpf(0),mp.mpf(0)
  zp=rp-(1-ep)-xp;y=z/(2*r0);return -2*r0*mp.log1p(-y)-z,-2*rp*mp.log1p(-y)+(zp-z*rp/r0)/(1-y)-zp
 md=-(1-ep) if m>0 or (w0==a and side=='left') else mp.mpf(0);E0,E0p=e(mp.mpf(0),mp.mpf(0));Em,Emp=e(v(m),md);B=E0+Em-v(m);Bp=E0p+Emp-md;delta=B/(2*v(rho));return delta,(Bp-B/q0)/(2*(1+delta))
def direct(q,h,eps,t,w):
 rho=t*q;R=S*q-h;a=(1-eps)*q;g=R-a;rh=v(rho);rm=v(R);am=v(a);gm=v(g);wm=v(w);cuts=sorted(set([-rh,rh,-rh+gm,rh-gm,max(-rh,min(rh,wm-rh-am)),max(-rh,min(rh,wm+rh+am))]))
 def f(z):
  ell=min(mp.mpf(1),(rm+rh+am-abs(z))/(2*rm));return (1/ell if abs(z-wm)<rh+am else 0)-1
 return mp.fsum(mp.quad(f,[x,y]) for x,y in zip(cuts,cuts[1:]) if y>x)/(2*rh)
def register():
 x={'batch':'R065 multiplicative capture kink','code_sha256':H(Path(__file__)),'q0':'2/3 + (909/1000-2/3)*(2*i+1)/256, i=0..127','s':'3/2','r':'s*q0 fixed during derivative','h_over_r':[str(x) for x in HF],'rho_over_r':[str(x) for x in RF],'epsilon':'factor*h/(s*q0), factors 0,1/4,1/2; fixed during derivative','R':'s*q-h','a0':'(1-epsilon)*q','g':'R-a0; derivative 1/2+epsilon','t':'rho/q0 fixed','points':'0, ±r, ±a0, ±(a0±2^-j) j=24,64,128,192, ±(r−g), ±(r−g/2)','states_per_chunk':4*6*3*25,'tasks':128,'max_states':128*4*6*3*25,'workers':15,'precision':'stable formula 160 decimal digits; direct quadrature/difference at 256 decimal digits (>512 bits)','direct_checks':'per chunk at h/r=2^-24 and 2^-128, rho/r=2^256, epsilon=half upper bound; w=a0−2^-192,a0,a0+2^-192,r, with difference step 2^-240','baseline':'at qb=(r+h)/s check uniform w mean J=1 via piecewise closed-form quadrature for same 2 parameter combinations per chunk','stop':'finite all states or first inequality violation or exception in affected chunk; exceptions explicitly retained','threshold':'abs(rho*dq log J)<=256; direct rho-scaled errors <1e-80 and derivative errors <1e-25','registered_unix':time.time()};(P/'derivative_registration.json').write_text(json.dumps(x,indent=2)+'\n');print(H(P/'derivative_registration.json'),flush=True)
def task(i):
 mp.mp.dps=160;start=time.monotonic();cpu=time.process_time();q=F(2,3)+(F(909,1000)-F(2,3))*F(2*i+1,256);r=S*q;states=0;maxv=mp.mpf(0);maxdata=None;checks=[];fail=None
 try:
  for hf in HF:
   h=r*hf
   for rf in RF:
    rho=rf*r;t=rho/q
    for ef in EF:
     eps=ef*h/(S*q);a=(1-eps)*q;g=S*q-h-a;assert eps*q<=h/S and g>0
     points=[('zero',F(0)),('r+',r),('r-',-r),('kink+',a),('kink-',-a),('pin+',r-g),('pin-',-r+g),('midpin+',r-g/2),('midpin-',-r+g/2)]
     for j in BITS:
      for z in (a-F(1,2**j),a+F(1,2**j)):points.extend([('near+',z),('near-',-z)])
     for label,w in points:
      assert abs(w)<=r
      for side in (('left','right') if abs(w)==a else (None,)):
       delta,z=calc(q,h,eps,t,w,side)
       if abs(z)>maxv:maxv=abs(z);maxdata={'h_over_r':str(hf),'rho_over_r':str(rf),'epsilon_factor':str(ef),'label':label,'w':str(w),'side':side,'scaled_derivative':mp.nstr(z,40)}
       if abs(z)>256:raise ArithmeticError('derivative bound violated')
      states+=1
     if hf in (HF[0],HF[-1]) and rf==RF[-1] and ef==EF[-1]:
      with mp.workdps(256):
       dq=F(1,2**240)
       for w in (a-F(1,2**192),a,a+F(1,2**192),r):
        z=direct(q,h,eps,t,w);ref,_=calc(q,h,eps,t,w);e0=abs(v(rho)*(z-ref));base=mp.log1p(z);dl=(base-mp.log1p(direct(q-dq,h,eps,t,w)))/v(dq);dr=(mp.log1p(direct(q+dq,h,eps,t,w))-base)/v(dq);left=calc(q,h,eps,t,w,'left' if abs(w)==a else None)[1];right=calc(q,h,eps,t,w,'right' if abs(w)==a else None)[1];e1=max(abs(v(rho)*dl-left),abs(v(rho)*dr-right));
        if e0>=mp.mpf('1e-80') or e1>=mp.mpf('1e-25'):raise ArithmeticError('independent numerical check failed: '+mp.nstr(e0,8)+' '+mp.nstr(e1,8))
        checks.append({'w':str(w),'at_kink':abs(w)==a,'h_over_r':str(hf),'direct_scaled_error':mp.nstr(e0,12),'derivative_scaled_error':mp.nstr(e1,12)})
       qb=(r+h)/S;ab=(1-eps)*qb;gb=r-ab;cuts=[-v(r),-v(ab),v(ab),v(r)]
       def ff(x):
        # integrate the stable formula at fixed qb without converting quadrature nodes to rational
        mm=max(mp.mpf(0),abs(x)-v(ab));R=v(r);g0=v(gb)
        def ee(z):
         y=g0-z
         return -2*R*mp.log1p(-y/(2*R))-y if y>0 else mp.mpf(0)
        return (ee(mp.mpf(0))+ee(mm)-mm)/2
       avg=mp.fsum(mp.quad(ff,[x,y]) for x,y in zip(cuts,cuts[1:]))/(2*v(r))
       if abs(avg)>mp.mpf('1e-100'):raise ArithmeticError('baseline E J check failed '+mp.nstr(avg,12))
       checks.append({'baseline_q':str(qb),'h_over_r':str(hf),'rho_scaled_mean_delta':mp.nstr(avg,12)})
 except Exception as ex:fail={'error':repr(ex),'i':i,'q0':str(q),'h_over_r':str(hf),'rho_over_r':str(rf),'epsilon_factor':str(ef),'w':str(w) if 'w' in locals() else None,'label':label if 'label' in locals() else None}
 out={'i':i,'q0':str(q),'status':'exception_or_candidate_counterexample' if fail else 'completed','failure':fail,'states':states,'max_abs_scaled_derivative':mp.nstr(maxv,40),'max_parameters':maxdata,'independent_checks':checks,'pid':os.getpid(),'wall_seconds':time.monotonic()-start,'cpu_seconds':time.process_time()-cpu,'peak_rss_bytes':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss};f=P/f'derivative_{i:03}.json';f.write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({'i':i,'status':out['status'],'pid':out['pid'],'states':states,'checks':len(checks)}),flush=True);return {'file':f.name,'sha256':H(f),'status':out['status'],'states':states,'checks':len(checks),'cpu_seconds':out['cpu_seconds']}
def monitor(stop,parent):
 start=time.monotonic();delay=.5
 with (P/'derivative_monitor.jsonl').open('w') as f:
  while not stop.wait(delay):
   try:
    text=subprocess.check_output(['ps','-axo','pid=,ppid=,%cpu=,rss=,time=,command='],text=True);ws=[]
    for line in text.splitlines():
     z=line.strip().split(None,5)
     if len(z)==6 and int(z[1])==parent and 'spawn_main' in z[5]:ws.append({'pid':int(z[0]),'cpu_percent':float(z[2]),'rss_kib':int(z[3]),'cpu_time':z[4]})
    row={'wall_seconds':time.monotonic()-start,'workers':ws,'worker_count':len(ws),'total_cpu_percent':sum(x['cpu_percent'] for x in ws),'total_rss_kib':sum(x['rss_kib'] for x in ws)}
   except Exception as ex:row={'error':repr(ex)}
   f.write(json.dumps(row)+'\n');f.flush();delay=5
if __name__=='__main__':
 if sys.argv[1:]==['--register']:register();sys.exit()
 assert sys.argv[1:]==['--run'];reg=json.loads((P/'derivative_registration.json').read_text());assert reg['code_sha256']==H(Path(__file__));start=time.monotonic();stop=threading.Event();th=threading.Thread(target=monitor,args=(stop,os.getpid()),daemon=True);th.start();print(json.dumps({'event':'start','pid':os.getpid(),'registration_sha256':H(P/'derivative_registration.json'),'workers':15,'tasks':128}),flush=True)
 try:
  with concurrent.futures.ProcessPoolExecutor(max_workers=15) as ex:rs=list(ex.map(task,range(128)))
 finally:stop.set();th.join(timeout=10)
 out={'registration_sha256':H(P/'derivative_registration.json'),'wall_seconds':time.monotonic()-start,'receipts':rs};(P/'derivative_summary.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({'event':'complete','wall_seconds':out['wall_seconds'],'states':sum(x['states'] for x in rs),'checks':sum(x['checks'] for x in rs),'failures':sum(x['status']!='completed' for x in rs)}),flush=True)
