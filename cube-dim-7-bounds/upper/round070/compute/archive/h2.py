import os,sys,time,json,hashlib,resource,concurrent.futures,threading,subprocess,itertools
from pathlib import Path
from fractions import Fraction as F
os.environ.update({k:'1' for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS')})
import mpmath as mp
P=Path(__file__).resolve().parent;H=lambda p:hashlib.sha256(p.read_bytes()).hexdigest();SS=(F(6,5),F(7,5),F(3,2),F(9,5),F(2));HF=[F(1,2**x) for x in (24,48,80,96,128)];RF=[F(4),F(16),F(2**16),F(2**64),F(2**128),F(2**256)];EF=[F(0),F(1,4),F(1,2)];BITS=(24,64,128,192)
def v(x):return mp.mpf(x.numerator)/x.denominator
def h2(q,s,h,eps,t,w,side=None):
 rho=t*q;R=s*q-h;g=R-(1-eps)*q;a=(1-eps)*q;m=max(F(0),abs(w)-a);R0=v(R);g0=v(g);rp=v(s);ep=v(eps);gp=rp-1+ep
 def Q(x,xp):
  z=g0-x
  if z<=0:return mp.mpf(0),mp.mpf(0)
  return -mp.log1p(-z/(2*R0)),(gp-xp-z*rp/R0)/(2*R0-z)
 md=-(1-ep) if m>0 or (abs(w)==a and side=='left') else mp.mpf(0);q0,q0p=Q(mp.mpf(0),mp.mpf(0));q1,q1p=Q(v(m),md);T=q0**3+q1**3;A=R0*T/3;Ap=rp*T/3+R0*(q0*q0*q0p+q1*q1*q1p);return A/v(rho),Ap-A/v(q)
def direct(q,s,h,eps,t,w):
 rho=t*q;R=s*q-h;a=(1-eps)*q;g=R-a;rh=v(rho);rm=v(R);am=v(a);gm=v(g);wm=v(w);cuts=sorted(set([-rh,rh,-rh+gm,rh-gm,max(-rh,min(rh,wm-rh-am)),max(-rh,min(rh,wm+rh+am))]))
 def f(z):
  ell=min(mp.mpf(1),(rm+rh+am-abs(z))/(2*rm));return ((-mp.log(ell))**2/ell if abs(z-wm)<rh+am else 0)
 return mp.fsum(mp.quad(f,[x,y]) for x,y in zip(cuts,cuts[1:]) if y>x)/(2*rh)
def register():
 d={'batch':'R070 second target profile H2','scope':'finite numerical diagnostic, no Cn theorem','source_sha256':H(Path(__file__)),'s':[str(s) for s in SS],'q':'1/s+(97/100−1/s)*(2*i+1)/512, i=0..255','r':'s*q fixed during derivative','h_over_r':[str(h) for h in HF],'rho_over_r':[str(z) for z in RF],'epsilon':'factor*h/(s*q), factor 0,1/4,1/2','R':'s*q-h','a0':'(1-epsilon)*q','g':'R-a0','w':'0, ±r, ±a0, ±(a0±2^-j), j=24,64,128,192, ±(r-g), ±(r-g/2)','states_per_chunk':11250,'tasks':256,'pilot_indices':'0..13','production_indices':'14..255','total_states':2880000,'workers':14,'precision':'closed 160 decimal, cross-check 256 decimal (>512 bits)','crosscheck':'each chunk one of five s values rotating, h/r=2^-80, rho/r=2^256, eps factor 1/2; direct source integral and one-sided derivatives at kink and perturbations 2^-192; baseline mean and second moment/variance','threshold':'quadrature scaled error <1e-80; derivative scaled error <1e-25; baseline mean error <1e-100','stop':'pilot failure blocks production; in production exception or numerical disagreement stops affected chunk and records exact point; finite counts otherwise','registered_unix':time.time()};(P/'registration.json').write_text(json.dumps(d,indent=2)+'\n');print(H(P/'registration.json'),flush=True)
def task(i):
 mp.mp.dps=160;start=time.monotonic();cpu=time.process_time();states=0;checks=[];fail=None;maxval=mp.mpf(0);maxdata=None;maxratio=mp.mpf(0)
 try:
  for s in SS:
   q=1/s+(F(97,100)-1/s)*F(2*i+1,512);r=s*q
   for hf in HF:
    h=r*hf
    for rf in RF:
     rho=rf*r;t=rho/q
     for ef in EF:
      eps=ef*h/(s*q);a=(1-eps)*q;g=s*q-h-a;assert eps*q<=h/s and g>0
      points=[('zero',F(0)),('r+',r),('r-',-r),('kink+',a),('kink-',-a),('pin+',r-g),('pin-',-r+g),('midpin+',r-g/2),('midpin-',-r+g/2)]
      for j in BITS:
       for z in (a-F(1,2**j),a+F(1,2**j)):points.extend([('near+',z),('near-',-z)])
      for label,w in points:
       assert abs(w)<=r
       for side in (('left','right') if abs(w)==a else (None,)):
        value,der=h2(q,s,h,eps,t,w,side)
        if abs(der)>maxval:maxval=abs(der);maxdata={'s':str(s),'q':str(q),'h_over_r':str(hf),'rho_over_r':str(rf),'epsilon_factor':str(ef),'w':str(w),'label':label,'side':side,'scaled_derivative':mp.nstr(der,30)}
        assert value>=0
        R=v(s*q-h);G=v(g);M=v(max(F(0),abs(w)-a))
        def E(m):
         z=G-m
         return -2*R*mp.log1p(-z/(2*R))-z if z>0 else mp.mpf(0)
        B=E(mp.mpf(0))+E(M)-M;J=1+B/(2*v(rho));assert J>0;maxratio=max(maxratio,v(rho)*value/J)
       states+=1
      if s==SS[i%len(SS)] and hf==HF[2] and rf==RF[-1] and ef==EF[-1]:
       with mp.workdps(256):
        dq=F(1,2**240)
        for w in (a-F(1,2**192),a,a+F(1,2**192),r):
         z=direct(q,s,h,eps,t,w);ref,_=h2(q,s,h,eps,t,w);e0=abs(v(rho)*(z-ref));dl=(z-direct(q-dq,s,h,eps,t,w))/v(dq);dr=(direct(q+dq,s,h,eps,t,w)-z)/v(dq);left=h2(q,s,h,eps,t,w,'left' if abs(w)==a else None)[1];right=h2(q,s,h,eps,t,w,'right' if abs(w)==a else None)[1];e1=max(abs(v(rho)*dl-left),abs(v(rho)*dr-right));
         if e0>=mp.mpf('1e-80') or e1>=mp.mpf('1e-25'):raise ArithmeticError('independent precision crosscheck '+mp.nstr(e0,8)+' '+mp.nstr(e1,8))
         checks.append({'kind':'direct','w':str(w),'at_kink':abs(w)==a,'scaled_value_error':mp.nstr(e0,12),'scaled_derivative_error':mp.nstr(e1,12)})
        qb=(r+h)/s;ab=(1-eps)*qb;gb=r-ab;R=v(r);A=v(ab);G=v(gb);cuts=[-R,-A,A,R]
        def f(x):
         mm=max(mp.mpf(0),abs(x)-A);q0=-mp.log1p(-G/(2*R));z=G-mm;q1=-mp.log1p(-z/(2*R)) if z>0 else mp.mpf(0);return R*(q0**3+q1**3)/3
        mean=mp.fsum(mp.quad(f,[x,y]) for x,y in zip(cuts,cuts[1:]))/(2*R);m2=mp.fsum(mp.quad(lambda x:f(x)**2,[x,y]) for x,y in zip(cuts,cuts[1:]))/(2*R);pin=1-G/(2*R);Q=-mp.log(pin);expected=2*R*(2-pin*(Q**2+2*Q+2));err=abs(mean-expected);variance=m2-mean**2
        if err>mp.mpf('1e-100') or variance < -mp.mpf('1e-100'):raise ArithmeticError('baseline mean or variance check '+mp.nstr(err,8))
        checks.append({'kind':'baseline','q_baseline':str(qb),'s':str(s),'rho_scaled_mean_error':mp.nstr(err,12),'rho2_scaled_variance':mp.nstr(variance,25),'rho2_scaled_second_moment':mp.nstr(m2,25)})
 except Exception as ex:fail={'error':repr(ex),'i':i,'s':str(s),'q':str(q),'h_over_r':str(hf),'rho_over_r':str(rf),'epsilon_factor':str(ef),'w':str(w) if 'w' in locals() else None,'label':label if 'label' in locals() else None}
 out={'i':i,'status':'exception_or_candidate_counterexample' if fail else 'completed','failure':fail,'states':states,'max_abs_rho_derivative':mp.nstr(maxval,40),'max_parameters':maxdata,'max_rho_H2_over_J':mp.nstr(maxratio,40),'independent_checks':checks,'source_sha256':H(Path(__file__)),'pid':os.getpid(),'wall_seconds':time.monotonic()-start,'cpu_seconds':time.process_time()-cpu,'peak_rss_bytes':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss};f=P/f'chunk_{i:03}.json';f.write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({'i':i,'status':out['status'],'pid':out['pid'],'states':states,'checks':len(checks)}),flush=True);return {'file':f.name,'sha256':H(f),'status':out['status'],'states':states,'checks':len(checks),'cpu_seconds':out['cpu_seconds']}
def monitor(stop,parent,mode):
 start=time.monotonic();delay=.5
 with (P/(mode+'_monitor.jsonl')).open('w') as f:
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
 assert sys.argv[1:] in (['--pilot'],['--production']);mode=sys.argv[1][2:];reg=json.loads((P/'registration.json').read_text());assert reg['source_sha256']==H(Path(__file__))
 if mode=='production':
  pilot=json.loads((P/'pilot_summary.json').read_text());assert len(pilot['receipts'])==14 and all(x['status']=='completed' for x in pilot['receipts'])
 start=time.monotonic();stop=threading.Event();th=threading.Thread(target=monitor,args=(stop,os.getpid(),mode),daemon=True);th.start();indices=range(14) if mode=='pilot' else range(14,256);print(json.dumps({'event':'start','mode':mode,'pid':os.getpid(),'registration_sha256':H(P/'registration.json'),'workers':14,'tasks':len(indices)}),flush=True)
 try:
  with concurrent.futures.ProcessPoolExecutor(max_workers=14) as ex:rs=list(ex.map(task,indices))
 finally:stop.set();th.join(timeout=10)
 out={'mode':mode,'registration_sha256':H(P/'registration.json'),'wall_seconds':time.monotonic()-start,'receipts':rs};(P/(mode+'_summary.json')).write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({'event':'complete','mode':mode,'wall_seconds':out['wall_seconds'],'states':sum(x['states'] for x in rs),'checks':sum(x['checks'] for x in rs),'failures':sum(x['status']!='completed' for x in rs)}),flush=True)
