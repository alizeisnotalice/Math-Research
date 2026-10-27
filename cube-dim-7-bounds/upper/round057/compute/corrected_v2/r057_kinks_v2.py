#!/usr/bin/env python3
"""R057 finite inverse-source derivative/kink diagnostic, no operator theorem."""
import os,sys,time,json,hashlib,resource,concurrent.futures,threading,subprocess
from pathlib import Path
from fractions import Fraction as Q
os.environ.update({k:'1' for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS')})
import mpmath as mp
D=Path(__file__).resolve().parent;sha=lambda b:hashlib.sha256(b).hexdigest();DS=[Q(1,16)+(Q(1,3)-Q(1,16))*Q(i,14) for i in range(15)];FACT=(Q(4),Q(16),Q(256),Q(2**16),Q(2**32),Q(2**64),Q(2**128),Q(2**256));HS=(Q(1,128),Q(1,32),Q(1,8));GAMMAS=(Q(0),Q(1,4),Q(1,2));JS=(24,48,96,160)
def mq(v):return mp.mpf(v.numerator)/v.denominator

def calc(q,d,h,gamma,t,w,side=None):
 rho=t*q;R=q+d-h;g=d-h+gamma*q;a=q*(1-gamma);aa=abs(w);m=max(Q(0),aa-a)
 Rm=mq(R);gm=mq(g);mm=mq(m);qm=mq(q);gmprime=mq(gamma);amprime=1-gmprime
 def E(mv,mdot):
  v=gm-mv
  if v<=0:return mp.mpf(0),mp.mpf(0)
  vdot=gmprime-mdot;y=v/(2*Rm);f=-2*Rm*mp.log1p(-y)-v;df=-2*mp.log1p(-y)+(vdot-v/Rm)/(1-y)-vdot;return f,df
 E0,E0p=E(mp.mpf(0),mp.mpf(0));
 if m>0:mdot=-amprime
 elif aa==a and side=='left':mdot=-amprime
 else:mdot=mp.mpf(0)
 Em,Emp=E(mm,mdot);B=E0+Em-mm;Bp=E0p+Emp-mdot;delta=B/(2*mq(rho));scaled=(Bp-B/qm)/(2*(1+delta))
 return delta,scaled

def direct(q,d,h,gamma,t,w):
 rho=t*q;R=q+d-h;g=d-h+gamma*q;a=q*(1-gamma);rhoM=mq(rho);RM=mq(R);gm=mq(g);am=mq(a);wm=mq(w)
 cuts=sorted(set([-rhoM,rhoM,-rhoM+gm,rhoM-gm,max(-rhoM,min(rhoM,wm-rhoM-am)),max(-rhoM,min(rhoM,wm+rhoM+am))]))
 def f(z):
  ell=min(mp.mpf(1),(RM+rhoM+am-abs(z))/(2*RM));return (1/ell if abs(z-wm)<rhoM+am else 0)-1
 return mp.fsum(mp.quad(f,[a0,b0]) for a0,b0 in zip(cuts,cuts[1:]) if b0>a0)/(2*rhoM)

def register():
 grid={'d':[str(x) for x in DS],'q_rule':'1-3d/4+i*d/30, i=0..15','rho_over_r':[str(x) for x in FACT],'h_over_d':[str(x) for x in HS],'gamma_over_h':[str(x) for x in GAMMAS],'kink_perturb_bits':list(JS),'w':'0, ±r, ±(r-2d), ±a, ±(a±2^-j); a=q(1-gamma)','r':'q+d','t':'rho/q fixed during derivative','gamma':'zeta*t fixed during derivative'}
 out={'batch':'R057','objective':'high-precision derivative/kink finite counterexample diagnostics, not Cn theorem','code_sha256':sha(Path(__file__).read_bytes()),'grid':grid,'precision_decimal_digits':160,'independent_crosscheck_decimal_digits':256,'predecessor':'../registration.json and ../summary.json; first run failed precision checks; original retained','tasks':240,'workers':15,'states_per_task':8*3*3*23,'max_states':240*8*3*3*23,'independent_checks':'each chunk: rho/r in {4,2^256}, h=d/8, gamma=h/2, w in {a-2^-96,a,a+2^-96,r}; direct piecewise quad and one/two-sided finite differences','finite_difference_step':'2^-200','tolerances':'rho*closed-vs-direct < 1e-95; rho*derivative difference <1e-35','stop':'all finite states or first numerical violation/exception stops affected chunk with exact parameters','monitor':'0.5 sec then 5 sec','registered_unix_time':time.time()};(D/'registration.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({'event':'registered','sha256':sha((D/'registration.json').read_bytes()),'states_max':out['max_states'],'tasks':240}),flush=True)

def task(arg):
 mp.mp.dps=160;start=time.monotonic();cpu=time.process_time();di,qi=arg;d=DS[di];q=1-3*d/4+d*Q(qi,30);r=q+d;states=0;directchecks=[];maximum=mp.mpf(0);maxdata=None;worstJ=mp.mpf(0);fail=None
 try:
  for fac in FACT:
   rho=fac*r;t=rho/q
   for hf in HS:
    h=d*hf
    for gf in GAMMAS:
     gamma=h*gf;a=q*(1-gamma);assert rho>=4*r and h/d<=Q(1,8) and gamma<=h/2 and 0<a<r
     points=[('zero',Q(0)),('plus_r',r),('minus_r',-r),('plus_boundary',r-2*d),('minus_boundary',-(r-2*d)),('kink_plus',a),('kink_minus',-a)]
     for j in JS:
      e=Q(1,2**j)
      for z in (a-e,a+e):points.extend([('near_'+str(j)+'_plus',z),('near_'+str(j)+'_minus',-z)])
     for label,w in points:
      if abs(w)>r:raise AssertionError('target point outside r')
      sides=('left','right') if abs(w)==a else (None,)
      for side in sides:
       delta,v=calc(q,d,h,gamma,t,w,side);indicator=int(abs(w)>r-2*d);den=mq(d*d+indicator);score=abs(v)/den
       if score>maximum:maximum=score;maxdata={'rho_over_r':str(fac),'h_over_d':str(hf),'gamma_over_h':str(gf),'w':str(w),'label':label,'side':side,'score':mp.nstr(score,30)}
       rhoD=mq(rho)*delta;worstJ=max(worstJ,abs(rhoD));assert rhoD>=-3*mq(d)-mp.mpf('1e-120') and rhoD<=2*mq(d*d)+mp.mpf('1e-120') and score<=64+mp.mpf('1e-100')
      states+=1
     if hf==HS[-1] and gf==GAMMAS[-1] and fac in (FACT[0],FACT[-1]):
      with mp.workdps(256):
       ep=Q(1,2**96);dq=Q(1,2**200)
       for w in (a-ep,a,a+ep,r):
        d0=direct(q,d,h,gamma,t,w);ref,_=calc(q,d,h,gamma,t,w);assert abs(mq(rho)*(d0-ref))<mp.mpf('1e-95')
        base=mp.log1p(d0);forward=(mp.log1p(direct(q+dq,d,h,gamma,t,w))-base)/mq(dq);backward=(base-mp.log1p(direct(q-dq,d,h,gamma,t,w)))/mq(dq)
        left=calc(q,d,h,gamma,t,w,'left' if abs(w)==a else None)[1];right=calc(q,d,h,gamma,t,w,'right' if abs(w)==a else None)[1]
        err=max(abs(mq(rho)*backward-left),abs(mq(rho)*forward-right));assert err<mp.mpf('1e-35')
        directchecks.append({'rho_over_r':str(fac),'w':str(w),'at_kink':abs(w)==a,'scaled_derivative_max_error':mp.nstr(err,10)})
 except Exception as exc:fail={'d':str(d),'q':str(q),'r':str(r),'rho_over_r':str(fac),'h_over_d':str(hf),'gamma_over_h':str(gf),'label':label if 'label' in locals() else None,'w':str(w) if 'w' in locals() else None,'error':repr(exc)}
 out={'d':str(d),'q':str(q),'q_index':qi,'status':'exception_or_candidate_counterexample' if fail else 'completed','failure':fail,'states_completed':states,'max_scaled_derivative_ratio':mp.nstr(maximum,40),'max_parameters':maxdata,'max_abs_rho_Jminus1':mp.nstr(worstJ,30),'direct_checks':directchecks,'direct_quadrature_calls':len(directchecks)*3,'pid':os.getpid(),'wall_seconds':time.monotonic()-start,'cpu_seconds':time.process_time()-cpu,'peak_rss_bytes':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,'scope':'finite 160-digit numerical diagnostics; kink uses separate one-sided derivatives; no rigorous interval proof'}
 name='chunk_d'+str(di)+'_q'+str(qi)+'.json';(D/name).write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({'event':'chunk','d_index':di,'q_index':qi,'pid':out['pid'],'status':out['status'],'states':states,'max':mp.nstr(maximum,8),'direct':len(directchecks),'wall':out['wall_seconds']}),flush=True);return {'file':name,'sha256':sha((D/name).read_bytes()),'status':out['status'],'states':states,'direct_checks':len(directchecks),'pid':out['pid'],'wall_seconds':out['wall_seconds'],'cpu_seconds':out['cpu_seconds'],'peak_rss_bytes':out['peak_rss_bytes']}

def monitor(stop,parent):
 start=time.monotonic();delay=.5
 with (D/('monitor_'+sys.argv[1][2:]+'.jsonl')).open('w') as f:
  while not stop.wait(delay):
   try:
    raw=subprocess.check_output(['ps','-axo','pid=,ppid=,%cpu=,rss=,time=,command='],text=True);ws=[]
    for line in raw.splitlines():
     z=line.strip().split(None,5)
     if len(z)==6 and int(z[1])==parent and 'spawn_main' in z[5]:ws.append({'pid':int(z[0]),'cpu_percent':float(z[2]),'rss_kib':int(z[3]),'cpu_time':z[4]})
    row={'wall_seconds':time.monotonic()-start,'worker_count':len(ws),'total_cpu_percent':sum(v['cpu_percent'] for v in ws),'total_rss_kib':sum(v['rss_kib'] for v in ws),'workers':ws}
   except Exception as exc:row={'wall_seconds':time.monotonic()-start,'monitor_error':repr(exc)}
   f.write(json.dumps(row,separators=(',',':'))+'\n');f.flush();delay=5
if __name__=='__main__':
 if sys.argv[1:]==['--register']:register();sys.exit()
 assert sys.argv[1:] in (['--pilot'],['--run']);mode=sys.argv[1][2:];r=json.loads((D/'registration.json').read_bytes());assert r['code_sha256']==sha(Path(__file__).read_bytes());tasks=([(0,j) for j in range(16)] if mode=='pilot' else [(i,j) for i in range(1,15) for j in range(16)]);start=time.monotonic();stop=threading.Event();th=threading.Thread(target=monitor,args=(stop,os.getpid()),daemon=True);th.start();print(json.dumps({'event':'start','pid':os.getpid(),'workers':15,'tasks':len(tasks),'registration_sha256':sha((D/'registration.json').read_bytes()),'monitor':str(D/'monitor.jsonl')}),flush=True)
 try:
  with concurrent.futures.ProcessPoolExecutor(max_workers=15) as ex:results=list(ex.map(task,tasks))
 finally:stop.set();th.join(timeout=10)
 out={'batch':'R057','registration_sha256':sha((D/'registration.json').read_bytes()),'wall_seconds':time.monotonic()-start,'results':results};(D/('summary_'+mode+'.json')).write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({'event':'complete','chunks':len(results),'wall_seconds':out['wall_seconds'],'states':sum(x['states'] for x in results),'direct_checks':sum(x['direct_checks'] for x in results)}),flush=True)
