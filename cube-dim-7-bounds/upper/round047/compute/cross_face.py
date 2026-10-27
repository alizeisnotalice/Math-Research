#!/usr/bin/env python3
"""R047 full pinned-face pair moment, continuous free-face variables."""
import os,sys,time,json,hashlib,resource,concurrent.futures,threading,subprocess
from pathlib import Path
from fractions import Fraction as Q
os.environ.update({k:'1' for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS')})
import mpmath as mp
D=Path(__file__).resolve().parent;sha=lambda b:hashlib.sha256(b).hexdigest()
class Source:
 def __init__(self,R,rho,g):
  self.R=R;self.rho=rho;self.g=g;self.a=R-g;self.ell=1-g/(2*R);self.e0=self.E(0)
 def E(self,m):
  z=self.g-m
  return mp.mpf(0) if z<=0 else -2*self.R*mp.log1p(-z/(2*self.R))-z
 def delta(self,w):
  w=abs(w);m=max(mp.mpf(0),w-self.a);return (self.E(m)+self.e0-m)/(2*self.rho)
 def H(self,z):return (4*self.R**2-2*self.R*z)*mp.log1p(-z/(2*self.R))+2*self.R*z-z*z/2
 def primitive(self,w):
  a=self.a;g=self.g;R=self.R
  if w<=a:return w*self.e0/self.rho
  m=w-a
  if m<=g:return (2*a*self.e0+m*self.e0-m*m/2+self.H(g)-self.H(g-m))/(2*self.rho)
  return (2*a*self.e0+g*self.e0-g*g/2+self.H(g)+(m-g)*self.e0-(m*m-g*g)/2)/(2*self.rho)
 def s(self,w):return (1+(1 if abs(w)<self.a else 0))/self.ell

def aggreg(n,pin,free):
 # pin=(J1,J2,s1,s2); free=(mu,a=E[s1 J2],b=E[s2 J1],c=E[s1 s2]); supports toy enumeration.
 J1,J2,s1,s2=pin;mu,a,b,c=free;k=n-1
 return (s1*s2*mu**k+k*(s1*J2*b+s2*J1*a+J1*J2*c)*mu**(k-1)+k*(k-1)*J1*J2*a*b*mu**(k-2))/4

def micro():
 mp.mp.dps=80;n=4;r=mp.mpf('1.1');s1=Source(mp.mpf('1.02'),mp.mpf(8),mp.mpf('.03'));s2=Source(mp.mpf('1.04'),mp.mpf(9),mp.mpf('.025'));vals=[mp.mpf(0),r/2];prob=mp.mpf('.5')
 pin=(1+s1.delta(r),1+s2.delta(r),s1.s(r),s2.s(r));free=(sum((1+s1.delta(w))*(1+s2.delta(w)) for w in vals)*prob,sum(s1.s(w)*(1+s2.delta(w)) for w in vals)*prob,sum(s2.s(w)*(1+s1.delta(w)) for w in vals)*prob,sum(s1.s(w)*s2.s(w) for w in vals)*prob)
 import itertools
 exact=mp.mpf(0)
 for ys in itertools.product(vals,repeat=n-1):
  ws=[r,*ys];S=[]
  for src in (s1,s2):S.append(sum(src.s(w)*mp.fprod(1+src.delta(ws[k]) for k in range(n) if k!=i) for i,w in enumerate(ws))/2)
  exact+=S[0]*S[1]/2**(n-1)
 got=aggreg(n,pin,free);assert abs(got-exact)<mp.mpf('1e-70')
 out={'status':'PASS','n':n,'enumerated_configurations':2**(n-1),'enumerated':mp.nstr(exact,40),'aggregate':mp.nstr(got,40),'abs_difference':mp.nstr(abs(got-exact),10),'scope':'independent finite discrete algebra check; not a continuous theorem'};(D/'micro.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out),flush=True)

def register():
 # Complete finite seed lists: deterministic Cartesian product, no pseudorandom samples.
 grid={'u_over_d':[str(Q(1,4)+Q(i,30)) for i in range(16)],'theta1':[str(-1+Q(2*i,15)) for i in range(16)],'theta2':[str(-1+Q(2*i,15)) for i in range(16)]};cases=[{'m':m,'d':str(Q(1,2**k))} for m in (32,64,128,256,512) for k in (8,12,16)]
 out={'batch':'R047','model':'complete source pinned-face cross-radius second moment; numerical diagnostics, not a Cn theorem','code_sha256':sha(Path(__file__).read_bytes()),'micro_sha256':sha((D/'micro.json').read_bytes()),'cases':cases,'grid':grid,'cartesian_triples_per_case':4096,'workers':15,'precision_dps':80,'integration':'exact piecewise first-moment primitives; mpmath piecewise 1D quad of shared delta1*delta2; cache swapped radius pairs','limit':'4096 triples per task; stop that task on first numerical >4096 or invalid geometry; no new scans','monitor':'first 0.05s then 5s','registered_unix_time':time.time()};(D/'registration.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({'event':'registered','sha256':sha((D/'registration.json').read_bytes()),'tasks':15,'triples':61440,'code_sha256':out['code_sha256']}),flush=True)

def task(case,grid):
 mp.mp.dps=80;start=time.monotonic();cpu=time.process_time();m=case['m'];n=mp.mpf(2)**m;dQ=Q(case['d']);d=mp.mpf(dQ.numerator)/dQ.denominator;x=m*mp.log(2)/2;t=n/x;zeta=n**-3;h=mp.mpf(2)**-40/(x*d);assert 0<h<d
 us=[Q(v) for v in grid['u_over_d']];ths=[Q(v) for v in grid['theta1']];rows=[];calls=0;quad_parts=0;largest=mp.mpf(0);where=None;failure=None
 def qmp(fr):return mp.mpf(fr.numerator)/fr.denominator
 def mu_pow(mu_delta,offset):
  # n−offset may round to n at 80 digits; this form retains the correction.
  v=mp.log1p(mu_delta);return mp.exp(n*v-offset*v)
 for ui,uf in enumerate(us):
  r=1+d*qmp(uf);sources=[]
  for theta in ths:
   q=r-d+qmp(theta)*h;rho=t*q;R=q+d-h;g=d-h+zeta*rho
   if not (q>0 and rho>=4*R and 0<g<R/2 and r>R-g):
    failure={'u_over_d':str(uf),'theta':str(theta),'reason':'invalid geometry'};break
   sources.append(Source(R,rho,g))
  if failure:break
  ccache={};means=[s.primitive(r)/r for s in sources];invr=1/r;invn=1/n
  for i in range(16):
   s1=sources[i]
   for j in range(16):
    s2=sources[j];key=tuple(sorted((i,j)))
    if key not in ccache:
     cuts=sorted(set([mp.mpf(0),min(r,s1.a),min(r,s1.R),min(r,s2.a),min(r,s2.R),r]));parts=[(l,u) for l,u in zip(cuts,cuts[1:]) if u>l]
     cross=mp.fsum(mp.quad(lambda w:s1.delta(w)*s2.delta(w),[l,u]) for l,u in parts)*invr;calls+=1;quad_parts+=len(parts);ccache[key]=cross
    cross=ccache[key];mu_delta=means[i]+means[j]+cross
    a12=(1+s1.a*invr+means[j]+s2.primitive(s1.a)*invr)/s1.ell
    b21=(1+s2.a*invr+means[i]+s1.primitive(s2.a)*invr)/s2.ell
    c12=(1+(s1.a+s2.a+min(s1.a,s2.a))*invr)/(s1.ell*s2.ell)
    Jp1=1+s1.delta(r);Jp2=1+s2.delta(r);sp1=s1.s(r);sp2=s2.s(r)
    # exact full pin-pinned + pin-free + same-free + distinct-free terms, scaled by n^2
    value=(sp1*sp2*invn**2*mu_pow(mu_delta,1)+invn*(1-invn)*(sp1*Jp2*b21+sp2*Jp1*a12+Jp1*Jp2*c12)*mu_pow(mu_delta,2)+(1-invn)*(1-2*invn)*Jp1*Jp2*a12*b21*mu_pow(mu_delta,3))/4
    if value>largest:largest=value;where={'u_over_d':str(uf),'theta1':str(ths[i]),'theta2':str(ths[j])}
    rows.append({'u_over_d':str(uf),'theta1':str(ths[i]),'theta2':str(ths[j]),'S_cross_over_n2':mp.nstr(value,32)})
    if not mp.isfinite(value) or value<0 or value>4096:
     failure={'u_over_d':str(uf),'theta1':str(ths[i]),'theta2':str(ths[j]),'value':mp.nstr(value,50),'reason':'candidate counterexample or numerical error'};break
   if failure:break
  if failure:break
  if ui%4==3:print(json.dumps({'event':'progress','m':m,'d':str(dQ),'u_done':ui+1,'triples':len(rows),'pid':os.getpid()}),flush=True)
 out={'m':m,'n':'2^'+str(m),'d':str(dQ),'status':'counterexample_or_error' if failure else 'completed','failure':failure,'triples_completed':len(rows),'cross_integral_calls':calls,'quadrature_parts':quad_parts,'source_states':len(us)*len(ths),'max_S_cross_over_n2':mp.nstr(largest,40),'max_parameters':where,'results':rows,'precision_dps':80,'pid':os.getpid(),'wall_seconds':time.monotonic()-start,'cpu_seconds':time.process_time()-cpu,'peak_rss_bytes':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,'scope':'one-dimensional quadrature diagnostic; not an interval proof or operator claim'}
 name='case_m'+str(m)+'_d'+str(dQ.denominator)+'.json';(D/name).write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({'event':'case','m':m,'d':str(dQ),'status':out['status'],'triples':len(rows),'cross_calls':calls,'max':mp.nstr(largest,12),'pid':out['pid'],'wall':out['wall_seconds']}),flush=True);return {'file':name,'sha256':sha((D/name).read_bytes()),'status':out['status'],'pid':out['pid'],'wall_seconds':out['wall_seconds'],'cpu_seconds':out['cpu_seconds'],'peak_rss_bytes':out['peak_rss_bytes'],'cross_integral_calls':calls,'triples_completed':len(rows)}

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
 if sys.argv[1:]==['--micro']:micro();sys.exit()
 if sys.argv[1:]==['--register']:register();sys.exit()
 assert sys.argv[1:]==['--run'];r=json.loads((D/'registration.json').read_bytes());assert r['code_sha256']==sha(Path(__file__).read_bytes()) and r['micro_sha256']==sha((D/'micro.json').read_bytes());start=time.monotonic();stop=threading.Event();th=threading.Thread(target=monitor,args=(stop,os.getpid()),daemon=True);th.start();print(json.dumps({'event':'start','pid':os.getpid(),'workers':15,'registration_sha256':sha((D/'registration.json').read_bytes()),'monitor':str(D/'monitor.jsonl')}),flush=True)
 try:
  with concurrent.futures.ProcessPoolExecutor(max_workers=15) as ex:results=list(ex.map(task,r['cases'],[r['grid']]*15))
 finally:stop.set();th.join(timeout=10)
 out={'batch':'R047','registration_sha256':sha((D/'registration.json').read_bytes()),'wall_seconds':time.monotonic()-start,'results':results};(D/'summary.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({'event':'complete','cases':len(results),'wall_seconds':out['wall_seconds']}),flush=True)
