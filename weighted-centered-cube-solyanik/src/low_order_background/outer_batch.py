#!/usr/bin/env python3
"""15-process outer-halo checks for finite homogeneous singular mixtures.

H_alpha(E) is Fraction-exact and Lebesgue-only.  Background b is positive
piecewise constant (1<=b<=B); v is a finite set of separated |x-z|^(tau-1)
terms.  The output checks common-halo convexity and local scaled stability.
"""
import os
for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS'):
    os.environ[key]='1'
import argparse,hashlib,json,math,os as _os,platform,sys,time
from concurrent.futures import ProcessPoolExecutor,as_completed
from fractions import Fraction as F
from pathlib import Path

ROOT=Path(__file__).resolve().parent; CHECK=ROOT/'checkpoints'
PARAMETER_VERSION='outer-homogeneous-v3-bg5-cacheguard-20261105'

def merge(I):
 out=[]
 for a,b in sorted((a,b) for a,b in I if b>a):
  if out and a<=out[-1][1]:out[-1]=(out[-1][0],max(out[-1][1],b))
  else:out.append((a,b))
 return out

def halo_fraction(E,alpha):
 # Local side-effect-free reuse of round20 final Fraction halo / old round013 verifier.
 ends=sum(([a,b] for a,b in E),[]);ell=sum(b-a for a,b in E);lo=ends[0]-ell/(2*alpha);hi=ends[-1]+ell/(2*alpha);pieces=list(E)
 for e in ends:
  cuts=sorted({lo,hi,*((e+f)/2 for f in ends if lo<(e+f)/2<hi)})
  def G(x):
   a,b=min(e,2*x-e),max(e,2*x-e)
   return sum(max(F(0),min(b,d)-max(a,c)) for c,d in E)-2*alpha*abs(x-e)
  for a,b in zip(cuts,cuts[1:]):
   ga,gb=G(a),G(b)
   if ga>0 and gb>0:pieces.append((a,b))
   elif ga>0:pieces.append((a,a+(b-a)*ga/(ga-gb)))
   elif gb>0:pieces.append((a+(b-a)*(-ga)/(gb-ga),b))
 return merge(pieces)

def intersect(a,b,c,d):
 a,b=max(a,c),min(b,d);return (a,b) if b>a else None
def signedpow(x,t):return math.copysign(abs(x)**t,x) if x else 0.
def singular_mass(c,z,tau,I):return sum(c*(signedpow(float(b-z),tau)-signedpow(float(a-z),tau))/tau for a,b in I)
def bg_mass(steps,I):
 # b=1+sum h*1_(c,d), so it is globally positive even outside listed steps.
 total=sum(float(b-a) for a,b in I)
 for a,b in I:
  for c,d,h in steps:
   q=intersect(a,b,c,d)
   if q:total+=h*float(q[1]-q[0])
 return total
def component_masses(config,I):
 return [bg_mass(config['b'],I)]+[singular_mass(c,z,tau,I) for c,z,tau in config['v']]

def configs():
 baseE=[[(F(-6,5),F(-1,20)),(F(1,20),F(6,5))],[(F(-1),F(-1,10)),(F(1,15),F(9,10))],[(F(-7,8),F(-1,8)),(F(1,10),F(19,20))]]
 tau=(.15,.25,.4,.5,.75,.9);out=[]
 for k in range(15):
  z1,z2=(F(-3-k%3,5),F(2+k%4,5)) if k<11 else (F(-3-k%3,25),F(3+k%3,25))
  B=4+4*(k%4); b=[(F(-3),F(-1),1.),(F(-1),F(0),float(B-1)),(F(0),F(1),float(B//2)),(F(1),F(3),1.)]
  out.append({'id':k,'b':b,'v':[(.45+.05*(k%4),z1,tau[k%6]),(.60+.06*(k%3),z2,tau[(2*k+1)%6])],'E':baseE[k%3]})
 return out

def config_signature(c):
 payload={'id':c['id'],'E':[[str(a),str(b)] for a,b in c['E']],
          'b':[[str(a),str(b),repr(h)] for a,b,h in c['b']],
          'v':[[repr(q),str(z),repr(tau)] for q,z,tau in c['v']],
          'alphas':['1/4','2/5','3/4','9/10','99/100'],'local_alpha':'9/10'}
 encoded=json.dumps(payload,sort_keys=True,separators=(',',':'))
 return hashlib.sha256(encoded.encode()).hexdigest()

def local(config,j):
 c,z,tau=config['v'][j];other=config['v'][1-j][1];sep=abs(z-other);alpha=F(9,10);rows=[]
 for k in range(4,13):
  eps=F(1,2**k)
  if 5*eps>=sep:continue
  E=[(z-F(3,2)*eps,z-F(1,20)*eps),(z+F(1,20)*eps,z+F(3,2)*eps)];H=halo_fraction(E,alpha)
  den=component_masses(config,E);num=component_masses(config,H);pure=num[j+1]/den[j+1];mixed=sum(num)/sum(den);scale=float(eps)**(1-tau)
  rows.append({'k':k,'epsilon':float(eps),'pure':pure,'mixed':mixed,'difference':mixed-pure,'scaled_abs_difference':abs(mixed-pure)/scale,'halo_intervals':len(H)})
 return {'site':str(z),'tau':tau,'separation':str(sep),'rows':rows}

def task(c,barrier=None):
 if barrier is not None:barrier.wait(timeout=30)
 f=CHECK/f"config_{c['id']:02d}.json"
 sig=config_signature(c)
 if f.exists():
  cached=json.loads(f.read_text())
  if cached.get('parameter_version')==PARAMETER_VERSION and cached.get('config_signature')==sig:
   cached['cache_state']='reused';cached['checkpoint_worker_pid']=cached.get('worker_pid');cached['worker_pid']=_os.getpid()
   return cached
 start=time.time(); rows=[]
 Bmax=1+max(x[2] for x in c['b'])
 for alpha in (F(1,4),F(2,5),F(3,4),F(9,10),F(99,100)):
  H=halo_fraction(c['E'],alpha);de=component_masses(c,c['E']);dh=component_masses(c,H);ratios=[x/y for x,y in zip(dh,de)];mix=sum(dh)/sum(de)
  bg_bound=2*Bmax*float(1-alpha)/float(alpha)
  rows.append({'alpha':str(alpha),'halo_intervals':len(H),'ratios':ratios,'mixture':mix,'convex_upper':max(ratios),'convex_pass':mix<=max(ratios)+2e-12,
               'background_ratio_minus_one':ratios[0]-1,'background_BG5_bound':bg_bound,'background_BG5_pass':ratios[0]-1<=bg_bound+2e-12})
 loc=[local(c,0),local(c,1)]
 finite=all(math.isfinite(v) for x in rows for v in x['ratios']+[x['mixture']]) and all(math.isfinite(v) for q in loc for x in q['rows'] for v in (x['pure'],x['mixed'],x['difference']))
 valid=finite and all(x['convex_pass'] and x['background_BG5_pass'] for x in rows)
 r={'id':c['id'],'Bmax':Bmax,'rows':rows,'local':loc,'worker_pid':_os.getpid(),'checkpoint_worker_pid':_os.getpid(),
    'seconds':time.time()-start,'parameter_version':PARAMETER_VERSION,'config_signature':sig,'cache_state':'computed','run_validation':valid}
 tmp=f.with_suffix('.tmp');tmp.write_text(json.dumps(r,indent=2,allow_nan=False));tmp.replace(f);return r

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--workers',type=int,default=15);a=ap.parse_args();CHECK.mkdir(exist_ok=True);cs=configs();start=time.time();done=[];fail=[]
 from multiprocessing import Manager
 manager=Manager();barrier=manager.Barrier(15);print(json.dumps({'event':'start','tasks':15,'workers':min(15,a.workers),'platform':platform.platform(),'python':sys.version}),flush=True)
 with ProcessPoolExecutor(max_workers=min(15,a.workers)) as ex:
  fs={ex.submit(task,c,barrier if a.workers>=15 else None):c['id'] for c in cs}
  for f in as_completed(fs):
   try:
    r=f.result();done.append(r)
    if not r.get('run_validation',False):fail.append({'id':r['id'],'error':'run_validation=false'})
    print(json.dumps({'event':'done','id':r['id'],'pid':r['worker_pid'],'cache':r.get('cache_state'),'valid':r.get('run_validation',False)}),flush=True)
   except Exception as e:
    fail.append({'id':fs[f],'error':repr(e)})
    print(json.dumps({'event':'fail','id':fs[f],'error':repr(e)}),flush=True)
 manager.shutdown();out={'finished':len(done),'failures':fail,'seconds':time.time()-start,'records':sorted(done,key=lambda x:x['id']),'audit_full_pool':a.workers>=15,'parameter_version':PARAMETER_VERSION}
 (ROOT/'outer_summary.json').write_text(json.dumps(out,indent=2,allow_nan=False));print(json.dumps({'event':'finish','finished':len(done),'failures':len(fail),'seconds':out['seconds']}),flush=True)
 if fail: raise SystemExit(1)
if __name__=='__main__':main()
