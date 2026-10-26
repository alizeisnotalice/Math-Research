#!/usr/bin/env python3
"""Finite registered zeroth-DD positive-band engine; no global frontier."""
from pathlib import Path
from fractions import Fraction as F
from itertools import product
from math import prod
import argparse,gzip,hashlib,json,multiprocessing as mp,os,resource,signal,sys,time,traceback
import flint
from flint import arb,ctx,fmpq
if hasattr(sys,'set_int_max_str_digits'):sys.set_int_max_str_digits(0)
HERE=Path(__file__).resolve().parent

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def atomic(p,o):
 t=p.with_suffix('.tmp');t.write_text(json.dumps(o,sort_keys=True,separators=(',',':'))+'\n');os.replace(t,p)
def rq(x):
 a=F(str(x));return fmpq(a.numerator,a.denominator)
def ep(x,upper=False):return str((x.upper() if upper else x.lower()).fmpq())
def enc(x):return {'lower':ep(x),'upper':ep(x,True)}
def upper(x):return F(ep(x,True))
def quadratic(u):return sum(x*x for x in u)/4-sum(u)**2/36

def geometry(row):
 c=list(map(F,row['candidates'][row['selected']]['point']))[1:];indices=row['indices'];lo=[F(i,64)-x for i,x in zip(indices,c)];hi=[F(i+1,64)-x for i,x in zip(indices,c)]
 assert all(l<0<h for l,h in zip(lo,hi))
 mu=[(l+h)/2 for l,h in zip(lo,hi)];var=[(h-l)**2/12 for l,h in zip(lo,hi)]
 eta=quadratic(mu)+F(2,9)*sum(var);B2=max(quadratic(u) for u in product(*zip(lo,hi)))
 return {'lo':list(map(str,lo)),'hi':list(map(str,hi)),'mu':list(map(str,mu)),'eta':str(eta),'B2':str(B2),'h_inf':str(max(max(abs(l),abs(h)) for l,h in zip(lo,hi))),'lambda':str(max(map(abs,mu))),'volume':str(F(1,64)**4)}

def bounds(rows,g):
 B2,h,eta,lam=map(F,(g['B2'],g['h_inf'],g['eta'],g['lambda']));L=arb(0);U=arb(0);M=arb(0);bands=[]
 for i,row in enumerate(rows):
  a=F(row['t']);cur=row['H'];nxt=rows[i+1]['H'] if i+1<len(rows) else {'lower':'0','upper':'0'}
  ml=F(cur['lower'])-F(nxt['upper']);mu=F(cur['upper'])-F(nxt['lower']);assert 0<ml<=mu,('mass_not_positive',i,str(ml),str(mu))
  exp=arb(rq(h/a)).exp();s=arb(rq(B2/(20*a*a))).sqrt();old=((1+s*s).sqrt()+s)**10
  RU=min(upper(exp),upper(old));C=min(RU,1+lam/a+F(21,8)*eta*RU/(a*a));assert RU>=1 and C>=1
  L+=arb(rq(ml/RU));U+=arb(rq(mu*RU));M+=arb(rq(mu*C))
  bands.append({'mass_lower':str(ml),'mass_upper':str(mu),'R_upper':str(RU),'mean_factor':str(C),'exp_selected':upper(exp)<=upper(old),'taylor_selected':C<RU})
 l=F(ep(L));u=upper(U);m=min(u,upper(M));assert 0<l<=m<=u
 q=arb(5)/4;v=arb(rq(g['volume']));lp=arb(rq(l))**q;up=arb(rq(u))**q
 chord=v*(lp+(up-lp)*arb(rq((m-l)/(u-l))));baseline=v*arb(rq(rows[0]['H']['lower']))**q
 return {'H_lower':str(l),'H_upper':str(u),'mean_upper':str(m),'power_integral_lower':ep(v*lp),'chord_integral_upper':ep(chord,True),'constant_integral_upper':ep(v*up,True),'relative_center_baseline_upper':ep(chord/baseline,True),'positive_band_count':len(bands),'taylor_band_count':sum(b['taylor_selected'] for b in bands),'exp_band_count':sum(b['exp_selected'] for b in bands),'band_evidence_sha256':hashlib.sha256(json.dumps(bands,sort_keys=True,separators=(',',':')).encode()).hexdigest()}

def alarm_handler(sig,frame):raise TimeoutError('590 second hard alarm')

def worker(payload):
 row,reg,fingerprint,resume=payload;began=time.monotonic();signal.signal(signal.SIGALRM,alarm_handler);signal.alarm(590);ctx.prec=384
 out=HERE/(row['id']+'.json');cp=HERE/(row['id']+'_checkpoint.json')
 if out.exists():
  assert resume,'refuse replay existing result';r=json.loads(out.read_text());assert r['fingerprint']==fingerprint;return {'id':row['id'],'status':'reused_completed','pid':os.getpid()}
 state={'fingerprint':fingerprint,'id':row['id'],'stage':'stencil','next_cone':0,'weights':{},'endpoints':[]};runtime={'pid':os.getpid(),'began_unix':time.time(),'new_jet_count':0}
 try:
  if cp.exists():
   assert resume,'existing checkpoint requires resume';state=json.loads(cp.read_text());assert state['fingerprint']==fingerprint
  table=json.loads(gzip.decompress(Path(reg['sources']['table017']['path']).read_bytes()));verts=table['vertices'];cones=table['simplices'];assert len(verts)==280 and len(cones)==3504 and table['denominator']==144
  cinfo=row['candidates'][row['selected']];D=cinfo['D'];c=list(map(F,cinfo['point']));nodes=[sum(int(x*D)*z for x,z in zip(c,v)) for v in verts];assert len(set(nodes))==280
  weights={int(k):rq(v) for k,v in state['weights'].items()};stage=time.monotonic()
  for ci in range(state['next_cone'],3504):
   cone=cones[ci];ks=[nodes[k] for k in cone['vertices']];det=int(cone['det_numerator'])
   for j,k in enumerate(ks):
    w=-fmpq(det,2*prod(k-l for n,l in enumerate(ks) if n!=j));key=k*k;weights[key]=weights.get(key,fmpq(0))+w
   state['next_cone']=ci+1
   if (ci+1)%512==0 or ci==3503:
    state['weights']={str(k):str(v) for k,v in weights.items() if v};atomic(cp,state)
    if time.monotonic()-began>=540:return {'id':row['id'],'status':'checkpointed','pid':os.getpid()}
  assert len(weights)==140 and sum(weights.values(),fmpq(0))==0
  runtime['stencil_seconds']=time.monotonic()-stage;stage=time.monotonic();state['stage']='endpoints';factor=arb(D)**4/(144*(2*arb.pi())**5)
  for ts in reg['time_endpoints'][len(state['endpoints']):]:
   t=F(ts);val=factor*sum((arb(w)*arb(rq((144*D*t)**2+k2)).log() for k2,w in sorted(weights.items())),arb(0))
   assert val.is_finite() and val>0,('invalid_H',ts,str(val));state['endpoints'].append({'t':ts,'H':enc(val)});atomic(cp,state)
   if time.monotonic()-began>=540:return {'id':row['id'],'status':'checkpointed','pid':os.getpid()}
  runtime['endpoint_seconds']=time.monotonic()-stage;stage=time.monotonic();g=geometry(row);b=bounds(state['endpoints'],g);runtime['bounds_seconds']=time.monotonic()-stage
  # Exact vertex projection zero crossings are diagnostic; no guard required.
  crossings=0
  for z in verts:
   vlo=z[0]+sum(F(i if a>=0 else i+1,64)*a for i,a in zip(row['indices'],z[1:]));vhi=z[0]+sum(F(i+1 if a>=0 else i,64)*a for i,a in zip(row['indices'],z[1:]));crossings+=vlo<=0<=vhi
  runtime.update({'seconds':time.monotonic()-began,'peak_rss':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,'H_count':len(state['endpoints']),'cone_count':3504})
  r={'fingerprint':fingerprint,'id':row['id'],'class':row['class'],'indices':row['indices'],'weight':row['weight'],'anchor':cinfo['point'],'D':D,'geometry':g,'weights':{str(k):str(v) for k,v in sorted(weights.items())},'weight_sum':'0','endpoints':state['endpoints'],'bounds':b,'runtime':runtime,'zero_crossing_ray_count':crossings,'new_global_bound':False}
  atomic(out,r);cp.unlink();return {'id':row['id'],'status':'complete','pid':os.getpid(),'seconds':runtime['seconds']}
 except BaseException as exc:
  state['diagnostic']={'error':repr(exc),'traceback':traceback.format_exc(),'runtime':runtime,'elapsed_seconds':time.monotonic()-began};atomic(cp,state)
  return {'id':row['id'],'status':'failed','pid':os.getpid(),'error':repr(exc)}
 finally:signal.alarm(0)

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--resume',action='store_true');args=ap.parse_args();began=time.monotonic();reg=json.loads((HERE/'registration.json').read_text());assert sha(HERE/'mechanism.md')==reg['mechanism_sha256']
 for v in reg['sources'].values():assert sha(v['path'])==v['sha256'],v['path']
 fingerprint=hashlib.sha256((sha(__file__)+sha(HERE/'registration.json')).encode()).hexdigest()
 if (HERE/'run_receipt.json').exists():assert args.resume,'existing run requires resume'
 receipt={'fingerprint':fingerprint,'engine_sha256':sha(__file__),'registration_sha256':sha(HERE/'registration.json'),'pid':os.getpid(),'python':sys.executable,'flint':flint.__version__,'worker_count':12,'status':'running','jobs':[],'began_unix':time.time()};atomic(HERE/'run_receipt.json',receipt)
 signal.signal(signal.SIGALRM,alarm_handler);signal.alarm(590);pool=mp.get_context('spawn').Pool(12)
 try:
  for r in pool.imap_unordered(worker,[(row,reg,fingerprint,args.resume) for row in reg['cells']]):
   receipt['jobs'].append(r);atomic(HERE/'run_receipt.json',receipt);print(json.dumps(r),flush=True)
  pool.close();pool.join();receipt['workers_joined']=True
 except BaseException as e:
  pool.terminate();pool.join();receipt.update({'workers_joined':True,'exception':repr(e)})
 finally:signal.alarm(0)
 receipt.update({'status':'complete' if len(receipt['jobs'])==32 and all(r['status'] in ('complete','reused_completed') for r in receipt['jobs']) else 'incomplete','wall_seconds':time.monotonic()-began,'worker_pids':sorted({r['pid'] for r in receipt['jobs']}),'peak_parent_rss':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,'children_user_seconds':resource.getrusage(resource.RUSAGE_CHILDREN).ru_utime,'children_system_seconds':resource.getrusage(resource.RUSAGE_CHILDREN).ru_stime});atomic(HERE/'run_receipt.json',receipt)
 print(json.dumps(receipt),flush=True)
if __name__=='__main__':main()
