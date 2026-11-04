#!/usr/bin/env python3
"""High-precision cube diagnostics for a=(rho/|x1|)^gamma, 0<gamma<1.

All formulae are explicit positive/integrable shell calculations.  Output is
diagnostic evidence for the proposed unbounded-angular extension, not a proof
of its all-input maximal inequality.
"""
import os
for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS'): os.environ[k]='1'
import hashlib,json,os as _os,platform,sys,time
from concurrent.futures import ProcessPoolExecutor,as_completed
from decimal import Decimal as D,getcontext
from multiprocessing import Manager
from pathlib import Path
ROOT=Path(__file__).resolve().parent;CACHE=ROOT/'checkpoints';VERSION='round28-unbounded-cube-v3-20261105';PREC=120
NS=(2,3);GS=('.25','.5','.75','.9');TS=('.1','.01','.0001','1e-6')
def source_hash():return hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
def pw(x,e):
 if x==0:return D(0)
 assert x>0
 return (e*x.ln()).exp()
def angular_direct(n,g):
 """Average over cube: integrate first x1>=0 then exploit symmetry."""
 d=D(n-1);nn=D(n)
 return D(1)/nn+d/(d+g)*(D(1)/(1-g)-D(1)/nn)
def primitive_small(a,n,g,t,r):
 """∫_0^a x^-g ∫_{[-r,r]^(n-1)} max(x,||y||)^alpha dy dx, a<=r."""
 d=D(n-1);q=t+g-1;c=D(2)**(n-1)
 if q==0:
  # B(x)=2^d[1+d log(r/x)] after multiplying the y integral by x^-g.
  z=1-g
  return c*pw(a,z)*((1+d*(r/a).ln())/z+d/(z*z))
 return c*((1-d/q)*pw(a,t)/t+d*pw(r,q)/q*pw(a,1-g)/(1-g))
def primitive_x(L,n,g,t,r):
 """Integral from x=0 to L of the positive x1-slice contribution."""
 if L<=r:return primitive_small(L,n,g,t,r)
 q0=primitive_small(r,n,g,t,r);e=t-n+1
 return q0+(D(2)*r)**(n-1)*(pw(L,e)-pw(r,e))/e
def atom_mass(n,g,t):
 r=1+t
 # x1 ranges [-t,2+t]: reflect negative portion to [0,t].
 mu=primitive_x(2+t,n,g,t,r)+primitive_x(t,n,g,t,r)
 cost=(D(2)**n)*r**n
 return mu/cost,mu,cost
def local_a1_bound(n,g,t):
 # For Q=x+rK with r<=rho(x)/2, the radial factor varies by at most
 # 3^|n-t-g|.  A general (not necessarily symmetric) x1-interval crossing
 # zero needs the safe A1 bound 2/(1-g), not the symmetric 1/(1-g).
 return D(2)*pw(D(3),abs(D(n)-t-g))/(1-g)
def sig(n,gs,ts,run_mode):
 p={'version':VERSION,'source_sha256':source_hash(),'run_mode':run_mode,'n':n,'gamma':gs,'tau':ts,'precision':PREC,'python':sys.version,'platform':platform.platform(),'model':'cube angular exact + single atom positive x1 integration'}
 return hashlib.sha256(json.dumps(p,sort_keys=True).encode()).hexdigest(),p
def job(item,barrier=None,run_mode='full'):
 n,gs,ts=item;path=CACHE/f'{run_mode}_n{n}_g{gs.replace(".","p")}_t{ts.replace(".","p")}.json';s,p=sig(*item,run_mode)
 if barrier is not None:barrier.wait(timeout=60)
 if path.exists():
  old=json.loads(path.read_text())
  if old.get('signature')==s:old['cache_state']='reused';old['replay_pid']=_os.getpid();return old
  stale=CACHE/'stale';stale.mkdir(exist_ok=True);path.replace(stale/(path.stem+'_stale_'+str(time.time_ns())+'.json'))
 getcontext().prec=PREC;g=D(gs);t=D(ts);ea=angular_direct(n,g);expect=(D(n)-g)/(1-g);R,mu,cost=atom_mass(n,g,t);a1=local_a1_bound(n,g,t)
 out={'n':n,'gamma':gs,'tau':ts,'a_star':'1','Ea':str(ea),'Ea_prediction':str((D(n)-g)/(D(n)*(1-g))),'angular_error':str(ea-(D(n)-g)/(D(n)*(1-g))),'target_tauW':str(expect),'atom_mass':str(mu),'atom_cost':str(cost),'atom_R':str(R),'tau_atom_R':str(t*R),'target_error':str(t*R-expect),'local_A1_bound':str(a1),'worker_pid':_os.getpid(),'cache_state':'computed','signature':s,'signature_payload':p,'valid':ea>1 and R>0 and a1>1}
 tmp=path.with_suffix('.tmp');tmp.write_text(json.dumps(out,indent=2));tmp.replace(path);return out
def run(smoke=False):
 items=[(n,g,t) for n in NS for g in GS for t in TS]
 if smoke:items=[(2,'.9','.1')]
 CACHE.mkdir(exist_ok=True);rows=[];bad=[];workers=1 if smoke else 15
 print(json.dumps({'event':'start','tasks':len(items),'workers':workers,'precision':PREC,'smoke':smoke}),flush=True)
 if smoke:
  for x in items:
   try:rows.append(job(x,run_mode='smoke'));print(json.dumps({'event':'done','row':rows[-1]}),flush=True)
   except Exception as e:bad.append({'task':x,'error':repr(e)})
 else:
  with Manager() as manager:
   barrier=manager.Barrier(15)
   with ProcessPoolExecutor(max_workers=15) as ex:
    fs={ex.submit(job,x,barrier):x for x in items[:15]}
    for f in as_completed(fs):
     try:z=f.result();rows.append(z);print(json.dumps({'event':'done','n':z['n'],'gamma':z['gamma'],'tau':z['tau'],'pid':z['worker_pid'],'tauR':z['tau_atom_R'],'cache':z['cache_state']}),flush=True)
     except Exception as e:bad.append({'task':fs[f],'error':repr(e)});print(json.dumps({'event':'failure','task':fs[f],'error':repr(e)}),flush=True)
    fs={ex.submit(job,x):x for x in items[15:]}
    for f in as_completed(fs):
     try:z=f.result();rows.append(z);print(json.dumps({'event':'done','n':z['n'],'gamma':z['gamma'],'tau':z['tau'],'pid':z['worker_pid'],'tauR':z['tau_atom_R'],'cache':z['cache_state']}),flush=True)
     except Exception as e:bad.append({'task':fs[f],'error':repr(e)});print(json.dumps({'event':'failure','task':fs[f],'error':repr(e)}),flush=True)
 rows.sort(key=lambda x:(x['n'],D(x['gamma']),D(x['tau'])));out={'version':VERSION,'precision':PREC,'rows':rows,'failures':bad,'unique_inputs':len(items),'configured_workers':workers,'actual_worker_pids':sorted(set(x['worker_pid'] for x in rows))}
 suffix='smoke' if smoke else 'summary';(ROOT/f'unbounded_cube_{suffix}.json').write_text(json.dumps(out,indent=2));print(json.dumps({'event':'finish','rows':len(rows),'failures':len(bad)}),flush=True)
 if bad:raise SystemExit(1)
if __name__=='__main__':
 import argparse
 ap=argparse.ArgumentParser();ap.add_argument('--smoke',action='store_true');run(ap.parse_args().smoke)
