#!/usr/bin/env python3
"""120-digit scalar diagnostics for the supercritical angular hole.

The finite table checks the stated scalar formula only; it is not a proof of
the all-input maximal-operator theorem. Generated cache and full tables are
deliberately not versioned.
"""
import os
for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS'): os.environ[key]='1'
import hashlib,json,platform,sys
from concurrent.futures import ProcessPoolExecutor,as_completed
from decimal import Decimal as D,getcontext
from multiprocessing import Manager
from pathlib import Path
ROOT=Path(__file__).resolve().parent; CACHE=ROOT/'runtime_cache'; PREC=120
VERSION='supercritical-angular-round31-v2-log-domain-20261105'
NS=(2,3); GS=('.25','.5','.75','.9'); CS=('.01','.05','.1','.25','.5','.75','.95'); BASE_K=(16,64,256,1024); EXTRA_K=(4096,16384)
def sha(): return hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
def pw(x,y): return (y*x.ln()).exp()
def parameters(n,gs,cs):
 gamma=D(gs); beta=1-gamma; d=D(n)*gamma*D(cs); tau=D(n)*beta+d; e=D(n)-tau; P=d/gamma; tail_power=d/e; p=D(n)*gamma/e; coeff=D(n-1)/P
 mass1=tau*gamma/e; mass2=tau*D(n-1)/(P*(D(n)-gamma)); A=tau/(D(n)-gamma)*pw(coeff,gamma)
 assert 0<d and D(n)*beta<tau<D(n) and all(x>0 for x in (e,P,tail_power,p,coeff,mass1,mass2,A))
 return gamma,beta,d,tau,e,P,tail_power,p,coeff,mass1,mass2,A
NEG_INF=D('-Infinity')
def logadd(x,y):
 if x==NEG_INF:return y
 if y==NEG_INF:return x
 high,low=(x,y) if x>=y else (y,x); return high+(D(1)+(low-high).exp()).ln()
def log1mexp(x): return NEG_INF if x==0 else (D(1)-(-x).exp()).ln()
def log_volume(L,p,coeff,tail_power): return logadd(-p*L,coeff.ln()-L+log1mexp(tail_power*L))
def log_mass(L,beta,mass1,mass2,tail_power): return logadd(-mass1*L,mass2.ln()-beta*L+log1mexp(tail_power*L))
def derivative(a,p,coeff,tail_power): return p*pw(a,p-1)+coeff*(1-(1+tail_power)*pw(a,tail_power))
def solve_L(H,p,coeff,tail_power):
 lo=D(0); hi=H+coeff.ln()+10; assert log_volume(lo,p,coeff,tail_power)>-H and log_volume(hi,p,coeff,tail_power)<-H
 for _ in range(800):
  mid=(lo+hi)/2
  if log_volume(mid,p,coeff,tail_power)>-H:lo=mid
  else:hi=mid
 L=(lo+hi)/2; residual=log_volume(L,p,coeff,tail_power)+H; assert abs(residual)<D('1e-100')
 return L,residual
def choose_delta(d,tau): return (d/tau,'interior_d_over_tau') if d/tau<D(1)/3 else (D(1)/4,'conservative_quarter')
def data(n,gs,cs,k):
 getcontext().prec=PREC; gamma,beta,d,tau,e,P,tail_power,p,coeff,mass1,mass2,A=parameters(n,gs,cs)
 H=D(k)*D(2).ln(); L,res=solve_L(H,p,coeff,tail_power); normalized=(log_mass(L,beta,mass1,mass2,tail_power)+beta*H).exp()
 delta,policy=choose_delta(d,tau); B=A*pw(delta,d)*pw(1-delta,D(n)*beta)
 sampled_min=min(derivative(pw(D(10),-D(j)/2),p,coeff,tail_power) for j in range(81)); assert sampled_min>0 and 0<delta<D(1)/3
 return {'n':n,'gamma':gs,'c':cs,'k':k,'d':str(d),'tau':str(tau),'P':str(P),'log_inverse_residual':str(res),'A':str(A),'normalized_F':str(normalized),'normalization_error':str(normalized-A),'delta':str(delta),'delta_policy':policy,'B_delta':str(B),'sampled_derivative_min':str(sampled_min),'valid':True}
def signature(item): return hashlib.sha256(json.dumps({'version':VERSION,'source_sha256':sha(),'item':item,'precision':PREC,'python':sys.version,'platform':platform.platform()},sort_keys=True).encode()).hexdigest()
def job(item,barrier=None):
 n,g,c,k=item; target=CACHE/f'n{n}_g{g.replace(".","p")}_c{c.replace(".","p")}_k{k}.json'; sig=signature(item)
 if barrier is not None:barrier.wait(timeout=60)
 if target.exists():
  old=json.loads(target.read_text())
  if old.get('signature')==sig:old['cache']='reused';old['replay_pid']=os.getpid();return old
 out=data(n,g,c,k)|{'signature':sig,'pid':os.getpid(),'cache':'computed'}; tmp=target.with_suffix('.tmp');tmp.write_text(json.dumps(out));tmp.replace(target);return out
def main():
 CACHE.mkdir(exist_ok=True); items=[(n,g,c,k) for n in NS for g in GS for c in CS for k in BASE_K+(EXTRA_K if c=='.01' else ())];rows=[];failures=[]
 with Manager() as manager:
  barrier=manager.Barrier(15)
  with ProcessPoolExecutor(max_workers=15) as executor:
   first={executor.submit(job,x,barrier):x for x in items[:15]}
   for future in as_completed(first):
    try:rows.append(future.result())
    except Exception as exc:failures.append({'item':first[future],'error':repr(exc)})
   rest={executor.submit(job,x):x for x in items[15:]}
   for future in as_completed(rest):
    try:rows.append(future.result())
    except Exception as exc:failures.append({'item':rest[future],'error':repr(exc)})
 rows.sort(key=lambda r:(r['n'],D(r['gamma']),D(r['c']),r['k']));result={'version':VERSION,'tasks':len(items),'precision':PREC,'rows':rows,'failures':failures,'worker_pids':sorted(set(r['pid'] for r in rows))}
 (ROOT/'round31_summary.json').write_text(json.dumps(result,indent=2));print(json.dumps({'tasks':len(rows),'failures':len(failures),'workers':result['worker_pids']}))
 if failures:raise SystemExit(1)
if __name__=='__main__':main()
