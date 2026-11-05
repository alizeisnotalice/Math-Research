#!/usr/bin/env python3
"""120+ digit numerical sign-bracket diagnostics for flat concentration.

Finite brackets are numerical diagnostics, not interval certificates or a
proof of the full translated-cube concentration theorem.
"""
import os
for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS'):os.environ[key]='1'
import hashlib,json,platform,sys
from concurrent.futures import ProcessPoolExecutor,as_completed
from decimal import Decimal as D,getcontext,ROUND_CEILING
from multiprocessing import Manager
from pathlib import Path
ROOT=Path(__file__).resolve().parent;CACHE=ROOT/'runtime_cache';PREC=120;VERSION='flat-concentration-round32-v2-adaptive-20261105'
BETAS=('.1','.25','.5','.75','.9');KS=(4,8,16,32,64,128,256)
def sha():return hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
def pw(x,p):return (p*x.ln()).exp()
def quantities(a,u,beta):
 M=pw(a,beta)+pw(u-a,beta);Dn=pw(a,beta)+pw(1-a,beta);Mp=beta*(pw(a,beta-1)-pw(u-a,beta-1));Dp=beta*(pw(a,beta-1)-pw(1-a,beta-1));return M,Dn,Mp*Dn-M*Dp
def solve_b(u,beta):
 lo=D(10)**(-D(PREC-20));hi=D(1)/2-D(10)**(-D(PREC-20));_,_,initial_lo=quantities(u*lo,u,beta);_,_,initial_hi=quantities(u*hi,u,beta);assert initial_lo>0 and initial_hi<0
 stagnated=False;zero_mid=False;iterations=0
 for iterations in range(700):
  mid=(lo+hi)/2
  if mid==lo or mid==hi:stagnated=True;break
  _,_,value=quantities(u*mid,u,beta)
  if value>0:lo=mid
  elif value<0:hi=mid
  else:zero_mid=True;break # retain previous strict numerical bracket
 _,_,final_lo=quantities(u*lo,u,beta);_,_,final_hi=quantities(u*hi,u,beta);assert final_lo>0 and final_hi<0
 return lo,hi,initial_lo,initial_hi,final_lo,final_hi,stagnated,zero_mid,iterations
def one(beta_s,kind,k):
 loss=int((D(k)*D(2).ln()/D(10).ln()).to_integral_value(rounding=ROUND_CEILING));workprec=PREC if kind=='small' else max(PREC,PREC+loss+20);getcontext().prec=workprec
 beta=D(beta_s);eps=(-D(k)*D(2).ln()).exp();u=eps if kind=='small' else D(1)-eps;lo,hi,KiLo,KiHi,KfLo,KfHi,stagnated,zero_mid,iterations=solve_b(u,beta);b=(lo+hi)/2;a=u*b;M,Dn,Kmid=quantities(a,u,beta);omega=M/Dn
 expansion_status='not_applicable';remainder=None;scaled=None
 if kind=='small':
  leading=D(2)**(D(1)-beta)*pw(u,beta);second=D(2)**(D(1)-2*beta)*pw(u,2*beta);remainder=omega-leading+second;scale=max(abs(omega),abs(leading),abs(second));expansion_status='unresolved_at_120_digits' if abs(remainder)<=scale*D(10)**(-D(PREC-20)) else 'resolved';scaled=None if expansion_status!='resolved' else remainder/pw(u,2*beta)
 half={}
 if beta==D('.5'):
  r=(D(2)*u).sqrt();aclosed=u/(D(2)*(D(1)+r+u));oclosed=(u+r)/(D(1)+r);closedb=aclosed/u;btol=D(4)*(hi-lo)+D(10)**(-D(90));otol=min(oclosed,D(1)-oclosed)*D(10)**(-D(80));assert abs(b-closedb)<=btol and abs(omega-oclosed)<=otol
  half={'half_b_closed':str(closedb),'half_omega_closed':str(oclosed),'half_b_discrepancy':str(b-closedb),'half_omega_discrepancy':str(omega-oclosed),'half_b_tolerance':str(btol),'half_omega_tolerance':str(otol)}
 scaleK=(D(1)-u) if kind=='near_one' else D(1)
 return {'beta':beta_s,'u_kind':kind,'k':k,'working_precision':workprec,'u':str(u),'b_lo':str(lo),'b_hi':str(hi),'b_width':str(hi-lo),'K_initial_lo':str(KiLo),'K_initial_hi':str(KiHi),'K_final_lo':str(KfLo),'K_final_hi':str(KfHi),'K_final_lo_scaled':str(KfLo/scaleK),'K_final_hi_scaled':str(KfHi/scaleK),'K_mid':str(Kmid),'iterations':iterations,'rounding_stagnation':stagnated,'zero_mid_evaluation':zero_mid,'omega':str(omega),'expansion_status':expansion_status,'expansion_remainder':None if remainder is None else str(remainder),'expansion_scaled':None if scaled is None else str(scaled)}|half
def signature(item):return hashlib.sha256(json.dumps({'version':VERSION,'source_sha256':sha(),'item':item,'precision':PREC,'python':sys.version,'platform':platform.platform()},sort_keys=True).encode()).hexdigest()
def job(item,barrier=None):
 beta,kind,k=item;path=CACHE/f'b{beta.replace(".","p")}_{kind}_k{k}.json';sig=signature(item)
 if barrier is not None:barrier.wait(timeout=60)
 if path.exists():
  old=json.loads(path.read_text())
  if old.get('signature')==sig:old['cache']='reused';old['replay_pid']=os.getpid();return old
 out=one(beta,kind,k)|{'signature':sig,'pid':os.getpid(),'cache':'computed'};tmp=path.with_suffix('.tmp');tmp.write_text(json.dumps(out));tmp.replace(path);return out
def main():
 CACHE.mkdir(exist_ok=True);items=[(b,kind,k) for b in BETAS for kind in ('small','near_one') for k in KS];rows=[];bad=[]
 with Manager() as manager:
  barrier=manager.Barrier(15)
  with ProcessPoolExecutor(max_workers=15) as ex:
   first={ex.submit(job,x,barrier):x for x in items[:15]}
   for f in as_completed(first):
    try:rows.append(f.result())
    except Exception as exc:bad.append({'item':first[f],'error':repr(exc)})
   rest={ex.submit(job,x):x for x in items[15:]}
   for f in as_completed(rest):
    try:rows.append(f.result())
    except Exception as exc:bad.append({'item':rest[f],'error':repr(exc)})
 rows.sort(key=lambda r:(D(r['beta']),r['u_kind'],r['k']));out={'version':VERSION,'tasks':len(items),'precision':PREC,'rows':rows,'failures':bad,'worker_pids':sorted(set(r['pid'] for r in rows))};(ROOT/'round32_summary.json').write_text(json.dumps(out,indent=2));print(json.dumps({'tasks':len(rows),'failures':len(bad),'workers':out['worker_pids']}));
 if bad:raise SystemExit(1)
if __name__=='__main__':main()
