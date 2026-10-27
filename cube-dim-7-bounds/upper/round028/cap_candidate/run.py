#!/usr/bin/env python3
from fractions import Fraction as Q
from pathlib import Path
from math import comb,log,sqrt,isfinite
from flint import arb,fmpq,ctx
from scipy.optimize import minimize_scalar
import hashlib,json,time,datetime,resource,sys,flint,scipy
ctx.prec=512;D=Path(__file__).resolve().parent
now=lambda:datetime.datetime.now(datetime.timezone.utc).isoformat()
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def fr(p):return Q(int(p[0]))*Q(2)**int(p[1])
def a(x):return arb(fmpq(x.numerator,x.denominator)) if isinstance(x,Q) else arb(x)
def ep(x):
 l=x.lower().man_exp();u=x.upper().man_exp();return {'lower':[str(l[0]),int(l[1])],'upper':[str(u[0]),int(u[1])],'display':x.str(45)}
def grid(iv):
 u=fr(iv['upper']);n=u.numerator*10**12//u.denominator+1;assert u<Q(n,10**12);return str(Q(n,10**12))
def coeff(t):
 d=1+1/(2*t);h=(1+1/t)/2;ts=[Q(comb(6,k))*d**(6-k)*Q(-1,2)**k*(1-t**(-k))/k for k in range(1,7)];return 7-Q(7,2)*d**6,1-h**7-Q(7,2)*sum(ts)
def main():
 assert not (D/'trace.jsonl').exists();start=now();w=time.perf_counter();c=time.process_time();status='failed';err=None;new=0;hits=0;cb=0;cert=0
 try:
  s=D/'sources';mf=json.loads((s/'manifest.json').read_text());assert all(sha(s/p)==mf['files'][p] for p in ('result.json','terminal.jsonl'));assert sha(s/'author.md')=='6ecbf72fa9213677dc99e2e2afaa510df240b3d1b08047366a8d0103d010c80f';assert sha(s/'audit.md')=='00555fc5ca9dff8cdb3e467ef0bf4230b5eed8b603a744256f9a5fb66224db6a';old=json.loads((s/'result.json').read_text());rows=[json.loads(x) for x in (s/'terminal.jsonl').read_text().splitlines()];assert old['status']=='completed' and len(rows)==5 and all(x['status']=='completed' for x in rows) and old['winner']==rows[3] and rows[3]['t']=='6'
  cache={float(Q(r['t'])):(float((fr(r['F']['lower'])+fr(r['F']['upper']))/2),r) for r in rows};assert len(cache)==5
  (D/'cache.json').write_text(json.dumps([{'t':r['t'],'float_t_hex':float(Q(r['t'])).hex(),'midpoint_float_hex':cache[float(Q(r['t']))][0].hex(),'interval':r['F']} for r in rows],indent=2)+'\n')
  def objective(x):
   nonlocal new,hits,cb
   x=float(x);cb+=1;hit=x in cache
   if hit:y=cache[x][0];hits+=1
   else:
    if new>=64:raise RuntimeError('64-new-evaluation budget exhausted')
    new+=1;t=Q.from_float(x);k,z=coeff(t);B=float(z)+float(k)*log(x)
    if not isfinite(B) or B<=0:raise ArithmeticError('nonpositive/nonfinite exploratory B')
    y=(sqrt(float((1+1/t)**7))+sqrt(B))**2
   if not isfinite(y):raise ArithmeticError('nonfinite exploratory F')
   with (D/'trace.jsonl').open('a') as f:f.write(json.dumps({'call':cb,'x_hex':x.hex(),'y_hex':float(y).hex(),'cache_hit':hit,'new_eval_count':new})+'\n');f.flush()
   return y
  opt=minimize_scalar(objective,bounds=(4,10),method='bounded',options={'maxiter':64,'xatol':1e-10});assert opt.success and new<=64 and 4<=opt.x<=10
  raw=float(opt.x);q=round(Q.from_float(raw)*2**20);t=Q(q,2**20);assert 4<=t<=10
  if float(t) in cache and t==Q.from_float(float(t)):
   row=cache[float(t)][1];iv=row['F'];cert_source='R027 exact old row';coef,const=coeff(t)
  else:
   cert+=1;coef,const=coeff(t);L=a(t).log();B=a(const)+a(coef)*L;gap=1+7*L-B;assert B>0 and gap>0;v=(a((1+1/t)**7).sqrt()+B.sqrt())**2;assert v>0 and v.is_finite();iv=ep(v);cert_source='new 512-bit Arb row'
  result=dict(status='completed',opt_success=bool(opt.success),opt_message=str(opt.message),opt_nfev=int(opt.nfev),callback_calls=cb,new_double_evaluations=new,cache_hits=hits,raw_t_hex=raw.hex(),raw_value_hex=float(opt.fun).hex(),integer_q=q,denominator=2**20,exact_t=str(t),cert_source=cert_source,new_arb_rows=cert,log_coefficient=str(coef),rational_constant=str(const),candidate_interval=iv,candidate_grid_upper=grid(iv),strict_below_R027=fr(iv['upper'])<fr(old['winner']['F']['lower']),R027_winner_interval=old['winner']['F'],scope='local bounded numerical candidate only; no claim of a global continuous optimum or C7 lower bound')
  (D/'result.json').write_text(json.dumps(result,indent=2)+'\n');status='completed'
 except BaseException as e:err=repr(e)
 rc=dict(status=status,error=err,start_utc=start,finish_utc=now(),wall_seconds=time.perf_counter()-w,cpu_seconds=time.process_time()-c,peak_rss_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,callback_calls=cb,new_double_evaluations=new,cache_hits=hits,new_arb_rows=cert,python=sys.version,flint=flint.__version__,scipy=scipy.__version__,precision_bits=512)
 (D/'receipt.json').write_text(json.dumps(rc,indent=2)+'\n');print(json.dumps(rc,indent=2));sys.exit(0 if status=='completed' else 1)
if __name__=='__main__':main()
