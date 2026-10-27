#!/usr/bin/env python3
from flint import arb,fmpq,ctx
from fractions import Fraction
from pathlib import Path
import hashlib,json,time,datetime,resource,sys,flint
ctx.prec=512;D=Path(__file__).resolve().parent;now=lambda:datetime.datetime.now(datetime.timezone.utc).isoformat()
def A(x):return arb(fmpq(x.numerator,x.denominator))
def ep(x):
 l=x.lower().man_exp();u=x.upper().man_exp();return {'lower':[str(l[0]),int(l[1])],'upper':[str(u[0]),int(u[1])],'display':x.str(45)}
def fr(x):return Fraction(int(x[0]))*Fraction(2)**int(x[1])
def grid(e):
 u=fr(e['upper']);N=10**12;return str(Fraction(u.numerator*N//u.denominator+1,N))
def main():
 assert not (D/'terminal.jsonl').exists();start=now();t=time.perf_counter();cpu=time.process_time();rows=[];error=None;status='failed'
 def record(r):
  rows.append(r)
  with (D/'terminal.jsonl').open('a') as f:f.write(json.dumps(r)+'\n');f.flush()
 try:
  old=json.loads((D/'r024-result.json').read_text());mf=json.loads((D/'r024-manifest.json').read_text());assert hashlib.sha256((D/'r024-result.json').read_bytes()).hexdigest()==mf['files']['result.json'] and old['status']=='completed';oldlo=fr(old['F_candidate']['lower']);logs={}
  for i,x in enumerate([Fraction(5),Fraction(11,2),Fraction(6)]):
   logs[x]=A(x).log();base=A(1+1/x);b=1+7*logs[x];assert base>0 and b>0;rootA=base**3*base.sqrt();rootB=b.sqrt();v=(rootA+rootB)**2;e=ep(v);assert v>0 and v.is_finite();row=dict(id=i,kind='optimized_tau',t=str(x),status='completed',interval=e,grid_upper=grid(e),strict_below_R024=fr(e['upper'])<oldlo);record(row)
  x=Fraction(5);tau=Fraction(20,3);a=A(1+1/x)**7;b=1+7*logs[x];v=(1+a/A(tau))*(A(tau)+b);e=ep(v);assert v>0 and v.is_finite();record(dict(id=3,kind='rational_tau_calibration',t='5',tau='20/3',status='completed',interval=e));assert fr(rows[0]['interval']['upper'])<fr(e['lower']) and fr(e['upper'])<Fraction(2233,75)
  win=min(rows[:3],key=lambda r:Fraction(r['grid_upper']));u=fr(win['interval']['upper']);assert all(u<fr(r['interval']['lower']) for r in rows[:3] if r is not win) and u<oldlo
  out=dict(status='completed',dimension=7,precision_bits=512,scientific_evaluations=4,optimized_rows=3,calibration_rows=1,log_calls=3,sqrt_calls=6,winner=dict(t=win['t'],interval=win['interval'],grid_upper=win['grid_upper']),winner_strictly_best_registered=True,strictly_below_R024=True,R024_candidate_interval=old['F_candidate'],continuous_t_optimum_claim=False)
  (D/'result.json').write_text(json.dumps(out,indent=2)+'\n');status='completed'
 except BaseException as e:error=repr(e)
 rc=dict(status=status,error=error,start_utc=start,finish_utc=now(),wall_seconds=time.perf_counter()-t,cpu_seconds=time.process_time()-cpu,peak_rss_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,terminal_rows=len(rows),python=sys.version,flint=flint.__version__,precision_bits=512);(D/'receipt.json').write_text(json.dumps(rc,indent=2)+'\n');print(json.dumps(rc,indent=2));sys.exit(0 if status=='completed' else 1)
if __name__=='__main__':main()
