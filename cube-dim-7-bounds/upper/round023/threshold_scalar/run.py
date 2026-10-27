#!/usr/bin/env python3
from flint import arb,fmpq,ctx
from fractions import Fraction
from pathlib import Path
import json,hashlib,time,datetime,resource,sys,flint
ctx.prec=512;D=Path(__file__).resolve().parent;now=lambda:datetime.datetime.now(datetime.timezone.utc).isoformat()
def ep(x):
 l=x.lower().man_exp();u=x.upper().man_exp();return {'lower':[str(l[0]),int(l[1])],'upper':[str(u[0]),int(u[1])],'display':x.str(45)}
def frac(v):return Fraction(int(v[0]))*Fraction(2)**int(v[1])
def main():
 assert not (D/'terminal.jsonl').exists();t=time.perf_counter();cpu=time.process_time();start=now();old=json.loads((D/'r022-result.json').read_text());manifest=json.loads((D/'r022-manifest.json').read_text());assert hashlib.sha256((D/'r022-result.json').read_bytes()).hexdigest()==manifest['files']['result.json'] and old['status']=='completed';oldlo=Fraction(old['winner_lower_exact']);oldhi=Fraction(old['winner_upper_exact']);assert oldlo<oldhi
 out=[]
 try:
  for a in (32,36,40,64):
   A=arb(fmpq(a+2,a))**7;B=arb(fmpq(a+1,a))**7*(1+7*arb(a).log())
   for tau in (6,7):
    U=(1+A/tau)*(tau+B);alt=A+B+tau+A*B/tau;assert U>0 and U.is_finite() and U.overlaps(alt);e=ep(U);upper=frac(e['upper']);N=10**12;grid=Fraction(upper.numerator*N//upper.denominator+1,N);assert upper<grid;row=dict(a=a,tau=tau,status='completed',interval=e,grid_upper=str(grid),algebraic_overlap=True,strictly_below_R022=upper<oldlo)
    out.append(row)
    with (D/'terminal.jsonl').open('a') as f:f.write(json.dumps(row)+'\n');f.flush()
  assert frac(next(x for x in out if (x['a'],x['tau'])==(64,6))['interval']['upper'])<Fraction(1040858,21375)
  winner=min(out,key=lambda r:Fraction(r['grid_upper']));wu=frac(winner['interval']['upper']);assert all(wu<frac(x['interval']['lower']) for x in out if x is not winner) and wu<oldlo
  N=10**12;k=wu.numerator*N//wu.denominator+1;assert wu<Fraction(k,N)
  r=dict(status='completed',dimension=7,precision_bits=512,evaluations=8,winner=dict(a=winner['a'],tau=winner['tau'],interval=winner['interval']),winner_strict_best_of_registered_8=True,strictly_below_R022_full_interval=True,old_R022_a=old['winner_a'],old_R022_lower=str(oldlo),old_R022_upper=str(oldhi),grid_upper=str(Fraction(k,N)),grid_denominator=N,continuous_optimum_claim=False)
  (D/'result.json').write_text(json.dumps(r,indent=2)+'\n');status='completed';error=None
 except BaseException as e:status='failed';error=repr(e)
 rc=dict(status=status,error=error,start_utc=start,finish_utc=now(),wall_seconds=time.perf_counter()-t,cpu_seconds=time.process_time()-cpu,peak_rss_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,terminal_pairs=len(out),python=sys.version,flint=flint.__version__,precision_bits=512);(D/'receipt.json').write_text(json.dumps(rc,indent=2)+'\n');print(json.dumps(rc,indent=2));sys.exit(0 if status=='completed' else 1)
if __name__=='__main__':main()
