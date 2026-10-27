#!/usr/bin/env python3
from flint import arb,fmpq,ctx
from fractions import Fraction
from pathlib import Path
import hashlib,json,time,datetime,resource,sys,flint
ctx.prec=512;D=Path(__file__).resolve().parent;now=lambda:datetime.datetime.now(datetime.timezone.utc).isoformat();logs={};calls={'log':0,'sqrt':0}
def A(x):return arb(fmpq(x.numerator,x.denominator)) if isinstance(x,Fraction) else arb(x)
def lg(x):
 if x not in logs:calls['log']+=1;logs[x]=A(x).log()
 return logs[x]
def sq(x):
 assert x>0;calls['sqrt']+=1;return x.sqrt()
def ep(x):
 l=x.lower().man_exp();u=x.upper().man_exp();return {'lower':[str(l[0]),int(l[1])],'upper':[str(u[0]),int(u[1])],'display':x.str(45)}
def fr(a):return Fraction(int(a[0]))*Fraction(2)**int(a[1])
def grid(e):
 u=fr(e['upper']);N=10**12;n=u.numerator*N//u.denominator+1;assert u<Fraction(n,N);return str(Fraction(n,N))
def evalF(x,fixed=False):
 z=lg(x)
 if fixed:
  assert x==36;r=A(Fraction(19,18));s=A(Fraction(37,36));v=(r**3*sq(r)+s**3*sq(s)*sq(1+7*z))**2
 else:
  aa=A(1+2/x)**7;bb=A(1+1/x)**7*(1+7*z);v=(sq(aa)+sq(bb))**2;assert v.overlaps(aa+bb+2*sq(aa*bb))
 assert v>0 and v.is_finite();return v
def evalR(x):
 z=lg(x);r=A((x+1)/(x+2));v=(r*r*sq(r))*(A(x)-7*z)/(2*sq(1+7*z))-1;assert v.is_finite();return v
def main():
 assert not (D/'terminal.jsonl').exists();start=now();t=time.perf_counter();cpu=time.process_time();out=[];status='failed';error=None
 def record(r):
  out.append(r)
  with (D/'terminal.jsonl').open('a') as f:f.write(json.dumps(r)+'\n');f.flush()
 try:
  mf=json.loads((D/'r023-manifest.json').read_text());assert hashlib.sha256((D/'r023-result.json').read_bytes()).hexdigest()==mf['files']['result.json'];old=json.loads((D/'r023-result.json').read_text());assert old['status']=='completed';oldlo=fr(old['winner']['interval']['lower']);oldhi=fr(old['winner']['interval']['upper']);assert oldlo<oldhi
  F36=ep(evalF(Fraction(36),True));record(dict(id=0,kind='F_fixed',a='36',status='completed',interval=F36,grid_upper=grid(F36),strictly_below_R023=fr(F36['upper'])<oldlo))
  lo=Fraction(35);hi=Fraction(37);rlo=evalR(lo);rhi=evalR(hi);record(dict(id=1,kind='R_endpoint',a='35',status='completed',interval=ep(rlo),sign=-1 if rlo<0 else 0));record(dict(id=2,kind='R_endpoint',a='37',status='completed',interval=ep(rhi),sign=1 if rhi>0 else 0));assert rlo<0 and rhi>0
  for i in range(80):
   mid=(lo+hi)/2;v=evalR(mid);sign=-1 if v<0 else 1 if v>0 else 0
   if sign==0:record(dict(id=i+3,kind='R_midpoint',step=i+1,a=str(mid),status='unresolved_sign',interval=ep(v),sign=0));raise ArithmeticError('R interval contains zero')
   if sign<0:lo=mid
   else:hi=mid
   record(dict(id=i+3,kind='R_midpoint',step=i+1,a=str(mid),status='completed',interval=ep(v),sign=sign,bracket_lower=str(lo),bracket_upper=str(hi)))
  assert hi-lo==Fraction(1,2**79);candidate=(lo+hi)/2;Fc=ep(evalF(candidate));record(dict(id=83,kind='F_candidate',a=str(candidate),status='completed',interval=Fc,grid_upper=grid(Fc),strictly_below_R023=fr(Fc['upper'])<oldlo))
  assert len(out)==84;result=dict(status='completed',precision_bits=512,function_evaluations=84,F36=F36,F36_grid=grid(F36),F36_strict_below_R023=fr(F36['upper'])<oldlo,root_bracket_lower=str(lo),root_bracket_upper=str(hi),root_bracket_width=str(hi-lo),candidate=str(candidate),F_candidate=Fc,candidate_grid_upper=grid(Fc),candidate_strict_below_R023=fr(Fc['upper'])<oldlo,candidate_strict_below_F36=fr(Fc['upper'])<fr(F36['lower']),R023_old_interval=old['winner']['interval'],log_calls=calls['log'],sqrt_calls=calls['sqrt'],scope='unique minimum only of proved parameter family; F(candidate) is an upper expression, not a lower bound on C7 or a family-minimum enclosure')
  (D/'result.json').write_text(json.dumps(result,indent=2)+'\n');status='completed'
 except BaseException as e:error=repr(e)
 rc=dict(status=status,error=error,start_utc=start,finish_utc=now(),wall_seconds=time.perf_counter()-t,cpu_seconds=time.process_time()-cpu,peak_rss_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,terminal_evaluations=len(out),log_calls=calls['log'],sqrt_calls=calls['sqrt'],python=sys.version,flint=flint.__version__,precision_bits=512);(D/'receipt.json').write_text(json.dumps(rc,indent=2)+'\n');print(json.dumps(rc,indent=2));sys.exit(0 if status=='completed' else 1)
if __name__=='__main__':main()
