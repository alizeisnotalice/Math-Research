#!/usr/bin/env python3
from flint import arb,fmpq,ctx
from fractions import Fraction as Q
from pathlib import Path
import hashlib,json,time,datetime,resource,sys,flint
ctx.prec=512; D=Path(__file__).resolve().parent
now=lambda:datetime.datetime.now(datetime.timezone.utc).isoformat()
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def A(x):return arb(fmpq(x.numerator,x.denominator)) if isinstance(x,Q) else arb(x)
def ep(x):
 l=x.lower().man_exp();u=x.upper().man_exp();return {'lower':[str(l[0]),int(l[1])],'upper':[str(u[0]),int(u[1])],'display':x.str(45)}
def fr(v):return Q(int(v[0]))*Q(2)**int(v[1])
def grid(v):
 u=fr(v['upper']); n=u.numerator*10**12//u.denominator+1;assert u<Q(n,10**12);return str(Q(n,10**12))
def G(x):return A(x**7/(x+1)**5-1)-7*A(x).log()
def F(x):
 r=A(1+1/x);L=1+7*A(x).log(); assert r>0 and L>0
 return (r**3*r.sqrt()+L.sqrt())**2
def main():
 assert not (D/'terminal.jsonl').exists();start=now();wall=time.perf_counter();cpu=time.process_time();out=[];status='failed';err=None;calls={'G':0,'F':0,'log':0,'sqrt':0}
 def rec(o):
  out.append(o)
  with (D/'terminal.jsonl').open('a') as f:f.write(json.dumps(o)+'\n');f.flush()
 try:
  s=D/'sources';assert sha(s/'proof.md')=='322fb89887e6612b83bf04dd29cf1735407fe667de9d06235c24ab1135774a6e';assert sha(s/'parameters.json')=='e0ca97979452e512faacd9f3073f357fdb5574c19f48eb4396b5623d71517c50'
  mf=json.loads((s/'manifest.json').read_text());assert all(sha(s/name)==mf['files'][name] for name in ['result.json','terminal.jsonl']);old=json.loads((s/'result.json').read_text());oldrow=[json.loads(x) for x in (s/'terminal.jsonl').read_text().splitlines() if json.loads(x)['id']==1];assert old['status']=='completed' and old['winner']['t']=='11/2' and len(oldrow)==1 and oldrow[0]['status']=='completed' and oldrow[0]['t']=='11/2' and oldrow[0]['interval']==old['winner']['interval'] and oldrow[0]['grid_upper']==old['winner']['grid_upper'];oldlo=fr(old['winner']['interval']['lower'])
  lo=Q(5);hi=Q(6)
  for i,x in enumerate((lo,hi)):
   calls['G']+=1;calls['log']+=1;v=G(x);sgn=-1 if v<0 else 1 if v>0 else 0;rec(dict(id=i,kind='G_endpoint',t=str(x),interval=ep(v),sign=sgn,status='completed' if sgn==(-1 if i==0 else 1) else 'failed'));assert sgn==(-1 if i==0 else 1)
  for i in range(80):
   mid=(lo+hi)/2;calls['G']+=1;calls['log']+=1;v=G(mid);sgn=-1 if v<0 else 1 if v>0 else 0
   if sgn==0:rec(dict(id=i+2,kind='G_midpoint',step=i+1,t=str(mid),interval=ep(v),sign=0,status='unresolved_sign'));raise ArithmeticError('unresolved sign')
   if sgn<0:lo=mid
   else:hi=mid
   rec(dict(id=i+2,kind='G_midpoint',step=i+1,t=str(mid),interval=ep(v),sign=sgn,bracket_lower=str(lo),bracket_upper=str(hi),status='completed'))
  assert hi-lo==Q(1,2**80);candidate=(lo+hi)/2;calls['F']+=1;calls['log']+=1;calls['sqrt']+=2;v=F(candidate);assert v>0 and v.is_finite();iv=ep(v);rec(dict(id=82,kind='F_candidate',t=str(candidate),interval=iv,grid_upper=grid(iv),status='completed'));assert len(out)==83
  result=dict(status='completed',precision_bits=512,new_function_evaluations=83,root_bracket_lower=str(lo),root_bracket_upper=str(hi),root_bracket_width=str(hi-lo),candidate=str(candidate),candidate_distance_bound='1/2417851639229258349412352',candidate_interval=iv,candidate_grid_upper=grid(iv),reused_R025_row=old['winner'],candidate_strictly_below_R025=fr(iv['upper'])<oldlo,calls=calls,scope='root uniquely minimizes this proved parameter family; candidate interval is not a family-minimum enclosure or C7 lower bound')
  (D/'result.json').write_text(json.dumps(result,indent=2)+'\n');status='completed'
 except BaseException as e:err=repr(e)
 rc=dict(status=status,error=err,start_utc=start,finish_utc=now(),wall_seconds=time.perf_counter()-wall,cpu_seconds=time.process_time()-cpu,peak_rss_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,terminal_rows=len(out),calls=calls,python=sys.version,flint=flint.__version__,precision_bits=512)
 (D/'receipt.json').write_text(json.dumps(rc,indent=2)+'\n');print(json.dumps(rc,indent=2));sys.exit(0 if status=='completed' else 1)
if __name__=='__main__':main()
