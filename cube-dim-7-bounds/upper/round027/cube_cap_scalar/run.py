#!/usr/bin/env python3
from flint import arb,fmpq,ctx
from fractions import Fraction as Q
from pathlib import Path
from math import comb
import hashlib,json,time,datetime,resource,sys,flint
ctx.prec=512;D=Path(__file__).resolve().parent
now=lambda:datetime.datetime.now(datetime.timezone.utc).isoformat()
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def a(q):return arb(fmpq(q.numerator,q.denominator)) if isinstance(q,Q) else arb(q)
def ep(x):
 l=x.lower().man_exp();u=x.upper().man_exp();return {'lower':[str(l[0]),int(l[1])],'upper':[str(u[0]),int(u[1])],'display':x.str(45)}
def fr(p):return Q(int(p[0]))*Q(2)**int(p[1])
def grid(iv):
 u=fr(iv['upper']);n=u.numerator*10**12//u.denominator+1;assert u<Q(n,10**12);return str(Q(n,10**12))
def compute(t):
 d=1+1/(2*t);h=(1+1/t)/2
 terms=[Q(comb(6,k))*d**(6-k)*Q(-1,2)**k*(1-t**(-k))/k for k in range(1,7)]
 coeff=Q(7)-Q(7,2)*d**6;constant=1-h**7-Q(7,2)*sum(terms);L=a(t).log();B=a(constant)+a(coeff)*L;B0=1+7*L;assert B>0 and B0-B>0
 z=a((1+1/t)**7);assert z>0;F=(z.sqrt()+B.sqrt())**2;assert F>0 and F.is_finite();return coeff,constant,terms,B,B0,F
def main():
 assert not (D/'terminal.jsonl').exists();start=now();wall=time.perf_counter();cpu=time.process_time();out=[];status='failed';err=None;calls={'rows_started':0,'log':0,'sqrt':0}
 def rec(o):
  out.append(o)
  with (D/'terminal.jsonl').open('a') as f:f.write(json.dumps(o)+'\n');f.flush()
 try:
  s=D/'sources';assert sha(s/'author.md')=='6ecbf72fa9213677dc99e2e2afaa510df240b3d1b08047366a8d0103d010c80f';assert sha(s/'audit.md')=='00555fc5ca9dff8cdb3e467ef0bf4230b5eed8b603a744256f9a5fb66224db6a';mf=json.loads((s/'manifest.json').read_text());assert sha(s/'result.json')==mf['files']['result.json'];old=json.loads((s/'result.json').read_text());assert old['status']=='completed';oldlo=fr(old['candidate_interval']['lower'])
  for i,t in enumerate((Q(4),Q(5),Q(11,2),Q(6),Q(8))):
   calls['rows_started']+=1;calls['log']+=1
   try:coef,const,terms,B,B0,F=compute(t);calls['sqrt']+=2
   except BaseException as e:rec(dict(id=i,t=str(t),status='failed',error=repr(e)));raise
   iv=ep(F);rec(dict(id=i,t=str(t),status='completed',log_coefficient=str(coef),rational_constant=str(const),integral_terms=[str(x) for x in terms],Bcap=ep(B),B0_minus_Bcap=ep(B0-B),F=iv,grid_upper=grid(iv),strict_below_R026=fr(iv['upper'])<oldlo))
  win=min(out,key=lambda r:fr(r['F']['upper']));assert all(fr(win['F']['upper'])<fr(r['F']['lower']) for r in out if r is not win)
  result=dict(status='completed',precision_bits=512,scientific_rows=5,winner=win,strict_below_R026=win['strict_below_R026'],R026_old_candidate_interval=old['candidate_interval'],calls=calls,scope='best of five registered t values; not a continuous optimum or lower bound on C7')
  (D/'result.json').write_text(json.dumps(result,indent=2)+'\n');status='completed'
 except BaseException as e:err=repr(e)
 rc=dict(status=status,error=err,start_utc=start,finish_utc=now(),wall_seconds=time.perf_counter()-wall,cpu_seconds=time.process_time()-cpu,peak_rss_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,terminal_rows=len(out),calls=calls,python=sys.version,flint=flint.__version__,precision_bits=512)
 (D/'receipt.json').write_text(json.dumps(rc,indent=2)+'\n');print(json.dumps(rc,indent=2));sys.exit(0 if status=='completed' else 1)
if __name__=='__main__':main()
