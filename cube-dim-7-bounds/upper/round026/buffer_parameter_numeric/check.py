#!/usr/bin/env python3
from fractions import Fraction as Q
from pathlib import Path
from flint import arb, fmpq, ctx
import hashlib,json,time,resource
ctx.prec=512;D=Path(__file__).resolve().parent;start=time.perf_counter();cpu=time.process_time()
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
for line in (D/'prerun-sha256.txt').read_text().splitlines():
 h,p=line.split('  ',1);assert sha(D/p)==h
s=D/'sources';mf=json.loads((s/'manifest.json').read_text())
for p in ('result.json','terminal.jsonl'):assert sha(s/p)==mf['files'][p]
old=json.loads((s/'result.json').read_text());oldrows=[json.loads(x) for x in (s/'terminal.jsonl').read_text().splitlines()];o=oldrows[1];assert o['id']==1 and o['status']=='completed' and o['t']=='11/2' and old['status']=='completed' and old['winner']['t']=='11/2' and o['interval']==old['winner']['interval'] and o['grid_upper']==old['winner']['grid_upper']
def a(x):return arb(fmpq(x.numerator,x.denominator))
def fr(v):return Q(int(v[0]))*Q(2)**int(v[1])
def verify(iv,v):
 l=fr(iv['lower']);u=fr(iv['upper']);lv,le=v.lower().man_exp();uv,ue=v.upper().man_exp();vl=Q(int(lv))*Q(2)**int(le);vu=Q(int(uv))*Q(2)**int(ue);assert l<=vl<=vu<=u
r=json.loads((D/'result.json').read_text());rows=[json.loads(x) for x in (D/'terminal.jsonl').read_text().splitlines()];assert len(rows)==83 and [x['id'] for x in rows]==list(range(83)) and all(x['status']=='completed' for x in rows)
lo=Q(5);hi=Q(6)
for i in range(82):
 row=rows[i];x=(lo if i==0 else hi if i==1 else (lo+hi)/2);assert Q(row['t'])==x
 v=a(x**7/(x+1)**5-1)-7*a(x).log();verify(row['interval'],v);sgn=-1 if v<0 else 1 if v>0 else 0;assert sgn==row['sign'] and sgn!=0
 if i<2:assert sgn==(-1 if i==0 else 1)
 else:
  assert row['step']==i-1
  if sgn<0:lo=x
  else:hi=x
  assert Q(row['bracket_lower'])==lo and Q(row['bracket_upper'])==hi
assert hi-lo==Q(1,2**80) and Q(r['root_bracket_lower'])==lo and Q(r['root_bracket_upper'])==hi
x=(lo+hi)/2;assert Q(r['candidate'])==x and Q(rows[82]['t'])==x
rr=a(1+1/x);ll=1+7*a(x).log();v=(rr**3*rr.sqrt()+ll.sqrt())**2;verify(rows[82]['interval'],v);assert rows[82]['interval']==r['candidate_interval'];u=fr(r['candidate_interval']['upper']);g=Q(u.numerator*10**12//u.denominator+1,10**12);assert Q(r['candidate_grid_upper'])==g and Q(rows[82]['grid_upper'])==g and u<g and u<fr(o['interval']['lower']) and r['candidate_strictly_below_R025'] is True
print(json.dumps({'status':'PASS','scientific_rows_rechecked':83,'verification_log_calls':83,'verification_sqrt_calls':2,'cpu_seconds':time.process_time()-cpu,'wall_seconds':time.perf_counter()-start,'peak_rss_bytes':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss},indent=2))
