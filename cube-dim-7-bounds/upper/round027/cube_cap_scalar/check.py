#!/usr/bin/env python3
from fractions import Fraction as Q
from pathlib import Path
from math import comb
from flint import arb,fmpq,ctx
import hashlib,json,time,resource
ctx.prec=512;D=Path(__file__).resolve().parent;wall=time.perf_counter();cpu=time.process_time()
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
for line in (D/'prerun-sha256.txt').read_text().splitlines():
 h,p=line.split('  ',1);assert sha(D/p)==h
s=D/'sources';mf=json.loads((s/'manifest.json').read_text());assert sha(s/'result.json')==mf['files']['result.json'];old=json.loads((s/'result.json').read_text());assert old['status']=='completed'
def a(x):return arb(fmpq(x.numerator,x.denominator)) if isinstance(x,Q) else arb(x)
def fr(p):return Q(int(p[0]))*Q(2)**int(p[1])
def overlap(iv,v):
 m,e=v.lower().man_exp();l=Q(int(m))*Q(2)**int(e);m,e=v.upper().man_exp();u=Q(int(m))*Q(2)**int(e);assert fr(iv['lower'])<=u and l<=fr(iv['upper'])
rows=[json.loads(x) for x in (D/'terminal.jsonl').read_text().splitlines()];result=json.loads((D/'result.json').read_text());assert len(rows)==5 and [r['id'] for r in rows]==list(range(5)) and all(r['status']=='completed' for r in rows)
for row,t in zip(rows,[Q(4),Q(5),Q(11,2),Q(6),Q(8)]):
 assert Q(row['t'])==t;d=1+1/(2*t);h=(1+1/t)/2;terms=[Q(comb(6,k))*d**(6-k)*Q(-1,2)**k*(1-t**(-k))/k for k in range(1,7)];assert [Q(i) for i in row['integral_terms']]==terms
 coef=7-Q(7,2)*d**6;const=1-h**7-Q(7,2)*sum(terms);assert Q(row['log_coefficient'])==coef and Q(row['rational_constant'])==const
 L=a(t).log();B=a(const)+a(coef)*L;gap=1+7*L-B;assert B>0 and gap>0;overlap(row['Bcap'],B);overlap(row['B0_minus_Bcap'],gap);F=(a((1+1/t)**7).sqrt()+B.sqrt())**2;overlap(row['F'],F);u=fr(row['F']['upper']);assert Q(row['grid_upper'])==Q(u.numerator*10**12//u.denominator+1,10**12);assert row['strict_below_R026']==(u<fr(old['candidate_interval']['lower']))
w=result['winner'];assert w==rows[w['id']] and all(fr(w['F']['upper'])<fr(r['F']['lower']) for r in rows if r['id']!=w['id']);assert result['strict_below_R026'] and w['strict_below_R026']
print(json.dumps({'status':'PASS','verification_rows':5,'verification_log_calls':5,'verification_sqrt_calls':10,'cpu_seconds':time.process_time()-cpu,'wall_seconds':time.perf_counter()-wall,'peak_rss_bytes':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss},indent=2))
