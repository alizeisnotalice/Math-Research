#!/usr/bin/env python3
from fractions import Fraction as Q
from pathlib import Path
from math import comb
from flint import arb,fmpq,ctx
import hashlib,json,time,resource
ctx.prec=512;D=Path(__file__).resolve().parent;w=time.perf_counter();c=time.process_time()
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
for line in (D/'prerun-sha256.txt').read_text().splitlines():
 h,p=line.split('  ',1);assert sha(D/p)==h
s=D/'sources';mf=json.loads((s/'manifest.json').read_text());assert all(sha(s/p)==mf['files'][p] for p in ('result.json','terminal.jsonl'));old=json.loads((s/'result.json').read_text());rows=[json.loads(x) for x in (s/'terminal.jsonl').read_text().splitlines()];r=json.loads((D/'result.json').read_text());trace=[json.loads(x) for x in (D/'trace.jsonl').read_text().splitlines()];cache=json.loads((D/'cache.json').read_text())
def fr(v):return Q(int(v[0]))*Q(2)**int(v[1])
def a(x):return arb(fmpq(x.numerator,x.denominator)) if isinstance(x,Q) else arb(x)
for c0,o in zip(cache,rows):
 assert c0['t']==o['t'] and c0['interval']==o['F'] and float.fromhex(c0['float_t_hex'])==float(Q(o['t'])) and float.fromhex(c0['midpoint_float_hex'])==float((fr(o['F']['lower'])+fr(o['F']['upper']))/2)
assert len(cache)==5 and len(trace)==r['callback_calls']==r['opt_nfev'] and r['new_double_evaluations']<=64 and r['new_double_evaluations']+r['cache_hits']==len(trace)
for i,z in enumerate(trace):
 assert z['call']==i+1 and 4<=float.fromhex(z['x_hex'])<=10 and z['cache_hit']==(z['x_hex'] in [v['float_t_hex'] for v in cache]);assert z['new_eval_count']==sum(not x['cache_hit'] for x in trace[:i+1]);float.fromhex(z['y_hex'])
raw=float.fromhex(r['raw_t_hex']);assert any(z['x_hex']==r['raw_t_hex'] and z['y_hex']==r['raw_value_hex'] for z in trace);assert round(Q.from_float(raw)*2**20)==r['integer_q'];t=Q(r['integer_q'],2**20);assert Q(r['exact_t'])==t and r['new_arb_rows']==1
d=1+1/(2*t);h=(1+1/t)/2;terms=[Q(comb(6,k))*d**(6-k)*Q(-1,2)**k*(1-t**(-k))/k for k in range(1,7)];coef=7-Q(7,2)*d**6;const=1-h**7-Q(7,2)*sum(terms);assert Q(r['log_coefficient'])==coef and Q(r['rational_constant'])==const
L=a(t).log();B=a(const)+a(coef)*L;assert B>0 and 1+7*L-B>0;F=(a((1+1/t)**7).sqrt()+B.sqrt())**2
l,e=F.lower().man_exp();u,k=F.upper().man_exp();lo=Q(int(l))*Q(2)**int(e);hi=Q(int(u))*Q(2)**int(k);iv=r['candidate_interval'];assert fr(iv['lower'])<=lo<=hi<=fr(iv['upper']);upper=fr(iv['upper']);assert Q(r['candidate_grid_upper'])==Q(upper.numerator*10**12//upper.denominator+1,10**12);assert upper<fr(old['winner']['F']['lower']) and r['strict_below_R027']
print(json.dumps({'status':'PASS','verified_arb_rows':1,'verification_log_calls':1,'verification_sqrt_calls':2,'exploratory_points_recomputed':0,'cpu_seconds':time.process_time()-c,'wall_seconds':time.perf_counter()-w,'peak_rss_bytes':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss},indent=2))
