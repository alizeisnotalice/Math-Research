#!/usr/bin/env python3
"""Independent one-point rational audit, no original math modules imported."""
from fractions import Fraction as Q
from pathlib import Path
from math import comb,isqrt
import json,hashlib,time,resource,sys,datetime,traceback
sys.set_int_max_str_digits(0)
D=Path(__file__).resolve().parent

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def check_sources(lock):
 root=Path(lock['root']);assert all(sha(root/n)==h for n,h in lock['files'].items())
def ep(pair):return Q(int(pair[0]))*Q(2)**int(pair[1])
def sqrtenclose(a):
 assert a>0;power=700;den=1<<power;n=isqrt((a.numerator<<(2*power))//a.denominator)
 assert n*n*a.denominator<=a.numerator*(den*den)<(n+1)**2*a.denominator
 return Q(n,den),Q(n+1,den)
def outward(x,direction,places=80):
 d=10**places;k=x.numerator*d//x.denominator
 return Q(k+(direction=='upper'),d)
def main():
 w=time.perf_counter();c=time.process_time();start=datetime.datetime.now(datetime.timezone.utc).isoformat();lock=json.loads((D/'sourcehash.json').read_text());check_sources(lock);root=Path(lock['root']);new=root/'upper/round028/cap_candidate';old=root/'upper/round027/cube_cap_scalar';lower=root/'lower/round006';audit=root/'audit/round006/lower'
 result=json.loads((new/'result.json').read_text());receipt=json.loads((new/'receipt.json').read_text());trace=[json.loads(s) for s in (new/'trace.jsonl').read_text().splitlines()];cache=json.loads((new/'cache.json').read_text());prior=json.loads((old/'result.json').read_text());rows=[json.loads(s) for s in (old/'terminal.jsonl').read_text().splitlines()]
 assert (new/'sources/result.json').read_bytes()==(old/'result.json').read_bytes() and (new/'sources/terminal.jsonl').read_bytes()==(old/'terminal.jsonl').read_bytes()
 assert len(rows)==len(cache)==5 and prior['winner']==rows[3] and rows[3]['t']=='6'
 assert all(Q(r['t'])==Q(k['t']) and float(Q(r['t'])).hex()==k['float_t_hex'] and r['F']==k['interval'] and float((ep(r['F']['lower'])+ep(r['F']['upper']))/2).hex()==k['midpoint_float_hex'] for r,k in zip(rows,cache))
 assert len(trace)==receipt['callback_calls']==result['callback_calls']==result['new_double_evaluations']==15 and receipt['cache_hits']==result['cache_hits']==0
 assert len({r['x_hex'] for r in trace})==15 and all(r['call']==i and r['new_eval_count']==i and not r['cache_hit'] for i,r in enumerate(trace,1))
 assert not {r['x_hex'] for r in trace}&{k['float_t_hex'] for k in cache}
 raw=result['raw_t_hex'];matches=[r for r in trace if r['x_hex']==raw];assert len(matches)==1 and matches[0]['y_hex']==result['raw_value_hex']
 rawq=Q.from_float(float.fromhex(raw));integer=round(rawq*2**20);t=Q(integer,2**20)
 assert integer==result['integer_q']==6783498 and result['denominator']==2**20 and t==Q(result['exact_t'])==Q(3391749,524288)
 ledger=json.loads((root/'results/bounds_ledger.json').read_text());L=ledger['lower'];lr=json.loads((lower/'results.json').read_text());la=json.loads((audit/'exactresults.json').read_text())
 assert 'RESEARCH_PAUSE_REQUESTED_AFTER_R028' in (root/'STATE.md').read_text() and ledger['research_pause']['requested_by_user'] and ledger['research_pause']['lower_remains_paused'] and not ledger['research_pause']['new_research_authorized']
 assert L['certificate_sha256']==sha(lower/'results.json') and L['exact_interval']==lr['full_construction_comparison_interval'] and Q(L['finite_L1_exact_lower'])==Q(lr['finite_L1']['exact_ratio_lower'])==Q(la['finite_L1_exact_lower']) and Q(L['simple_rational_strict_lower'])==Q(2136578727666,10**12)
 U=Q(4910992105193,200000000000);assert Q(result['candidate_grid_upper'])==U==Q(ledger['upper']['exact'])
 x=(t-1)/(t+1);assert x==Q(2867461,3916037) and 0<x<1
 N=512;term=x;total=Q(0)
 for k in range(N):total+=term/(2*k+1);term*=x*x
 loglo=2*total;loghi=loglo+2*term/((2*N+1)*(1-x*x));assert 0<loglo<loghi
 q=1+1/(2*t);h=(1+1/t)/2
 parts=[Q(comb(6,j))*q**(6-j)*Q(-1,2)**j*(1-t**(-j))/j for j in range(1,7)]
 constant=1-h**7-Q(7,2)*sum(parts,Q(0));coeff=7-Q(7,2)*q**6
 assert str(coeff)==result['log_coefficient'] and str(constant)==result['rational_constant']
 vals=(constant+coeff*loglo,constant+coeff*loghi);B_lo=min(vals);B_hi=max(vals);assert 0<B_lo<=B_hi
 r=(1+1/t)**7;s_lo,s_hi=sqrtenclose(r);b_lo,_=sqrtenclose(B_lo);_,b_hi=sqrtenclose(B_hi)
 F_lo=(s_lo+b_lo)**2;F_hi=(s_hi+b_hi)**2;oldlo=ep(prior['winner']['F']['lower']);oldhi=ep(prior['winner']['F']['upper'])
 assert F_lo<F_hi<U<oldlo<oldhi and result['R027_winner_interval']==prior['winner']['F']
 assert F_lo<ep(result['candidate_interval']['upper']) and F_hi>ep(result['candidate_interval']['lower'])
 check_sources(lock)
 return {'status':'PASS','scope':'independent exact Fraction recertification at one frozen t; no theory or lower geometry review','t':str(t),'x_atanh':str(x),'log_terms':N,'sqrt_dyadic_bits':700,'log_bounds_outward_80_places':[str(outward(loglo,'lower')),str(outward(loghi,'upper'))],'B_bounds_outward_80_places':[str(outward(B_lo,'lower')),str(outward(B_hi,'upper'))],'F_bounds_outward_80_places':[str(outward(F_lo,'lower')),str(outward(F_hi,'upper'))],'U':str(U),'B_strictly_positive':True,'F_upper_strict_below_U':True,'U_strict_below_R027_full_lower':True,'log_coefficient':str(coeff),'rational_constant':str(constant),'trace_unique_count':len(trace),'raw_to_integer_q':integer,'r027_cache_entries_verified':len(cache),'lower_hash_and_ledger_match':True,'research_paused':True,'source_hashes_unchanged':True,'no_optimizer_or_old_grid_replay':True,'started_utc':start,'wall_seconds':time.perf_counter()-w,'cpu_seconds':time.process_time()-c,'peak_rss_bytes':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,'executable':sys.executable}
if __name__=='__main__':
 try:print(json.dumps(main(),indent=2));sys.exit(0)
 except Exception as exc:traceback.print_exc(file=sys.stderr);print(json.dumps({'status':'FAIL','error':repr(exc)}));sys.exit(1)
