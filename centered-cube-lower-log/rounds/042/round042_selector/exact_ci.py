"""Exact integer certification of conservative rational Bernoulli KL endpoints.
Counts come from stored whole-query-cell envelopes; no claim about PRNG randomness.
"""
from pathlib import Path
from fractions import Fraction as F
import json,math,time
HERE=Path(__file__).resolve().parent
DEN=10**12

def certified(k,N,num,lower):
 assert 0<=num<=DEN
 assert num*N<=k*DEN if lower else num*N>=k*DEN
 # exp(-N*KL(k/N || num/DEN)) <= 1/720, using only integers.
 left=720*pow(num*N,k)*pow((DEN-num)*N,N-k)
 right=pow(DEN*k,k)*pow(DEN*(N-k),N-k)
 assert left<=right,(k,N,num,lower)
 return dict(numerator=num,denominator=DEN,exp_minus_N_KL_le_exact_alpha=True)

def main():
 start=time.time();target=HERE/'EXACT_CI.json';assert not target.exists()
 source=HERE/'CELL_RESULTS.json';d=json.loads(source.read_text());rows=[]
 for case in d['cases']:
  rr=[]
  for row in case['rows']:
   N=row['N'];kl=row['lower_successes'];ku=row['upper_successes']
   # Two rational grid steps outward; then independently prove the Chernoff inequality.
   l=max(0,math.floor(row['weak_CI'][0]*DEN)-2)
   u=min(DEN,math.ceil(row['weak_CI'][1]*DEN)+2)
   lo=certified(kl,N,l,True);hi=certified(ku,N,u,False)
   rr.append(dict(model=row['model'],lower_count=kl,upper_count=ku,N=N,lower=lo,upper=hi))
  priorgap=F(rr[0]['lower']['numerator']-rr[1]['upper']['numerator'],DEN)
  fixedgap=F(rr[0]['lower']['numerator']-max(r['upper']['numerator']for r in rr[2:]),DEN)
  rows.append(dict(n=case['n'],rows=rr,gated_lower_minus_prior_upper_exact=str(priorgap),gated_lower_minus_best_fixed_upper_exact=str(fixedgap)))
 import hashlib
 target.write_text(json.dumps(dict(cases=rows,source_SHA256=hashlib.sha256(source.read_bytes()).hexdigest(),
   script_SHA256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),endpoints=36,alpha_each_exact='1/720',
   total_failure_budget_exact='1/20',integer_Chernoff_comparisons=36,passed=True,
   scope='Exact departure certificates for rational bounds; uniformity via fixed-target monotone KL event. Requires valid cell envelopes and ideal iid53bit queries; neither formal PRNG proof nor deterministic continuum integral certificate.',seconds=time.time()-start),indent=2)+'\n')
 print(json.dumps(dict(passed=True,endpoints=36,seconds=time.time()-start)))
if __name__=='__main__':main()
