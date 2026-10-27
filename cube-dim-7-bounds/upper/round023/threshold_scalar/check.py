#!/usr/bin/env python3
"""Recheck the same eight Arb rows and exact ceilings, no additional parameters."""
from flint import arb,fmpq,ctx
from fractions import Fraction
from pathlib import Path
import hashlib,json
ctx.prec=512;p=Path(__file__).resolve().parent;sha=lambda q:hashlib.sha256(q.read_bytes()).hexdigest();mf=json.loads((p/'manifest.json').read_text())
for f,h in mf['files'].items():assert sha(p/f)==h,f
for line in (p/'prerun-sha256.txt').read_text().splitlines():
 h,f=line.split('  ',1);q=Path(f);assert sha(q if q.is_absolute() else p/q)==h,f
assert sha(p/'r022-result.json')==json.loads((p/'r022-manifest.json').read_text())['files']['result.json'];old=json.loads((p/'r022-result.json').read_text());assert old['status']=='completed' and old['winner_a']==63;oldlo=Fraction(old['winner_lower_exact']);oldhi=Fraction(old['winner_upper_exact']);assert oldlo<oldhi
T=[json.loads(x) for x in (p/'terminal.jsonl').read_text().splitlines()];R=json.loads((p/'result.json').read_text());P=json.loads((p/'parameters.json').read_text());assert R['status']=='completed' and len(T)==8 and [[x['a'],x['tau']] for x in T]==P['pair_order']
def fr(x):return Fraction(int(x[0]))*Fraction(2)**int(x[1])
def A(x):return arb(fmpq(x.numerator,x.denominator))
def ball(z):return A(fr(z['lower'])).union(A(fr(z['upper'])))
for x in T:
 a=x['a'];tau=x['tau'];alpha=arb(fmpq(a+2,a))**7;beta=arb(fmpq(a+1,a))**7*(1+7*arb(a).log());K=(1+alpha/tau)*(tau+beta);alt=alpha+beta+tau+alpha*beta/tau;assert K>0 and K.overlaps(alt) and K.overlaps(ball(x['interval']))
 u=fr(x['interval']['upper']);N=10**12;ceiling=Fraction(u.numerator*N//u.denominator+1,N);assert ceiling==Fraction(x['grid_upper']) and u<ceiling and x['strictly_below_R022']==(u<oldlo)
W=R['winner'];w=next(x for x in T if x['a']==W['a'] and x['tau']==W['tau']);assert W['interval']==w['interval'] and (W['a'],W['tau'])==(36,7)
assert all(fr(w['interval']['upper'])<fr(x['interval']['lower']) for x in T if x is not w) and Fraction(R['grid_upper'])==min(Fraction(x['grid_upper']) for x in T) and Fraction(R['grid_upper'])==Fraction(46650446233303,10**12) and fr(w['interval']['upper'])<oldlo
assert fr(next(x for x in T if (x['a'],x['tau'])==(64,6))['interval']['upper'])<Fraction(1040858,21375)<49
assert R['strictly_below_R022_full_interval'] and R['winner_strict_best_of_registered_8'] and not R['continuous_optimum_claim']
print('PASS: eight same-input 512-bit scalar checks, exact per-row ceilings, unique finite-grid winner, analytic tau=6 bound and full R022 interval comparison. Operator theorem depends on frozen proof; no continuous optimization.')
