#!/usr/bin/env python3
"""Recheck same 84 fixed function inputs and exact bracket, no new parameter row."""
from flint import arb,fmpq,ctx
from fractions import Fraction
from pathlib import Path
import hashlib,json
ctx.prec=512;p=Path(__file__).resolve().parent;sha=lambda q:hashlib.sha256(q.read_bytes()).hexdigest();m=json.loads((p/'manifest.json').read_text())
for f,h in m['files'].items():assert sha(p/f)==h,f
for line in (p/'prerun-sha256.txt').read_text().splitlines():
 h,f=line.split('  ',1);q=Path(f);assert sha(q if q.is_absolute() else p/q)==h,f
assert sha(p/'r023-result.json')==json.loads((p/'r023-manifest.json').read_text())['files']['result.json'];old=json.loads((p/'r023-result.json').read_text());R=json.loads((p/'result.json').read_text());T=[json.loads(x) for x in (p/'terminal.jsonl').read_text().splitlines()];assert len(T)==84 and [x['id'] for x in T]==list(range(84)) and all(x['status']=='completed' for x in T)
def A(x):return arb(fmpq(x.numerator,x.denominator))
def fr(x):return Fraction(int(x[0]))*Fraction(2)**int(x[1])
def ball(x):return A(fr(x['lower'])).union(A(fr(x['upper'])))
def F(x):
 aa=A(1+2/x)**7;bb=A(1+1/x)**7*(1+7*A(x).log());return aa+bb+2*(aa*bb).sqrt()
def S(x):
 r=A((x+1)/(x+2));z=A(x).log();return r*r*r.sqrt()*(A(x)-7*z)/(2*(1+7*z).sqrt())-1
assert F(Fraction(36)).overlaps(ball(T[0]['interval'])) and T[0]['kind']=='F_fixed';lo=Fraction(35);hi=Fraction(37);vlo=S(lo);vhi=S(hi);assert vlo<0 and vhi>0 and vlo.overlaps(ball(T[1]['interval'])) and vhi.overlaps(ball(T[2]['interval']))
for i,row in enumerate(T[3:83],1):
 mid=(lo+hi)/2;assert row['a']==str(mid) and row['step']==i;v=S(mid);sgn=-1 if v<0 else 1 if v>0 else 0;assert sgn!=0 and sgn==row['sign'] and v.overlaps(ball(row['interval']))
 if sgn<0:lo=mid
 else:hi=mid
 assert row['bracket_lower']==str(lo) and row['bracket_upper']==str(hi)
assert hi-lo==Fraction(1,2**79) and str(lo)==R['root_bracket_lower'] and str(hi)==R['root_bracket_upper'];candidate=(lo+hi)/2;assert str(candidate)==R['candidate']==T[83]['a'];assert F(candidate).overlaps(ball(T[83]['interval']))
oldlo=fr(old['winner']['interval']['lower']);assert R['R023_old_interval']==old['winner']['interval']
for row in (T[0],T[83]):
 u=fr(row['interval']['upper']);grid=Fraction((u.numerator*10**12)//u.denominator+1,10**12);assert Fraction(row['grid_upper'])==grid and u<grid and row['strictly_below_R023']==(u<oldlo)
assert R['F36']==T[0]['interval'] and R['F_candidate']==T[83]['interval'] and Fraction(R['candidate_grid_upper'])==Fraction(46644289559722,10**12) and R['candidate_strict_below_R023'] and R['candidate_strict_below_F36'] and fr(T[83]['interval']['upper'])<fr(T[0]['interval']['lower'])
print('PASS: same 84 fixed inputs, strict sign/bracket sequence, dyadic midpoint and both ceilings, full R023 comparison. The family minimizer bracket is not a C7 lower bound or minimum-value enclosure.')
