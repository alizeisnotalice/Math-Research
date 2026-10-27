#!/usr/bin/env python3
"""Recheck only the same four scalar inputs, their ceilings and the old interval."""
from flint import arb,fmpq,ctx
from fractions import Fraction
from pathlib import Path
import hashlib,json
ctx.prec=512;p=Path(__file__).resolve().parent;sha=lambda f:hashlib.sha256(f.read_bytes()).hexdigest();mf=json.loads((p/'manifest.json').read_text())
for f,h in mf['files'].items():assert sha(p/f)==h,f
for line in (p/'prerun-sha256.txt').read_text().splitlines():
 h,f=line.split('  ',1);q=Path(f);assert sha(q if q.is_absolute() else p/q)==h,f
assert sha(p/'r024-result.json')==json.loads((p/'r024-manifest.json').read_text())['files']['result.json'];old=json.loads((p/'r024-result.json').read_text());r=json.loads((p/'result.json').read_text());T=[json.loads(x) for x in (p/'terminal.jsonl').read_text().splitlines()];assert len(T)==4 and [x['id'] for x in T]==[0,1,2,3] and [x['t'] for x in T]==['5','11/2','6','5'] and all(x['status']=='completed' for x in T)
def A(v):return arb(fmpq(v.numerator,v.denominator))
def fr(v):return Fraction(int(v[0]))*Fraction(2)**int(v[1])
def box(v):return A(fr(v['lower'])).union(A(fr(v['upper'])))
logs={}
for row in T[:3]:
 t=Fraction(row['t']);z=A(t).log();logs[t]=z;aa=A(1+1/t)**7;bb=1+7*z;v=aa+bb+2*(aa*bb).sqrt();assert v.overlaps(box(row['interval']))
 u=fr(row['interval']['upper']);N=10**12;assert Fraction(u.numerator*N//u.denominator+1,N)==Fraction(row['grid_upper']) and u<Fraction(row['grid_upper'])
t=Fraction(5);tau=Fraction(20,3);aa=A(1+1/t)**7;bb=1+7*logs[t];v=(1+aa/A(tau))*(A(tau)+bb);assert v.overlaps(box(T[3]['interval'])) and fr(T[0]['interval']['upper'])<fr(T[3]['interval']['lower']) and fr(T[3]['interval']['upper'])<Fraction(2233,75)
w=next(x for x in T[:3] if x['t']==r['winner']['t']);assert w['t']=='11/2' and w['interval']==r['winner']['interval'] and all(fr(w['interval']['upper'])<fr(x['interval']['lower']) for x in T[:3] if x is not w)
assert Fraction(r['winner']['grid_upper'])==min(Fraction(x['grid_upper']) for x in T[:3])==Fraction(29059800754803,10**12) and fr(w['interval']['upper'])<fr(old['F_candidate']['lower']) and r['R024_candidate_interval']==old['F_candidate'] and not r['continuous_t_optimum_claim']
print('PASS: three same-input optimized rows, one rational-tau calibration, strict grid/unique finite winner, and full R024 interval comparison. No extra parameter or continuous optimum claim.')
