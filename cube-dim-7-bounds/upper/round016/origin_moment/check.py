#!/usr/bin/env python3
"""Recheck the same ten scalar inputs with Arb and exact rational comparisons; no phases or samples."""
from pathlib import Path
from fractions import Fraction
from flint import arb,fmpq,ctx
import json,hashlib,math
ctx.prec=512;p=Path(__file__).resolve().parent;sha=lambda f:hashlib.sha256(f.read_bytes()).hexdigest();mf=json.loads((p/'manifest.json').read_text())
for f,h in mf['files'].items():assert sha(p/f)==h,f
for line in (p/'prerun-sha256.txt').read_text().splitlines():
 h,f=line.split('  ',1);q=Path(f);assert sha(q if q.is_absolute() else p/q)==h,f
for f,m,k in [('r010-cells.jsonl','r010-manifest.json','cells.jsonl'),('r006-certificate.json','r006-manifest.json','certificate.json')]:assert sha(p/f)==json.loads((p/m).read_text())['files'][k]
cells=[json.loads(x) for x in (p/'r010-cells.jsonl').read_text().splitlines()];seen=set();S=Fraction()
for c in cells:
 for x in c['simplices']:
  key=tuple(x['vertices']);assert key not in seen;seen.add(key);q=Fraction(abs(x['detW']),math.prod(x['q']));assert q==Fraction(x['alpha']);S+=q
R=json.loads((p/'result.json').read_text());assert R['status']=='completed' and len(seen)==212 and str(S)==R['A_sum'] and 0<5040*S<=Fraction(1,512) and len(R['enclosures'])==10
def A(x):return arb(fmpq(x.numerator,x.denominator)) if isinstance(x,Fraction) else arb(x)
def dy(x):return Fraction(int(x[0]))*Fraction(2)**int(x[1])
def ball(d):
 def v(x):return A(dy(x))
 return v(d['lower']).union(v(d['upper']))
E=R['enclosures'];pi=arb.pi();p0=5040*A(S)/pi**7;sharp=6/(pi**4*(32*arb(2).sqrt())**7);M=2*pi.sqrt()*(p0.log()/7).exp();raw=(5040*A(S)).log()/6;Bfan=raw.exp()/(840*(pi.log()*arb(fmpq(5,3))).exp());Bsharp=(sharp.log()/6).exp()/(840*pi.sqrt());ratio=arb(fmpq(65,2)).gamma()/arb(32).gamma();b=arb(fmpq(65,2))*arb(fmpq(67,2))*arb(fmpq(69,2));eta=ratio/(b.log()/6).exp();BRB=eta*Bfan
for name,v in [('p0_fan',p0),('p0_sharp_upper',sharp),('M_raw_bound',M),('B_fan',Bfan),('B_sharp',Bsharp),('eta',eta),('B_RB',BRB)]:assert v.overlaps(ball(E[name])),name
via=(M.log()*arb(fmpq(7,6))).exp()/(105*(arb(2).log()*arb(fmpq(25,6))).exp()*(pi.log()*arb(fmpq(13,12))).exp());assert via.overlaps(Bfan) and (arb(fmpq(71,2)).gamma()/arb(32).gamma()).overlaps(ratio*b) and 0<eta<1 and p0<sharp
rawup=Fraction(R['B_fan_rational_upper']);rbup=Fraction(R['B_RB_rational_upper']);assert dy(E['B_fan']['upper'])<rawup and dy(E['B_RB']['upper'])<rbup<rawup
old=json.loads((p/'r006-certificate.json').read_text());assert rawup<dy(old['I_upper']['lower']);assert dy(E['B_sharp']['upper'])<Fraction(4940777144766,10**18)
g=Fraction(1,320000);req=7*(A(rbup/128).log()*arb(fmpq(6,7))).exp();assert req.overlaps(ball(E['g_required_B_RB'])) and req<A(g)
for name,B,k in [('conditional_cost_B_fan',rawup,'conditional_cost6_raw_exact'),('conditional_cost_B_RB',rbup,'conditional_cost6_RB_exact')]:
 c6=Fraction(7**7)*B**6/g**7;assert c6==Fraction(R[k]) and c6<128**6;cost=(A(c6).log()/6).exp();assert cost.overlaps(ball(E[name]))
assert not R['global_g_proved'] and not R['new_C7_bound']
print('PASS: 212 exact alpha sum, ten same-input 512-bit scalar rechecks, raw/Jensen normalization, gamma/RB eta, strict rational budgets and conditional sixth-power cost. No phase, sample or global g claim.')
