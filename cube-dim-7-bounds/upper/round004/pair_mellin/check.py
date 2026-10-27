#!/usr/bin/env python3
"""Validate hashes, 512 terminal panels and independently reconstruct certified scalar intervals; no coefficient recomputation."""
from flint import arb,fmpq,ctx
from fractions import Fraction
from pathlib import Path
import json,hashlib
ctx.prec=512
p=Path(__file__).resolve().parent;mfest=json.loads((p/'manifest.json').read_text())
for f,h in mfest['files'].items():assert hashlib.sha256((p/f).read_bytes()).hexdigest()==h,f
P=[json.loads(x) for x in (p/'panels.jsonl').read_text().splitlines()];assert len(P)==512 and {x['panel'] for x in P}==set(range(512)) and all(x['status']=='completed' for x in P)
def dy(x):return arb(int(x[0]))*arb(2)**x[1]
def frac(x):return Fraction(int(x[0]))*Fraction(2)**int(x[1])
def box(d):return dy(d['lower']).union(dy(d['upper']))
def m(x):return ((-4*x/7)*arb(2).log()).exp()*(4+x).gamma()*(4+8*x/7).rgamma()*((arb(fmpq(1,2))+x/7).gamma()*arb(fmpq(1,2)).rgamma())**8
c=json.loads((p/'certificate.json').read_text());a=arb(fmpq(1,16));r=arb(fmpq(3,64));s=arb(fmpq(1,8));norm=arb.pi()*s.gamma();lo=arb(0);hi=arb(0)
for x in P:
 l=x['integral']['lower'];u=x['integral']['upper'];assert frac(l)>0 and frac(l)<=frac(u);lo+=dy(l);hi+=dy(u)
err=32*(a-r).gamma()*(a+r).gamma()*m(a-r)*m(a+r)*(arb(fmpq(2,3))**129)*3
tail=m(a)**2*(a+16).gamma()**2/(norm*31*arb(32)**31)
qlo=((lo-err)/norm).lower();qhi=((hi+err)/norm+tail).upper();assert frac(c['q']['lower'])<=frac(qlo.man_exp())
assert box(c['q']).contains(qlo) and box(c['q']).contains(qhi)
q=box(c['q']);M=q**4/2;I=(M.log()*arb(fmpq(7,6))).exp()/(105*(arb(2).log()*arb(fmpq(25,6))).exp()*(arb.pi().log()*arb(fmpq(13,12))).exp())
assert box(c['M_pair']).overlaps(M) and box(c['I_pair']).overlaps(I)
assert frac(c['I_pair']['upper'])<Fraction(c['rational_I_upper_numerator'],c['rational_I_upper_denominator'])<frac(c['I_H_R002']['lower'])
print('PASS: manifest, 512 terminal positive panel intervals, independently assembled error/tail, scalar endpoint and strict improvement. Mathematical proof depends on the cited theory; panel coefficients are not recomputed.')
