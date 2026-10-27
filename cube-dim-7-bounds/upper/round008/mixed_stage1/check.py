#!/usr/bin/env python3
"""Check signed panel accounting and reassemble scalar bounds; no panel rerun."""
from flint import arb,fmpq,ctx
from fractions import Fraction
from pathlib import Path
import json,hashlib
ctx.prec=512;p=Path(__file__).resolve().parent;sha=lambda f:hashlib.sha256((p/f).read_bytes()).hexdigest();m=json.loads((p/'manifest.json').read_text())
for f,h in m['files'].items():assert sha(f)==h,f
P=[json.loads(x) for x in (p/'panels.jsonl').read_text().splitlines()];assert len(P)==256 and {(x['group'],x['panel']) for x in P}=={(j,k) for j in range(4) for k in range(64)} and all(x['status']=='completed' for x in P)
assert sha('r004-certificate.json')==json.loads((p/'r004-manifest.json').read_text())['files']['certificate.json'];old=json.loads((p/'r004-certificate.json').read_text());C=json.loads((p/'certificate.json').read_text())
def dy(d):return arb(int(d[0]))*arb(2)**int(d[1])
def fr(d):return Fraction(int(d[0]))*Fraction(2)**int(d[1])
def F(z,j):
 e=arb(fmpq(1,7));return ((-4*z/7-2*j*e)*arb(2).log()).exp()*(4+z-j*e).gamma()*(4+8*z/7-j*e).rgamma()*(arb(fmpq(1,2))+z/7-e).gamma()**j*(arb(fmpq(1,2))+z/7).gamma()**(8-j)*arb(fmpq(1,2)).rgamma()**8
for row in P:
 terms=row['signed_integrated_even_terms'];assert len(terms)==65;sl=sum((fr(t['lower']) for t in terms),Fraction());su=sum((fr(t['upper']) for t in terms),Fraction());assert sl<=fr(row['integral']['upper']) and fr(row['integral']['lower'])<=su
assert any(fr(t['upper'])<0 for x in P for t in x['signed_integrated_even_terms'])
assert any(fr(x['constant_imag']['lower'])>0 or fr(x['constant_imag']['upper'])<0 for x in P)
a=arb(fmpq(9,16));r=arb(fmpq(27,64));h=arb(fmpq(1,4));b=arb(fmpq(9,8));norm=arb.pi()*b.gamma();I=[]
for j in range(4):
 lo=sum((dy(x['integral']['lower']) for x in P if x['group']==j),arb(0));hi=sum((dy(x['integral']['upper']) for x in P if x['group']==j),arb(0));c1=(a-r).gamma()*F(a-r,j)*(a+r).gamma()*F(a+r,7-j);c2=(a+r).gamma()*F(a+r,j)*(a-r).gamma()*F(a-r,7-j);u1=c1.upper();u2=c2.upper();cs=u1 if u1>u2 else u2
 err=32*cs*(h/r)**129/(1-h/r);tail=F(a,j)*F(a,7-j)*(a+16).gamma()**2/(norm*31*arb(32)**31);l=((lo-err)/norm-tail).lower();u=((hi+err)/norm+tail).upper()
 if l<0:l=arb(0)
 row=C['mixed_expectations'][j];assert fr(row['expectation_lower']['lower'])<=fr(l.man_exp()) and fr(row['expectation_upper']['upper'])>=fr(u.man_exp());I.append((dy(row['expectation_lower']['lower']),dy(row['expectation_upper']['upper'])))
co=[1,7,21,35];sl=2*sum((co[j]*I[j][0]**7 for j in range(4)),arb(0));su=2*sum((co[j]*I[j][1]**7 for j in range(4)),arb(0));qlo=dy(old['q']['lower']);qhi=dy(old['q']['upper']);LHlo=(8*(sl.log()/7).exp()/qhi).lower();LHhi=(8*(su.log()/7).exp()/qlo).upper();assert fr(C['LH_lower']['lower'])<=fr(LHlo.man_exp()) and fr(C['LH_upper']['upper'])>=fr(LHhi.man_exp())
assert LHhi<1 and C['status']=='mu_not_positive' and C['mu_positive'] is False and C['stage_II_authorized'] is False
formal=(dy(old['M_pair']['lower'])/LHhi.sqrt()).lower();assert formal>arb(fmpq(31,500))
print('PASS: 256 terminal signed panels, all 65 terms per panel, nonzero imaginary kernel, crossed Cauchy bounds, two-sided tails and LH<1; Stage II remains unauthorized. The formal negative-mu expression is not a feasible Radau bound.')
