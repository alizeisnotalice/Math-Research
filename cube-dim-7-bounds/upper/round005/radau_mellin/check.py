#!/usr/bin/env python3
"""Reassemble scalar certificate from saved panel endpoints; does not rerun Taylor panels."""
from flint import arb,fmpq,ctx
from fractions import Fraction
from pathlib import Path
import json,hashlib
ctx.prec=512;p=Path(__file__).resolve().parent;h=lambda f:hashlib.sha256((p/f).read_bytes()).hexdigest();mf=json.loads((p/'manifest.json').read_text())
for f,d in mf['files'].items():assert h(f)==d,f
P=[json.loads(x) for x in (p/'panels.jsonl').read_text().splitlines()];N=[1024,512,256];assert len(P)==1792 and {(x['group'],x['panel']) for x in P}=={(g,j) for g,n in enumerate(N) for j in range(n)} and all(x['status']=='completed' for x in P)
old=json.loads((p/'r004-certificate.json').read_text());assert h('r004-certificate.json')==json.loads((p/'r004-manifest.json').read_text())['files']['certificate.json'];C=json.loads((p/'certificate.json').read_text())
def dy(d):return arb(int(d[0]))*arb(2)**int(d[1])
def fr(d):return Fraction(int(d[0]))*Fraction(2)**int(d[1])
def box(d):return dy(d['lower']).union(dy(d['upper']))
def m(z):return ((-4*z/7)*arb(2).log()).exp()*(4+z).gamma()*(4+8*z/7).rgamma()*((arb(fmpq(1,2))+z/7).gamma()*arb(fmpq(1,2)).rgamma())**8
qs=[]
for g,(sn,sd,rn,rd) in enumerate([(1,16,3,128),(3,16,3,64),(1,4,3,32)]):
 lo=arb(0);hi=arb(0)
 for item in P:
  if item['group']==g:
   L=item['integral']['lower'];U=item['integral']['upper'];assert 0<fr(L)<=fr(U);lo+=dy(L);hi+=dy(U)
 s=arb(fmpq(sn,sd));a=s/2;r=arb(fmpq(rn,rd));h=arb(fmpq(16,N[g]));norm=arb.pi()*s.gamma();assert 0<h<r<a
 E=32*(a-r).gamma()*(a+r).gamma()*m(a-r)*m(a+r)*(h/r)**129/(1-h/r);tail=m(a)**2*(a+16).gamma()**2/(norm*31*arb(32)**31)
 lower=((lo-E)/norm).lower();upper=((hi+E)/norm+tail).upper();cert=C['new_q'][g]['q'];assert fr(cert['lower'])<=fr(lower.man_exp()) and fr(cert['upper'])>=fr(upper.man_exp());qs.append(box(cert))
A=box(old['q']);R1=qs[0]*qs[1]/A**2;R2=qs[2]/A**2;mu=R1-1;nu=R2/2-2*R1+arb(fmpq(3,2));assert mu.lower()>0 and nu.lower()>mu.upper()**2
mul=mu.lower();nuu=nu.upper();x=1+nuu/mul;t2=2*x*x-1;t4=2*t2*t2-1;t8=2*t4*t4-1;gg=(-t8.log()/4).exp();delta=mul**2/nuu*(1-gg);dl=dy(C['delta_lower']['lower']);assert dl>0 and dl<1 and delta.lower()>=dl
M=dy(old['M_pair']['upper'])*(1-dl)**2;I=(M.log()*arb(fmpq(7,6))).exp()/(105*(arb(2).log()*arb(fmpq(25,6))).exp()*(arb.pi().log()*arb(fmpq(13,12))).exp());assert I.overlaps(box(C['I_new_upper']))
assert fr(C['I_new_upper']['upper'])<Fraction(C['rational_I_upper_numerator'],C['rational_I_upper_denominator'])<fr(old['I_pair']['lower'])
print('PASS: hashes, 1792 unique panel terminals, independent analytic tails and scalar interval assembly, reused R004 hash, Radau positivity and strict I improvement; no panel series recomputed.')
