#!/usr/bin/env python3
"""Independent scalar reassembly and integrity check, without repeating panels."""
from flint import arb,fmpq,ctx
from fractions import Fraction
from pathlib import Path
import json,hashlib
ctx.prec=512;p=Path(__file__).resolve().parent;sha=lambda n:hashlib.sha256((p/n).read_bytes()).hexdigest();manifest=json.loads((p/'manifest.json').read_text())
for n,h in manifest['files'].items():assert sha(n)==h,n
P=[json.loads(x) for x in (p/'panels.jsonl').read_text().splitlines()];N=[2048,1024,512,512];assert len(P)==4096 and {(x['group'],x['panel']) for x in P}=={(g,i) for g,n in enumerate(N) for i in range(n)} and all(x['status']=='completed' for x in P)
for k in (4,5):assert sha(f'r00{k}-certificate.json')==json.loads((p/f'r00{k}-manifest.json').read_text())['files']['certificate.json']
O4=json.loads((p/'r004-certificate.json').read_text());O5=json.loads((p/'r005-certificate.json').read_text());C=json.loads((p/'certificate.json').read_text())
def dy(d):return arb(int(d[0]))*arb(2)**int(d[1])
def fr(d):return Fraction(int(d[0]))*Fraction(2)**int(d[1])
def box(d):return dy(d['lower']).union(dy(d['upper']))
def m(z):return ((-4*z/7)*arb(2).log()).exp()*(4+z).gamma()*(4+8*z/7).rgamma()*((arb(fmpq(1,2))+z/7).gamma()*arb(fmpq(1,2)).rgamma())**8
q=[]
for g,num in enumerate([1,3,5,7]):
 lo=arb(0);hi=arb(0)
 for x in P:
  if x['group']==g:
   l=x['integral']['lower'];u=x['integral']['upper'];assert 0<fr(l)<=fr(u);lo+=dy(l);hi+=dy(u)
 s=arb(fmpq(num,32));a=s/2;r=3*a/4;h=arb(fmpq(16,N[g]));norm=arb.pi()*s.gamma();assert 0<h<r<a
 E=32*(a-r).gamma()*(a+r).gamma()*m(a-r)*m(a+r)*(h/r)**129/(1-h/r);tail=m(a)**2*(a+16).gamma()**2/(norm*31*arb(32)**31)
 l=((lo-E)/norm).lower();u=((hi+E)/norm+tail).upper();c=C['new_q'][g]['q'];assert fr(c['lower'])<=fr(l.man_exp()) and fr(c['upper'])>=fr(u.man_exp());q.append(box(c))
A=box(O4['q']);old={v['s']:box(v['q']) for v in O5['new_q']};R1=q[1]*q[2]/A**2;R2=old['1/16']*old['3/16']/A**2;R3=q[0]*q[3]/A**2;R4=old['1/4']/A**2
M1=R1-1;M2=R2/2-2*R1+arb(fmpq(3,2));M3=R3/4-3*R2/2+15*R1/4-arb(fmpq(5,2));M4=R4/8-R3+7*R2/2-7*R1+arb(fmpq(35,8));D=M1*M3-M2*M2;assert D>0
alpha=(M2*M3-M1*M4)/D;beta=(M2*M4-M3*M3)/D;disc=alpha*alpha-4*beta;assert alpha<0 and beta>0 and disc>0
b=(-alpha+disc.sqrt())/2;a=beta/b;assert a>0 and a.upper()<b.lower();w1=(M1*b-M2)/(a*(b-a));w2=(M2-M1*a)/(b*(b-a));w0=1-w1-w2;assert w0>0 and w1>0 and w2>0 and (w0+w1+w2).contains(1)
for k,M in enumerate([arb(1),M1,M2,M3,M4]):assert ((w0 if k==0 else 0)+w1*a**k+w2*b**k-M).contains(0)
def corr(d,n):
 x=1+d;v=d
 for _ in range(n):v=2*v*(x+1);x=2*x*x-1
 root=x.sqrt().sqrt();return v/(root*(root**3+root**2+root+1))
Delta=w1.lower()*corr(a.lower(),4)+w2.lower()*corr(b.lower(),4);dl=dy(C['delta_lower']['lower']);assert 0<dl<1 and Delta.lower()>=dl
F=1/(105*(arb(2).log()*arb(fmpq(25,6))).exp()*(arb.pi().log()*arb(fmpq(13,12))).exp());MM=dy(O4['q']['upper'])**4*(1-dl)**2/2;II=(MM.log()*arb(fmpq(7,6))).exp()*F;assert II.overlaps(box(C['I_upper']))
mu=R2-1;nu=R4/2-2*R2+arb(fmpq(3,2));assert mu>0 and nu.lower()>mu.upper()**2;delta5=(mu.upper()**2/nu.lower()*corr(nu.lower()/mu.upper(),3)).upper();assert 0<delta5<1
MM5=dy(O4['q']['lower'])**4*(1-delta5)**2/2;II5=(MM5.log()*arb(fmpq(7,6))).exp()*F
assert fr(C['I_upper']['upper'])<Fraction(C['rational_I_upper_numerator'],C['rational_I_upper_denominator'])<fr(II5.lower().man_exp())
print('PASS: source hashes, 4096 unique terminal panels, independent q tails and 4-moment root/weight gates, correction direction, and strict comparison with exact R005 budget lower. No panel recomputation or C7 claim.')
