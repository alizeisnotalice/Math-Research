#!/usr/bin/env python3
from flint import arb,fmpq,ctx
from fractions import Fraction
from pathlib import Path
import json,hashlib
ctx.prec=512
D=Path(__file__).resolve().parent
P=[json.loads(s) for s in (D/'panels.jsonl').read_text().splitlines()]; ns=[1024,512,256];assert len(P)==sum(ns) and {(x['group'],x['panel']) for x in P}=={(g,j) for g,n in enumerate(ns) for j in range(n)} and all(x['status']=='completed' for x in P)
oldraw=(D/'r004-certificate.json').read_bytes();oldman=json.loads((D/'r004-manifest.json').read_text());assert hashlib.sha256(oldraw).hexdigest()==oldman['files']['certificate.json'];old=json.loads(oldraw)
def dy(x):return arb(int(x[0]))*arb(2)**x[1]
def box(d):return dy(d['lower']).union(dy(d['upper']))
def ep(a):
 l=a.lower().man_exp();u=a.upper().man_exp();return {'lower':[str(l[0]),int(l[1])],'upper':[str(u[0]),int(u[1])],'display':a.str(45)}
def m(z):return ((-4*z/7)*arb(2).log()).exp()*(4+z).gamma()/(4+8*z/7).gamma()*((arb(fmpq(1,2))+z/7).gamma()/arb(fmpq(1,2)).gamma())**8
pi=arb.pi(); qs=[];details=[]
for g,(sn,sd,rn,rd) in enumerate([(1,16,3,128),(3,16,3,64),(1,4,3,32)]):
 lo=arb(0);hi=arb(0)
 for p in P:
  if p['group']==g:lo+=dy(p['integral']['lower']);hi+=dy(p['integral']['upper'])
 s=arb(fmpq(sn,sd));a=s/2;r=arb(fmpq(rn,rd));h=arb(fmpq(16,ns[g]));norm=pi*s.gamma(); assert 0<h<r<a
 star=(a-r).gamma()*(a+r).gamma()*m(a-r)*m(a+r);err=32*star*(h/r)**129/(1-h/r);tail=m(a)**2*(a+16).gamma()**2/(norm*31*arb(32)**31)
 qlo=((lo-err)/norm).lower();qhi=((hi+err)/norm+tail).upper();assert qlo>0 and qhi>qlo;q=qlo.union(qhi);qs.append(q)
 details.append(dict(s=f'{sn}/{sd}',q=ep(q),Mstar=ep(star),quadrature_abs_error=ep(err),q_tail_bound=ep(tail),polynomial_lower=ep(lo),polynomial_upper=ep(hi)))
q0=box(old['q']);R1=qs[0]*qs[1]/q0**2;R2=qs[2]/q0**2;mu=R1-1;nu=R2/2-2*R1+arb(fmpq(3,2));assert mu.lower()>0 and nu.lower()>mu.upper()**2
mul=mu.lower();nuu=nu.upper();x=1+nuu/mul;t0=arb(1);t1=x
for _ in range(2,9):t0,t1=t1,2*x*t1-t0
assert t1>1;g=(-(t1.log())/4).exp();delta=(mul**2/nuu)*(1-g);dl=delta.lower();assert dl>0 and dl<1
Mold=dy(old['M_pair']['upper']);Mnew=Mold*(1-dl)**2;I=(Mnew.log()*arb(fmpq(7,6))).exp()/(105*(arb(2).log()*arb(fmpq(25,6))).exp()*(pi.log()*arb(fmpq(13,12))).exp());assert I<dy(old['I_pair']['lower'])
u=I.upper().man_exp();f=Fraction(int(u[0]))*Fraction(2)**int(u[1]);N=10**18;up=(f.numerator*N+f.denominator-1)//f.denominator;assert Fraction(up,N)>f
out=dict(new_q=details,reused_q_1_8=old['q'],R1=ep(R1),R2=ep(R2),mu=ep(mu),nu=ep(nu),nu_gt_mu_squared=True,delta_lower=ep(dl),M_new_upper=ep(Mnew),I_new_upper=ep(I),old_I=old['I_pair'],strict_improvement=True,rational_I_upper_numerator=up,rational_I_upper_denominator=N)
(D/'certificate.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
