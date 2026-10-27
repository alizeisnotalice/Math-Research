#!/usr/bin/env python3
from flint import arb,fmpq,ctx
from fractions import Fraction
from pathlib import Path
import json,hashlib
ctx.prec=512;D=Path(__file__).resolve().parent;N=[2048,1024,512,512]
P=[json.loads(s) for s in (D/'panels.jsonl').read_text().splitlines()];assert len(P)==sum(N) and {(x['group'],x['panel']) for x in P}=={(g,j) for g,n in enumerate(N) for j in range(n)} and all(x['status']=='completed' for x in P)
def old(k):
 b=(D/f'r00{k}-certificate.json').read_bytes();mf=json.loads((D/f'r00{k}-manifest.json').read_text());assert hashlib.sha256(b).hexdigest()==mf['files']['certificate.json'];return json.loads(b)
old4=old(4);old5=old(5)
def dy(x):return arb(int(x[0]))*arb(2)**int(x[1])
def box(d):return dy(d['lower']).union(dy(d['upper']))
def ep(a):
 l=a.lower().man_exp();u=a.upper().man_exp();return {'lower':[str(l[0]),int(l[1])],'upper':[str(u[0]),int(u[1])],'display':a.str(45)}
def m(z):return ((-4*z/7)*arb(2).log()).exp()*(4+z).gamma()/(4+8*z/7).gamma()*((arb(fmpq(1,2))+z/7).gamma()/arb(fmpq(1,2)).gamma())**8
pi=arb.pi();Q=[];details=[]
for g,sn in enumerate([1,3,5,7]):
 lo=arb(0);hi=arb(0)
 for row in P:
  if row['group']==g:lo+=dy(row['integral']['lower']);hi+=dy(row['integral']['upper'])
 s=arb(fmpq(sn,32));a=s/2;r=3*a/4;h=arb(fmpq(16,N[g]));norm=pi*s.gamma(); assert 0<h<r<a
 star=(a-r).gamma()*(a+r).gamma()*m(a-r)*m(a+r);err=32*star*(h/r)**129/(1-h/r);tail=m(a)**2*(a+16).gamma()**2/(norm*31*arb(32)**31)
 qlo=((lo-err)/norm).lower();qhi=((hi+err)/norm+tail).upper();assert qlo>0 and qhi>qlo;q=qlo.union(qhi);Q.append(q);details.append(dict(s=f'{sn}/32',q=ep(q),Mstar=ep(star),eta=ep(h/r),quadrature_abs_error=ep(err),q_tail_bound=ep(tail),polynomial_lower=ep(lo),polynomial_upper=ep(hi)))
A=box(old4['q']);prev={z['s']:box(z['q']) for z in old5['new_q']};R1=Q[1]*Q[2]/A**2;R2=prev['1/16']*prev['3/16']/A**2;R3=Q[0]*Q[3]/A**2;R4=prev['1/4']/A**2
m0=arb(1);m1=R1-1;m2=R2/2-2*R1+arb(fmpq(3,2));m3=R3/4-3*R2/2+15*R1/4-arb(fmpq(5,2));m4=R4/8-R3+7*R2/2-7*R1+arb(fmpq(35,8))
D2=m1*m3-m2*m2;assert D2>0
alpha=(m2*m3-m1*m4)/D2;beta=(m2*m4-m3*m3)/D2;disc=alpha*alpha-4*beta;assert alpha<0 and beta>0 and disc>0
root2=(-alpha+disc.sqrt())/2;root1=beta/root2;assert root1>0 and root2>0 and root1.upper()<root2.lower()
w1=(m1*root2-m2)/(root1*(root2-root1));w2=(m2-m1*root1)/(root2*(root2-root1));w0=1-w1-w2;assert w0>0 and w1>0 and w2>0 and (w0+w1+w2).contains(1)
residual=[]
for k,mom in enumerate([m0,m1,m2,m3,m4]):
 v=(w0 if k==0 else 0)+w1*root1**k+w2*root2**k-mom;assert v.contains(0);residual.append(ep(v))
def correction(d,n):
 t=1+d;v=d
 for _ in range(n):v=2*v*(t+1);t=2*t*t-1
 R=t.sqrt().sqrt();return v/(R*(R**3+R**2+R+1))
delta=w1.lower()*correction(root1.lower(),4)+w2.lower()*correction(root2.lower(),4);dl=delta.lower();assert 0<dl<1
B=1-dl;M=dy(old4['q']['upper'])**4*B**2/2;F=1/(105*(arb(2).log()*arb(fmpq(25,6))).exp()*(pi.log()*arb(fmpq(13,12))).exp());I=(M.log()*arb(fmpq(7,6))).exp()*F
# Lower endpoint of the actual R005 analytical budget expression, not an actual I lower bound.
mu5=R2-1;nu5=R4/2-2*R2+arb(fmpq(3,2));assert mu5>0 and nu5>0 and nu5.lower()>mu5.upper()**2
delta5=mu5.upper()**2/nu5.lower()*correction(nu5.lower()/mu5.upper(),3);d5u=delta5.upper();assert 0<d5u<1
M5low=dy(old4['q']['lower'])**4*(1-d5u)**2/2;I5low=(M5low.log()*arb(fmpq(7,6))).exp()*F;assert I.upper()<I5low.lower()
up=I.upper().man_exp();f=Fraction(int(up[0]))*Fraction(2)**int(up[1]);grid=10**18;n=(f.numerator*grid+f.denominator-1)//f.denominator;assert Fraction(n,grid)>f
out=dict(new_q=details,reused_q_1_8=old4['q'],reused_q_R005=old5['new_q'],R=[ep(v) for v in [R1,R2,R3,R4]],moments=[ep(v) for v in [m1,m2,m3,m4]],D=ep(D2),alpha=ep(alpha),beta=ep(beta),discriminant=ep(disc),roots=[ep(root1),ep(root2)],weights=[ep(w0),ep(w1),ep(w2)],moment_residuals=residual,delta_lower=ep(dl),B_upper=ep(B),M_upper=ep(M),I_upper=ep(I),R005_exact_budget_expression_lower=ep(I5low),strict_improvement=True,rational_I_upper_numerator=n,rational_I_upper_denominator=grid)
(D/'certificate.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
