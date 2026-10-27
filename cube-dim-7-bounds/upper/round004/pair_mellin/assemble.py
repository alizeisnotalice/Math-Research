#!/usr/bin/env python3
from flint import arb,fmpq,ctx
from fractions import Fraction
import json
from pathlib import Path
ctx.prec=512
D=Path(__file__).resolve().parent
panels=[json.loads(s) for s in (D/'panels.jsonl').read_text().splitlines()]
assert len(panels)==512 and {p['panel'] for p in panels}==set(range(512)) and all(p['status']=='completed' for p in panels)
def dy(x):return arb(int(x[0]))*arb(2)**x[1]
def ep(a):
 l=a.lower().man_exp();u=a.upper().man_exp();return {'lower':[str(l[0]),int(l[1])],'upper':[str(u[0]),int(u[1])],'display':a.str(45)}
def m(z):return ((-4*z/7)*arb(2).log()).exp()*(4+z).gamma()/(4+8*z/7).gamma()*((arb(fmpq(1,2))+z/7).gamma()/arb(fmpq(1,2)).gamma())**8
lo=arb(0);hi=arb(0)
for p in panels:lo+=dy(p['integral']['lower']);hi+=dy(p['integral']['upper'])
a=arb(fmpq(1,16));r=arb(fmpq(3,64));s=arb(fmpq(1,8));pi=arb.pi();norm=pi*s.gamma();star=(a-r).gamma()*(a+r).gamma()*m(a-r)*m(a+r);err=32*star*arb(fmpq(2,3))**129*3
# tail is already in q units
qtail=m(a)**2*(a+16).gamma()**2/(norm*31*arb(32)**31)
qlo=((lo-err)/norm).lower();qhi=((hi+err)/norm+qtail).upper();assert qlo>0 and qhi>qlo
q=qlo.union(qhi);M=q**4/2;I=(M.log()*arb(fmpq(7,6))).exp()/(105*(arb(2).log()*arb(fmpq(25,6))).exp()*(pi.log()*arb(fmpq(13,12))).exp())
old=(1/arb(8).sqrt())*(arb(2)**0 if False else ((-arb(fmpq(1,28)))*arb(2).log()).exp()*(arb(fmpq(65,16))).gamma()/(arb(fmpq(57,14))).gamma()*((arb(fmpq(57,112)).gamma()/arb(fmpq(1,2)).gamma())**8))**8
oldI=(old.log()*arb(fmpq(7,6))).exp()/(105*(arb(2).log()*arb(fmpq(25,6))).exp()*(pi.log()*arb(fmpq(13,12))).exp());assert old-I>0
u=I.upper().man_exp();f=Fraction(int(u[0]))*(Fraction(2)**int(u[1]));grid=10**18;rounded=(f.numerator*grid+f.denominator-1)//f.denominator;assert Fraction(rounded,grid)>f
out={'polynomial_lower':ep(lo),'polynomial_upper':ep(hi),'M_star':ep(star),'quadrature_absolute_error':ep(err),'q_tail_upper_bound':ep(qtail),'q':ep(q),'M_pair':ep(M),'I_pair':ep(I),'M_H_R002':ep(old),'I_H_R002':ep(oldI),'strict_I_improvement':bool(oldI-I>0),'rational_I_upper_numerator':rounded,'rational_I_upper_denominator':grid,'outward_rational_proved':True}
(D/'certificate.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
