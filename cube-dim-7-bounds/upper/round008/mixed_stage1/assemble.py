#!/usr/bin/env python3
from flint import arb,fmpq,ctx
from fractions import Fraction
from pathlib import Path
import json,hashlib
ctx.prec=512;D=Path(__file__).resolve().parent
P=[json.loads(s) for s in (D/'panels.jsonl').read_text().splitlines()];assert len(P)==256 and {(x['group'],x['panel']) for x in P}=={(j,k) for j in range(4) for k in range(64)} and all(x['status']=='completed' for x in P)
b=(D/'r004-certificate.json').read_bytes();assert hashlib.sha256(b).hexdigest()==json.loads((D/'r004-manifest.json').read_text())['files']['certificate.json'];old=json.loads(b)
def dy(d):return arb(int(d[0]))*arb(2)**int(d[1])
def ep(a):
 l=a.lower().man_exp();u=a.upper().man_exp();return {'lower':[str(l[0]),int(l[1])],'upper':[str(u[0]),int(u[1])],'display':a.str(45)}
def F(z,j):
 e=arb(fmpq(1,7));ln2=arb(2).log();return ((-4*z/7-2*j*e)*ln2).exp()*(4+z-j*e).gamma()*(4+8*z/7-j*e).rgamma()*(arb(fmpq(1,2))+z/7-e).gamma()**j*(arb(fmpq(1,2))+z/7).gamma()**(8-j)/arb(fmpq(1,2)).gamma()**8
a=arb(fmpq(9,16));rho=arb(fmpq(27,64));h=arb(fmpq(1,4));b=arb(fmpq(9,8));norm=arb.pi()*b.gamma();assert 0<h<rho<a and (a-rho)>7*arb(fmpq(1,7))-arb(fmpq(7,2))
I=[];details=[]
for j in range(4):
 lo=arb(0);hi=arb(0)
 for x in P:
  if x['group']==j:lo+=dy(x['integral']['lower']);hi+=dy(x['integral']['upper'])
 c1=(a-rho).gamma()*F(a-rho,j)*(a+rho).gamma()*F(a+rho,7-j)
 c2=(a+rho).gamma()*F(a+rho,j)*(a-rho).gamma()*F(a-rho,7-j)
 u1=c1.upper();u2=c2.upper();cstar=u1 if u1>u2 else u2
 err=32*cstar*(h/rho)**129/(1-h/rho);tail=F(a,j)*F(a,7-j)*(a+16).gamma()**2/(norm*31*arb(32)**31)
 low=((lo-err)/norm-tail).lower();high=((hi+err)/norm+tail).upper();clipped=bool(low<0)
 if clipped:low=arb(0)
 assert high>0 and high>=low;I.append((low,high));details.append(dict(j=j,polynomial_lower=ep(lo),polynomial_upper=ep(hi),cross1=ep(c1),cross2=ep(c2),cstar_upper=ep(cstar),raw_error=ep(err),normalised_tail_abs=ep(tail),expectation_lower=ep(low),expectation_upper=ep(high),lower_clipped_to_zero=clipped))
co=[1,7,21,35];lo=2*sum((co[j]*I[j][0]**7 for j in range(4)),arb(0));hi=2*sum((co[j]*I[j][1]**7 for j in range(4)),arb(0));qlo=dy(old['q']['lower']);qhi=dy(old['q']['upper']);assert qlo>0 and qhi>0 and hi>0
LHlo=arb(0) if not lo>0 else (8*(lo.log()/7).exp()/qhi).lower();LHhi=(8*(hi.log()/7).exp()/qlo).upper();mul=(LHlo-1).lower();optim=(dy(old['M_pair']['lower'])/LHhi.sqrt()).lower();target=arb(fmpq(31,500));pass1=bool(mul>0);pass2=bool(optim<target);status='both_gates_pass' if pass1 and pass2 else 'mu_not_positive' if not pass1 else 'surrogate_cannot_reach_target' if optim>=target else 'target_comparison_unresolved'
out=dict(status=status,mixed_expectations=details,reused_q=old['q'],LH_lower=ep(LHlo),LH_upper=ep(LHhi),mu_lower=ep(mul),mu_positive=pass1,optimistic_surrogate_M_lower=ep(optim),target='31/500',target_not_ruled_out_by_lower_bound=pass2,stage_II_authorized=False,scope='surrogate capability only; no actual M lower or strict I, g or Cn bound')
(D/'certificate.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
