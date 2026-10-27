#!/usr/bin/env python3
"""Read-only exact ledger and independent scalar assembly; no mesh/MC rerun."""
from pathlib import Path
from fractions import Fraction as Q
import json,hashlib,time
from flint import arb,acb,acb_series,fmpq,ctx
P=Path(__file__).resolve().parents[1]
def dy(x):return Q(int(x[0]))*Q(2)**int(x[1])
def A(x):x=Q(x);return arb(fmpq(x.numerator,x.denominator))
def ep(x):return [str(x.lower().fmpq()),str(x.upper().fmpq())]
def main():
 start=time.perf_counter();bindings={}
 for rel in ['upper/round004/pair_mellin','upper/round004/paired_cutoffs','audit/round004/pair_mellin']:
  b=P/rel;data=(b/'manifest.json').read_bytes();bindings[rel]=hashlib.sha256(data).hexdigest()
  for n,v in json.loads(data)['files'].items():
   h=v if isinstance(v,str)else v['sha256'];assert hashlib.sha256((b/n).read_bytes()).hexdigest()==h,n
 b=P/'upper/round004/strategy'
 for line in (b/'manifest.sha256').read_text().splitlines():
  h,n=line.split();assert hashlib.sha256((b/Path(n).name).read_bytes()).hexdigest()==h
 panels=[json.loads(x)for x in (P/'upper/round004/pair_mellin/panels.jsonl').read_text().splitlines()]
 assert len(panels)==512 and {x['panel']for x in panels}==set(range(512))
 L=Q(0);U=Q(0)
 for row in panels:
  assert row['status']=='completed'
  l,u=dy(row['integral']['lower']),dy(row['integral']['upper']);assert 0<l<=u;L+=l;U+=u
 ctx.prec=768
 def m(z):return ((-4*z/7)*A(2).log()+(4+z).lgamma()-(4+8*z/7).lgamma()+8*((A(Q(1,2))+z/7).lgamma()-A(Q(1,2)).lgamma())).exp()
 a=A(Q(1,16));r=A(Q(3,64));s=A(Q(1,8));norm=arb.pi()*s.gamma()
 star=(a-r).gamma()*(a+r).gamma()*m(a-r)*m(a+r)
 err=32*star*A(Q(2,3))**129*3
 tail=m(a)**2*(a+16).gamma()**2/(norm*31*A(32)**31)
 qlo=((A(L)-err)/norm).lower();qhi=((A(U)+err)/norm+tail).upper()
 Mlo=qlo**4/2;Mhi=qhi**4/2
 coeff=1/(105*A(2)**A(Q(25,6))*arb.pi()**A(Q(13,12)))
 ilo=Mlo**A(Q(7,6))*coeff;ihi=Mhi**A(Q(7,6))*coeff
 old=A(Q(8179950070308,10**18));assert ihi<old
 cert=json.loads((P/'upper/round004/pair_mellin/certificate.json').read_text())
 assert dy(cert['q']['lower'])<=Q(str(qlo.fmpq()))<Q(str(qhi.fmpq()))<=dy(cert['q']['upper'])
 audit=json.loads((P/'audit/round004/pair_mellin/results.json').read_text())
 oldIlo=Q(audit['scalar_enclosures']['R002_I_H'][0]);assert Q(ep(ihi)[1])<Q(8179950070308,10**18)<oldIlo
 # One newly registered t=0/order4 interface diagnostic, not a science panel.
 ctx.cap=4;z=acb_series([acb(a),acb(0,1)])
 rec=(z+1).gamma()-z*z.gamma();assert all(rec[k].contains(0)for k in range(4))
 f0=(a.gamma()*m(a))**2;assert f0>0
 F=z.gamma()*(4+z).gamma()*(4+8*z/7).rgamma()*(A(Q(1,2))+z/7).gamma()**8*((-4*z/7)*acb(A(2).log())).exp()/A(Q(1,2)).gamma()**8
 assert (F[0]*F[0].conjugate()-f0).contains(0)
 g=Q(1,208000);Iup=Q(ep(ihi)[1]);cost6=Q(7)**7*Iup**6/g**7
 assert cost6<Q(128)**6
 print(json.dumps(dict(status='PASS',bits=768,panels_recomputed=0,new_MC=0,unique_saved_panels=512,manifest_bindings=bindings,q=[str(qlo.fmpq()),str(qhi.fmpq())],M_pair=[ep(Mlo)[0],ep(Mhi)[1]],I_pair=[ep(ilo)[0],ep(ihi)[1]],new_t0_order4_Gamma_recurrence='PASS',new_t0_exact_formula=ep(f0),complete_I_rational_upper='8179950070308/1000000000000000000',conditional_g_target=str(g),conditional_cost_sixth_power_exact=str(cost6),conditional_below128=True,global_g_proved=False,new_C7_upper=False,elapsed_seconds=time.perf_counter()-start),indent=2))
if __name__=='__main__':main()
