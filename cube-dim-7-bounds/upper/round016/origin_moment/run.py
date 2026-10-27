#!/usr/bin/env python3
from flint import arb,fmpq,ctx
from fractions import Fraction
from pathlib import Path
import hashlib,json,math,time,datetime,resource,sys,flint
ctx.prec=512;D=Path(__file__).resolve().parent;now=lambda:datetime.datetime.now(datetime.timezone.utc).isoformat()
def ep(x):
 l=x.lower().man_exp();u=x.upper().man_exp();return {'lower':[str(l[0]),int(l[1])],'upper':[str(u[0]),int(u[1])],'display':x.str(45)}
def frac(a):return Fraction(int(a[0]))*Fraction(2)**int(a[1])
def A(x):return arb(fmpq(x.numerator,x.denominator)) if isinstance(x,Fraction) else arb(x)
def upgrid(x):
 upper=frac(x.upper().man_exp());N=10**18;n=(upper.numerator*N+upper.denominator-1)//upper.denominator;assert upper<Fraction(n,N);return Fraction(n,N)
def main():
 start=now();t=time.perf_counter();c=time.process_time()
 for f,m,k in [('r010-cells.jsonl','r010-manifest.json','cells.jsonl'),('r006-certificate.json','r006-manifest.json','certificate.json')]:assert hashlib.sha256((D/f).read_bytes()).hexdigest()==json.loads((D/m).read_text())['files'][k]
 cells=[json.loads(x) for x in (D/'r010-cells.jsonl').read_text().splitlines()];assert len(cells)==135;seen=set();S=Fraction()
 for cell in cells:
  for x in cell['simplices']:
   key=tuple(x['vertices']);assert key not in seen;seen.add(key);v=Fraction(abs(x['detW']),math.prod(x['q']));assert v==Fraction(x['alpha']);S+=v
 assert len(seen)==212 and 0<5040*S<=Fraction(1,512)
 sharp=json.loads(Path('/Users/zhengzhihao/Documents/Codex/2026-10-27/n-n-5-1-1-verify/work/r011-primary-euclidean-floor.json').read_text());assert sharp['max_squared_unitball_norm']=='1/2048' and sharp['psi_euclidean_lower_squared']==2048 and hashlib.sha256(Path(sharp['source']).read_bytes()).hexdigest()==sharp['source_sha256']
 pi=arb.pi();p0=5040*A(S)/pi**7;csharp=32*arb(2).sqrt();p0sharp=6/(pi**4*csharp**7);M=2*pi.sqrt()*(p0.log()/7).exp();via=(M.log()*arb(fmpq(7,6))).exp()/(105*(arb(2).log()*arb(fmpq(25,6))).exp()*(pi.log()*arb(fmpq(13,12))).exp());Bfan=(p0.log()/6).exp()/(840*pi.sqrt());Bsharp=(p0sharp.log()/6).exp()/(840*pi.sqrt())
 ratio=arb(fmpq(65,2)).gamma()/arb(32).gamma();b=arb(fmpq(65,2))*arb(fmpq(67,2))*arb(fmpq(69,2));eta=ratio/(b.log()/6).exp();assert (arb(fmpq(71,2)).gamma()/arb(32).gamma()).overlaps(ratio*b);BRB=eta*Bfan
 assert p0>0 and p0sharp-p0>0 and via.overlaps(Bfan) and eta<1 and eta>0 and Bfan-BRB>0
 rawup=upgrid(Bfan);rbup=upgrid(BRB);old=json.loads((D/'r006-certificate.json').read_text());oldlower=frac(old['I_upper']['lower']);assert rawup<=Fraction(4127002498830,10**18) and Bsharp<A(Fraction(4940777144766,10**18)) and rbup<rawup<oldlower
 req=7*(A(rbup/128).log()*arb(fmpq(6,7))).exp();g=Fraction(1,320000);assert req<A(g)
 rawcost6=Fraction(7**7)*rawup**6/g**7;rbcost6=Fraction(7**7)*rbup**6/g**7;assert rawcost6<128**6 and rbcost6<rawcost6;rawcost=(A(rawcost6).log()/6).exp();rbcost=(A(rbcost6).log()/6).exp()
 primary=json.loads(Path('/Users/zhengzhihao/Documents/Codex/2026-10-27/n-n-5-1-1-verify/work/r016-origin-moment-primary.json').read_text());assert primary['A_sum']==str(S) and hashlib.sha256(Path(primary['source']).read_bytes()).hexdigest()==primary['source_sha256'];primaryI=[Fraction(v) for v in primary['I_bound']];assert frac(Bfan.lower().man_exp())<=primaryI[1] and primaryI[0]<=frac(Bfan.upper().man_exp())
 balls={'p0_fan':ep(p0),'p0_sharp_upper':ep(p0sharp),'M_raw_bound':ep(M),'B_fan':ep(Bfan),'B_sharp':ep(Bsharp),'eta':ep(eta),'B_RB':ep(BRB),'g_required_B_RB':ep(req),'conditional_cost_B_RB':ep(rbcost),'conditional_cost_B_fan':ep(rawcost)};assert len(balls)==10
 out=dict(status='completed',A_sum=str(S),simplices=212,scalar_enclosures=len(balls),precision_bits=512,enclosures=balls,B_fan_rational_upper=str(rawup),B_RB_rational_upper=str(rbup),R006_deterministic_budget_expression_lower=old['I_upper']['lower'],p0_below_sharp=True,raw_M_interface_overlap=True,gamma_recurrence_overlap=True,eta_below_one=True,B_RB_below_B_fan=True,conditional_g='1/320000',conditional_cost6_raw_exact=str(rawcost6),conditional_cost6_RB_exact=str(rbcost6),conditional_costs_strictly_below_128=True,global_g_proved=False,new_C7_bound=False,primary_raw_scalar_I_overlap=True,raw_seventh_moment_identity='E[det(Sigma)^(-1/2)]=(4*pi)^(7/2)*p1(0); Z=det(Sigma)^(-1/14), Z^7=det(Sigma)^(-1/2)',radial_moment_identity='E Z^7=128 Gamma(71/2)/Gamma(32) E_P Z_P^7; not D^7')
 (D/'result.json').write_text(json.dumps(out,indent=2)+'\n');receipt=dict(status='completed',start_utc=start,finish_utc=now(),wall_seconds=time.perf_counter()-t,cpu_seconds=time.process_time()-c,peak_rss_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,python=sys.version,flint=flint.__version__,named_scalar_enclosures=10,MC=0,phase_evaluations=0,threads=1);(D/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt,indent=2));print(json.dumps({k:v['display'] for k,v in balls.items()},indent=2));print('B_fan_up',rawup,'B_RB_up',rbup)
if __name__=='__main__':main()
