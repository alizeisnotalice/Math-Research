"""Three-stage exact source-mass/threshold-budget pressure.
Mass vectors are realizable by disjoint unit boxes with these densities.
No original FIRST/early-exit or actual traffic claim.
"""
from fractions import Fraction as F
from decimal import Decimal,localcontext
from pathlib import Path
import json,hashlib

def dec(x):return Decimal(x.numerator)/Decimal(x.denominator)
def run():
 rows=[]
 with localcontext() as ctx:
  ctx.prec=80
  for stage,(K,d) in enumerate(((4,2),(16,4),(64,8)),1):
   H=F(2**(d+1));T=F(2**d)
   for name,p in [('uniform',[F(1,K)]*K),('geometric',[F(2**(K-1-i),2**K-1) for i in range(K)]),('heavy',[F(1,2)]+[F(1,2*(K-1))]*(K-1))]:
    for weight_name,w in [('proportional',p),('uniform',[F(1,K)]*K),('reversed',list(reversed(p)))]:
     m=[H*x for x in p];ratios=[x/y for x,y in zip(m,w)];shallow=[i for i,r in enumerate(ratios) if r<=T]
     paid=sum((m[i] for i in shallow),F(0));budget=sum((w[i] for i in shallow),F(0))
     assert paid<=T*budget<=T and paid<H
     lhs=sum((dec(x)*dec(r).ln() for x,r in zip(m,ratios)),Decimal(0))
     kl=sum((dec(x)*(dec(x)/dec(y)).ln() for x,y in zip(p,w)),Decimal(0))
     rhs=dec(H)*dec(H).ln()+dec(H)*kl
     assert abs(lhs-rhs)<Decimal('1e-68') and kl>=Decimal('-1e-70')
     rows.append(dict(stage=stage,K=K,depth=d,source=name,weights=weight_name,H=str(H),cutoff=str(T),shallow_mass=str(paid),shallow_threshold_weight=str(budget),uncapped_log_cost=str(lhs),KL=str(kl),identity_residual=str(abs(lhs-rhs))))
  # Boundary control: one source piece has mass exactly its allocated threshold;
  # its strict hard-heavy event is impossible. Global uncapped KL is not a
  # formula for capped/empty-fragment charges.
  m=[F(1,2),F(3,2)];w=[F(1,2),F(1,2)];B=Decimal(2).ln()
  capped=sum((dec(x)*(1+min(B,dec(x/y).ln())) if x>y else Decimal(0) for x,y in zip(m,w)),Decimal(0))
  baseline=Decimal(2)*(1+B)
  assert capped<baseline
  boundary={'masses':list(map(str,m)),'weights':list(map(str,w)),'capped_cost':str(capped),'unsplit_cost':str(baseline),'claim':'Capped/empty branches may improve cost; uncapped KL obstruction cannot be asserted globally.'}
 return rows,boundary
if __name__=='__main__':
 p=Path(__file__).with_name('fragmentation_pressure_results.json')
 if p.exists():raise SystemExit('Refuse overwrite')
 rows,boundary=run();r={'status':'PASS_MASS_BUDGET_ONLY','records':len(rows),'rounds':3,'boundary':boundary,'rows':rows,'seed':None,'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'scope':'Exact rational shallow-mass inequality and high precision log identity; not actual geom traffic.'};p.write_text(json.dumps(r,indent=2));print(json.dumps({k:v for k,v in r.items() if k!='rows'}))
