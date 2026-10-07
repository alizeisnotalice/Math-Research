"""Finite exact checks for the mass-truncated hard-kernel lemma only.
Not the original soft FIRST/early-exit traffic and not a dimension fit.
"""
from fractions import Fraction as F
from decimal import Decimal, localcontext
from pathlib import Path
import json,hashlib

def run():
    tau=F(49,32768);lam=F(1);alpha=2*tau*lam;rows=[]
    for stage,(n,N,na) in enumerate(((4,32,4),(16,128,8),(64,512,16)),1):
        lo=sum((F(1,N+i+1) for i in range(N)),F(0))
        hi=sum((F(1,N+i) for i in range(N)),F(0))
        with localcontext() as ctx:
            ctx.prec=60;ln2=Decimal(2).ln()
            assert Decimal(lo.numerator)/Decimal(lo.denominator)<ln2<Decimal(hi.numerator)/Decimal(hi.denominator)
        for depth in (0,1,n//2,n,n+1):
            mass=alpha*2**depth; weights=[mass*F(i+1,na*(na+1)//2) for i in range(na)]
            source=[tuple(F((i+1)*(j%3+1),2*na) for j in range(n)) for i in range(na)]
            cap_volume=min(F(2)**n,mass/alpha); checks=0;events=0;ties=0
            for probe in (tuple(F(0) for _ in range(n)),tuple(F(1,2) for _ in range(n)),tuple(F(1) for _ in range(n))):
                dist=[max(abs(x-y) for x,y in zip(probe,s)) for s in source]
                envelope=sum((w*(F(1) if d<=F(1,2) else (F(1)/(2*d)**n if (2*d)**n<=cap_volume else F(0))) for w,d in zip(weights,dist)),F(0))
                for side in (F(1),F(5,4),F(3,2),F(7,4),F(2)):
                    capture=sum((w for w,d in zip(weights,dist) if 2*d<=side),F(0));response=capture/side**n
                    ties+=int(response==alpha)
                    if response>alpha:
                        assert mass>alpha and side**n<mass/alpha
                        assert response<=envelope
                        events+=1
                    checks+=1
            q=min(depth,n)
            rows.append(dict(stage=stage,n=n,mass_depth=depth,atoms=na,rows=checks,strict_events=events,equality_rows=ties,cap_volume=str(cap_volume),norm_lower=str(1+q*lo),norm_upper=str(1+q*hi),norm_formula=f'1+{q}*log(2)',empty_when_mass_le_threshold=(depth==0)))
    return rows
if __name__=='__main__':
    p=Path(__file__).with_name('heavy_cutoff_pressure_results.json')
    if p.exists():raise SystemExit('Refuse overwrite')
    rows=run();r={'status':'PASS_FINITE_COMPONENT','rounds':3,'records':len(rows),'receiver_radius_rows':sum(x['rows'] for x in rows),'strict_events':sum(x['strict_events'] for x in rows),'seed':None,'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'rows':rows,'scope':'Positive atomic component kernel and rational shell enclosures. No original FIRST, CP/GP, LCA, early-jump, or spatial-volume certification.'};p.write_text(json.dumps(r,indent=2));print(json.dumps({k:v for k,v in r.items() if k!='rows'}))
