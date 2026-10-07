"""Exact finite positive-measure guards; NOT original kernel/FIRST samples."""
from fractions import Fraction as F
from random import Random
from pathlib import Path
import json
rng=Random(20261007019)
rounds=[]
for r,N in enumerate((8,32,128),1):
    cases=[]
    for trial in range(12):
        tau=F(1,4+trial)
        mu=[F(rng.randrange(1,33),32) for _ in range(N)]
        p=[F(rng.randrange(1,65),64) for _ in range(N)]
        # Correlated histories and hard gates, including Z=tau exactly.
        Z=[tau*F(rng.randrange(0,9),4) for _ in range(N)]
        Z[0]=tau
        k=[p[i]*Z[i] for i in range(N)]
        q=sum(mu[i]*p[i] for i in range(N))
        u=q*F(16,3)
        hard=[F(1,6),F(1,3),F(1,2)]
        lhs=F(0)
        raw=F(0)
        for i in range(N):
            split=[F(rng.randrange(0,5),12) for _ in range(3)]
            assert sum(split)<=1
            if Z[i]>=tau:
                raw+=mu[i]*p[i]
                for j in range(3):
                    for z in range(3):
                        gate=F((i+2*j+3*z+trial)%8,7)
                        lhs+=u*mu[i]*p[i]/q*split[j]*hard[z]*gate
        middle=F(16,3)*raw
        rhs=F(16,3)*sum(mu[i]*k[i] for i in range(N))/tau
        assert lhs<=middle<=rhs
        # Direct equality threshold belongs to paid side.
        assert Z[0]>=tau
        cases.append({'trial':trial,'tau':str(tau),'traffic':str(lhs),'posterior_bound':str(middle),'kernel_bound':str(rhs)})
    # Low row-average need not mean every original source has low ratio.
    rho=[F(1,100),F(99,100)]; zz=[F(2),F(1,1000)]
    assert sum(a*b for a,b in zip(rho,zz))<F(1,4) and zz[0]>=F(1,4)
    # Total subposterior mass alone does not imply original-source marginal domination.
    bad=[F(1),F(0)]
    assert sum(bad)==sum(rho) and bad[0]>rho[0]
    rounds.append({'round':r,'sources':N,'cases':cases,'low_average_high_source_guard':True,'mass_only_insufficiency_guard':True})
out={'scope':'Exact algebra only. Arbitrary correlated finite gates; not actual FIRST nor asymptotics.','seed':20261007019,'rounds':rounds,'pass':True}
Path(__file__).with_name('low_future_source_payment_guard_20261007_results.json').write_text(json.dumps(out,ensure_ascii=False,indent=2))
print(json.dumps({'pass':True,'rounds':3,'cases':36,'sources':[8,32,128]}))
