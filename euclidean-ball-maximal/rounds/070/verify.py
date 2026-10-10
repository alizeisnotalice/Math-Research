#!/usr/bin/env python3
"""Full-line Fraction checks for separated continuous densities and small tails."""
from fractions import Fraction as F
from pathlib import Path
import runpy, json, hashlib, random, sys
HERE=Path(__file__).resolve().parent
PREV=HERE.parent/'069'/'verify.py'
old=runpy.run_path(str(PREV))
levels,mass,length,merge=(old[k] for k in ('levels','mass','length','merge'))

def canonical(data):
    ends=sorted({x for a,b,h in data for x in (a,b)})
    return [(a,b,sum((h for l,r,h in data if l<(a+b)/2<r),F(0)))
            for a,b in zip(ends,ends[1:])]

def capmass(data,t):
    return sum(((b-a)*min(h,t) for a,b,h in canonical(data)),F(0))

def main():
    batches=[]
    for batch,(m,N) in enumerate(((1,2),(2,3),(3,4))):
        rng=random.Random(700000+batch)
        base=[(F(2*j),F(2*j)+F(rng.randint(1,3),4),F(rng.randint(1,9))) for j in range(m)]
        W=mass(base);base=[(a,b,h/W) for a,b,h in base]
        replicated=extracted=0
        for alpha in (F(1,8),F(1,2),F(2)):
            spacing=F(100*m)+10/alpha
            blocks=[[(a/N+i*spacing,b/N+i*spacing,h) for a,b,h in base] for i in range(N)]
            full=[v for block in blocks for v in block]
            assert mass(full)==1
            individual=[levels(block,alpha) for block in blocks]
            got=levels(full,alpha)
            assert got==merge([v for lev in individual for v in lev])
            assert length(got)==length(levels(base,alpha))
            assert sum((length(lev) for lev in individual),F(0))==length(got)
            for L in (F(1,4),F(1),F(8)):
                assert capmass(full,L*alpha)==capmass(base,L*alpha)
            assert sum((b-a for a,b,h in full),F(0))==sum((b-a for a,b,h in base),F(0))
            replicated+=1
            for eps,delta in ((F(1,100),F(1,10)),(F(1,4),F(1,2))):
                tau=1-delta
                weights=[F(rng.randint(1,10**(batch+1))) for _ in blocks]
                qs=[(1-eps)*w/sum(weights) for w in weights]
                scaled=[[(a,b,N*q*h) for a,b,h in block] for block,q in zip(blocks,qs)]
                # Test a diffuse bridge and a high-density residual in the gap.
                if eps<F(1,10):
                    left=F(-1);right=blocks[-1][-1][1]+1
                    tail=[(left,right,eps/(right-left))]
                else:
                    left=spacing/2
                    tail=[(left,left+eps/(10*alpha),10*alpha)]
                for first,second in zip(scaled,scaled[1:]):
                    assert second[0][0]-first[-1][1] >= (1-eps)/(tau*alpha)
                data=canonical([v for block in scaled for v in block]+tail)
                assert mass(data)==1
                flev=levels(data,alpha); hlev=levels(tail,delta*alpha)
                ilev=[levels(block,tau*alpha) for block in scaled]
                low=levels([v for block in scaled for v in block],tau*alpha)
                assert low==merge([v for lev in ilev for v in lev])
                assert length(low)==sum((length(lev) for lev in ilev),F(0))
                assert length(merge(low+hlev+flev))==length(merge(low+hlev))
                R=alpha*length(flev)
                Ri=[tau*alpha*length(lev)/q for lev,q in zip(ilev,qs)]
                Rh=delta*alpha*length(hlev)/eps
                # Finite envelope only: NOT the true optimal constant C_1.
                C=max([F(1),R,Rh,*Ri]);eta=1-R/C
                D=sum((q*(1-v/C) for q,v in zip(qs,Ri)),F(0))
                bound=delta+tau*eta+eps*(tau/delta-1)
                assert D<=bound<=eta+2*delta  # eps=delta^2
                for kappa in (F(1,10),F(1,2)):
                    bad=sum((q for q,v in zip(qs,Ri) if v/C<1-kappa),F(0))
                    assert kappa*bad<=D
                for block,lev,q in zip(scaled,ilev,qs):
                    normalized=[(a,b,h/q) for a,b,h in block]
                    assert levels(normalized,tau*alpha/q)==lev
                extracted+=1
        # Dimension-dependent replication scaling is checked without roots:
        # choose N=k^n, so its spatial contraction is exactly 1/k.
        dims=((1,2,3),(8,16,32),(64,128,256))[batch]
        certificates=0
        counterexamples=[]
        for k in ((6,8,10),(12,16,20),(24,32,48))[batch]:
            eps=F(1,2**k)
            pair=[(c-eps,c+eps,F(1,4)/eps) for c in (F(-2,5),F(2,5))]
            whole=levels(pair,F(1))
            separate=merge([v for block in pair for v in levels([block],F(1))])
            target=[[F(-1,20),F(1,20)]]
            assert length(merge(whole+target))==length(whole)
            assert length(merge(separate+target))==length(separate)+F(1,10)
            missed=length(merge(whole+separate))-length(separate)
            assert missed>=F(1,10)
            counterexamples.append(dict(epsilon_power=k,missed_volume=str(missed),carrier_volume=str(4*eps)))
        for n in dims:
            for k in (2,3):
                count=k**n
                assert count*F(1,k)**n==1
                assert count*F(1,count)==1
                certificates+=1
        batches.append(dict(batch=batch,seed=700000+batch,source_intervals=m,copies=N,
                            replication_levels=replicated,tail_extractions=extracted,
                            dimensions=list(dims),jacobian_certificates=certificates,counterexamples=counterexamples))
    print(json.dumps(dict(round=70,status='passed',python=sys.version.split()[0],batches=batches,
        sha256={str(p.relative_to(HERE.parent.parent)):hashlib.sha256(p.read_bytes()).hexdigest() for p in (Path(__file__).resolve(),PREV)},
        scope='Exact full-line continuous-density levels; abstract finite envelope is not C_1; high-dimensional checks are scaling identities only.',
        main_theorem_proved=False),indent=2))
if __name__=='__main__':main()
