"""Exact small different-channel diagonal/Gram and common-space L1 checks."""
from pathlib import Path
from fractions import Fraction as F
import itertools,math,json
HERE=Path(__file__).resolve().parent
def pb(ps):
    out=[F(1)]
    for p in ps:
        nxt=[F(0)]*(len(out)+1)
        for k,v in enumerate(out):nxt[k]+=v*(1-p);nxt[k+1]+=v*p
        out=nxt
    return out
def integrate_gram(n,weights,ps,prior):
    gram=[[F(0)for s in range(n+1)]for t in range(n+1)]
    for obs in itertools.product(range(4),repeat=n):
        w=math.prod(weights[o]for o in obs)
        L=[v/pi for v,pi in zip(pb([ps[o]for o in obs]),prior)]
        for s in range(n+1):
            for t in range(n+1):gram[s][t]+=w*L[s]*L[t]
    return gram
out=[]
for eta in [F(1,2),F(2,5)]:
    flat=[F(1,4)]*4;peaked=[F(1,8),F(3,8),F(3,8),F(1,8)]
    channels=[(flat,[eta-F(1,20)]*2+[eta+F(1,20)]*2),
              (peaked,[eta-F(1,10),eta,eta,eta+F(1,10)]),
              (flat,[eta-F(1,10)]*2+[eta+F(1,10)]*2),
              (peaked,[eta-F(1,8),eta,eta,eta+F(1,8)])]
    fs=[]
    for w,p in channels:
        assert sum(a*b for a,b in zip(w,p))==eta
        fs.append(sum(a*(b-eta)**2 for a,b in zip(w,p))/(eta*(1-eta)))
    assert fs[0]==fs[1] and max(fs)==fs[2]
    for n in [3,5]:
        prior=pb([eta]*n);R=4;cumulative=list(itertools.accumulate(prior))
        upper=[next(s for s,c in enumerate(cumulative)if R*c>=b)for b in range(1,R)]+[n]
        lower=[0]+[x+1 for x in upper[:-1]]
        bins=list(zip(lower,upper));pib=[sum(prior[a:z+1],F(0))if a<=z else F(0)for a,z in bins]
        grams=[integrate_gram(n,w,p,prior)for w,p in channels]
        assert grams[0]==grams[1]
        variances=[];star=[]
        for b,(a,z)in enumerate(bins):
            if not pib[b]:variances.append(F(0));star.append(F(0));continue
            def variance(g):
                return sum((prior[s]*prior[t]*g[s][t]for s in range(a,z+1)for t in range(a,z+1)),F(0))/pib[b]**2-1
            vb,vs=variance(grams[b]),variance(grams[2])
            assert vb<=vs
            variances.append(vb);star.append(vs)
        budget=sum((p*v for p,v in zip(pib,variances)),F(0))
        startrace=sum((p*v for p,v in zip(pib,star)),F(0))
        bound=fs[2]/(1-fs[2]);assert budget<=startrace<=bound
        record={'eta':str(eta),'n':n,'R':R,'f_by_table':list(map(str,fs)),
                'same_f_different_distribution_full_Gram_equal':True,
                'bins':bins,'pi':list(map(str,pib)),'variance_by_chosen_table':list(map(str,variances)),
                'weighted_actual_diagonals':str(budget),'common_star_coarsening_trace':str(startrace),
                'geometric_upper':str(bound),'monotonicity_and_trace_exact':True}
        if n==3:
            # Common base measure uniform on four output atoms. Each table has
            # its own true query density 4*w and posterior p, so A_b differ.
            l1=F(0)
            for obs in itertools.product(range(4),repeat=n):
                G=P=F(0)
                for b,((w,p),(a,z))in enumerate(zip(channels,bins)):
                    A=math.prod(4*w[o]for o in obs)
                    posterior=pb([p[o]for o in obs])
                    T=sum(posterior[a:z+1],F(0))if a<=z else F(0)
                    G+=A*T;P+=A*pib[b]
                l1+=abs(G-P)/4**n
            assert l1*l1<=budget
            record['actual_common_space_L1']=str(l1)
            record['L1_squared_le_weighted_diagonal_budget_exact']=True
        out.append(record)
target=HERE/'multitable_results.json'
if target.exists():raise RuntimeError('Refusing replacement')
result={'scope':'Fraction deterministic multitable binary product channels; no cube/MC substitution',
        'cases':out}
target.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
