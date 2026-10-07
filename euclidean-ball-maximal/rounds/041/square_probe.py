#!/usr/bin/env python3
"""Exact actual-band weighted square search. No stop simulation or artificial masks."""
from fractions import Fraction as F
from collections import defaultdict
from bisect import bisect_left
from pathlib import Path
import sys,json,time,random
sys.dont_write_bytecode=True
RUNTIME=Path(__file__).resolve().parents[2]/'runtime/work/general_input_20261003'
sys.path.insert(0,str(RUNTIME))
from stage24_clock_replacement_probe import observer
from stage30_transition_weight_probe import build_h,Linear
from stage40_band_search import chain
import argparse
_parser=argparse.ArgumentParser()
_parser.add_argument('--output',type=Path,required=True)
OUT=_parser.parse_args().output

def finish(raw,alpha=F(1)):
    d=defaultdict(F)
    for x,w in raw:d[x/alpha]+=w
    assert 0<sum(d.values())<1
    d[F(10)/alpha]+=1-sum(d.values())
    return sorted(d.items())

def tree(depth,q=8,weight=F(9,8),position=F(3,2),jitter=False,seed=0,leaf=F(1,16),skew=False):
    rng=random.Random(seed);raw=[]
    def visit(x,r,n):
        w=r*(leaf if n==0 else weight)
        if jitter:w*=F(rng.randint(48,80),64)
        raw.append((x,w))
        if n:
            visit(x-position*r,r/q,n-1)
            visit(x+position*r,r/(q*(2 if skew else 1)),n-1)
    visit(F(0),F(1,8),depth)
    return raw

def cases():
    # Benchmarks only; the three new batches below are the substantive search.
    yield 0,'old40_optimized_L4',chain(4,den=16,weight=F(27,16),position=F(9,8)),F(1),{}
    yield 0,'old40_chain_L24',chain(24,den=16),F(1),{}
    for depth,q,skew in [(2,8,False),(3,8,False),(4,8,False),(5,8,False),(4,6,False),(4,8,True)]:
        yield 1,f'nested_d{depth}_q{q}_skew{skew}',finish(tree(depth,q,skew=skew)),F(1),dict(depth=depth,q=q,skew=skew)
    for depth,alpha,position,weight in [(3,F(3,4),F(9,8),F(27,16)),(4,F(5,4),F(9,8),F(27,16)),(4,F(1),F(3,2),F(9,8))]:
        raw=tree(depth,position=position,weight=weight,jitter=True,seed=410727+depth)
        yield 2,f'mirror_tree_d{depth}_a{alpha}',finish(raw,alpha),alpha,dict(depth=depth,seed=410727+depth)
    for L,den in [(12,4),(24,16),(40,16)]:
        yield 2,f'mirrored_chain_L{L}_den{den}',chain(L,den=den,mirror=True,jitter=True,seed=410741+L),F(1),dict(L=L,den=den)
    for depth,seed in [(3,410801),(4,410803),(5,410807)]:
        raw=tree(depth,q=8,weight=F(27,16),position=F(9,8),jitter=True,seed=seed,leaf=F(1,4))
        yield 3,f'random_weights_d{depth}_s{seed}',finish(raw),F(1),dict(depth=depth,seed=seed)
    for L,position,weight in [(12,F(1),F(1)),(24,F(129,128),F(129,128)),(40,F(127,128),F(127,128)),(48,F(9,8),F(27,16))]:
        yield 3,f'threshold_L{L}_p{position}_w{weight}',chain(L,den=4,position=position,weight=weight),F(1),dict(L=L,position=str(position),weight=str(weight))

def evaluate(atoms,alpha,label):
    saved,en,R,H=observer(atoms,alpha,'band',4)
    assert F(saved['entrance_a'])==alpha/8 and F(saved['terminal_b'])==alpha/2 and saved['b_over_a']==4
    cells=[dict(v,radius=str(R[v['K']])) for v in saved['observer_cells']]
    hs=build_h(cells,F(saved['entrance_a']),'eligible');events=defaultdict(F)
    for h in hs.values():
        for x,d in h.events.items():events[x]+=d
    A=Linear(events);ks=sorted(hs);vectors=[(y,w,[hs[k].value(y) for k in ks]) for y,w in atoms]
    assert all(0<=v<=1 for y,w,vs in vectors for v in vs)
    M=sum((w*sum(vs) for y,w,vs in vectors),F(0))
    Q=sum((w*sum(vs)**2 for y,w,vs in vectors),F(0))
    diagonal=sum((w*sum(v*v for v in vs) for y,w,vs in vectors),F(0))
    gaps=defaultdict(F)
    for y,w,vs in vectors:
        for i,k in enumerate(ks):
            for j in range(i+1,len(ks)):gaps[ks[j]-k]+=2*w*vs[i]*vs[j]
    V=F(saved['eligible_volume']);X=alpha*V
    assert Q==diagonal+sum(gaps.values())
    # Tonelli terminal averages furnish a useful built-in independent mass ledger.
    integral=sum((F(v['hi'])-F(v['lo']))*sum(w for y,w in atoms if abs((F(v['lo'])+F(v['hi']))/2-y)<R[v['K']])/(2*R[v['K']]) for v in cells)
    assert M==integral and X/2<=M<=X
    summary=dict(label=label,N=len(atoms),D=saved['D'],cell_count=len(cells),K_count=len(ks),
        alpha=str(alpha),a=saved['entrance_a'],b=saved['terminal_b'],b_over_a=saved['b_over_a'],
        M=str(M),Q=str(Q),X=str(X),Q_over_X=str(Q/X) if X else None,
        diagonal=str(diagonal),offdiagonal=str(Q-diagonal),
        maxA=str(max((A.value(x) for x in A.knots),default=F(0))),
        max_source_A=str(max((sum(vs) for y,w,vs in vectors),default=F(0))),
        gap_mass={str(d):str(v) for d,v in sorted(gaps.items()) if v},
        gap_mass_over_X={str(d):str(v/X) for d,v in sorted(gaps.items()) if v},
        original_volume=saved['original_volume'],eligible_volume=str(V),
        displays=dict(M=float(M),Q=float(Q),X=float(X),Q_over_X=float(Q/X) if X else None,
                      maxA=float(max((A.value(x) for x in A.knots),default=F(0)))))
    return summary,dict(source=[dict(location=str(y),mass=str(w)) for y,w in atoms],alpha=str(alpha),observer_cells=cells,
        original_intervals=saved['original_intervals'],A_knots=A.save(),summary=summary)

def independent(cert):
    """Rebuild full band and clock cells from pair superlevel endpoints, then direct overlaps."""
    atoms=[(F(v['location']),F(v['mass'])) for v in cert['source']];alpha=F(cert['alpha']);s=cert['summary'];D=s['D']+1
    radii={j:F(4,2**j)/alpha for j in range(1,D+1)}
    cuts={y for y,w in atoms}|{y+sign*r for y,w in atoms for r in radii.values() for sign in [-1,1]}
    for i in range(len(atoms)):
        m=F(0)
        for j in range(i,len(atoms)):
            m+=atoms[j][1]
            for beta in [alpha,2*alpha]:cuts.update([atoms[j][0]-m/(2*beta),atoms[i][0]+m/(2*beta)])
    def distance_mass(x):
        ds=defaultdict(F)
        for y,w in atoms:ds[abs(y-x)]+=w
        assert min(ds)>0
        points=sorted(ds);prefix=[F(0)];v=F(0)
        for r in points:prefix.append(prefix[-1]+ds[r]);v=max(v,prefix[-1]/(2*r))
        return v,points,prefix
    cells=[];volume=F(0)
    for l,r in zip(sorted(cuts),sorted(cuts)[1:]):
        x=(l+r)/2
        mx,points,prefix=distance_mass(x)
        if not alpha<mx<=2*alpha:continue
        volume+=r-l
        g={j:prefix[bisect_left(points,rad)]/(2*rad) for j,rad in radii.items()}
        J=next(j for j in g if g[j]>alpha/8);K=next(j for j in g if g[j]>alpha/2)
        if J>=2:cells.append((l,r,K))
    def A(x):return sum(max(F(0),min(r,x+radii[k])-max(l,x-radii[k]))/(2*radii[k]) for l,r,k in cells)
    M=sum(w*A(y) for y,w in atoms);Q=sum(w*A(y)**2 for y,w in atoms);X=alpha*sum(r-l for l,r,k in cells)
    assert str(M)==s['M'] and str(Q)==s['Q'] and str(X)==s['X'] and volume==F(s['original_volume'])
    return dict(status='passed',independent_full_band_and_J_K=True,audit_D=D,M=str(M),Q=str(Q),X=str(X),Q_over_X=str(Q/X))

def main():
    start=time.monotonic();out=dict(status='running',dimension=1,arithmetic='Fraction exact',scope='full actual original band; eligible J>=2',
        candidate='integral A^2 dP <= C alpha V, dimension independent',stopping_simulated=False,cases=[],retained={})
    best=None
    for batch,label,atoms,alpha,params in cases():
        t=time.monotonic();s,c=evaluate(atoms,alpha,label);s.update(batch=batch,parameters=params,elapsed_seconds=time.monotonic()-t)
        s['source']=c['source'];out['cases'].append(s)
        if best is None or F(s['Q_over_X'])>F(best['summary']['Q_over_X']):best=c
        if batch==0:out['retained'][label]=c
        out['retained']['maximum_Q_over_X']=best
        OUT.write_text(json.dumps(out)+'\n')
        print(batch,label,s['displays'],flush=True)
    out['independent_maximum_verification']=independent(best)
    out['gap_diagnostics']=[dict(label=s['label'],offdiagonal_over_X=str(F(s['offdiagonal'])/F(s['X'])),
        gap_tail_over_X={str(d):str(sum((F(v) for k,v in s['gap_mass'].items() if int(k)>=d),F(0))/F(s['X'])) for d in [1,3,6,12,24,48,96]},
        max_two_to_gap_mass_over_X=str(max((2**int(k)*F(v)/F(s['X']) for k,v in s['gap_mass'].items()),default=F(0)))) for s in out['cases']]
    out.update(status='passed',new_batches=3,new_case_count=sum(s['batch']>0 for s in out['cases']),
        elapsed_seconds=time.monotonic()-start,uniform_bound_proved=False,counterexample_to_uniform_bound_found=False,
        maximum_Q_over_X=best['summary'])
    OUT.write_text(json.dumps(out,indent=2)+'\n');print('saved',OUT,'seconds',out['elapsed_seconds'],flush=True)
if __name__=='__main__':main()
