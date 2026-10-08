#!/usr/bin/env python3
"""Exact first-exit restriction and tight-cut checks; no uniform weak bound."""
from pathlib import Path
from fractions import Fraction as F
from collections import defaultdict
import importlib.util, hashlib, json, sys, tempfile
sys.dont_write_bytecode=True
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]

def module(path,name):
    s=importlib.util.spec_from_file_location(name,path)
    m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m

def maximal(atoms,x):
    distances=sorted({abs(y-x) for y,w in atoms})
    assert distances and distances[0]>0
    return max(sum((w for y,w in atoms if abs(y-x)<=r),F(0))/(2*r) for r in distances)

def capture_data(flow,record,old,atoms,alpha):
    saved,radii,cells,_=flow.context(old,atoms,alpha)
    groups=defaultdict(F);tasks=[]
    for a,b,J,trace in cells:
        x=(a+b)/2
        for k,(lo,hi) in record.intervals(trace,alpha/4,alpha/2).items():
            S=frozenset(i for i,(y,w) in enumerate(atoms) if abs(y-x)<radii[k])
            demand=(b-a)*(hi-lo)/(alpha/4)*trace[k-1]
            assert demand>0;groups[S]+=demand
            tasks.append((x,J,trace,k,lo,hi,S))
    return saved,radii,groups,tasks

def cut_audit(groups,weights,allocation,C):
    N=len(weights)
    def mass(U):return sum((weights[i] for i in U),F(0))
    def demand(U):return sum((d for S,d in groups.items() if S<=U),F(0))
    sets=[frozenset(i for i in range(N) if mask>>i&1) for mask in range(1<<N)]
    tight=[U for U in sets if demand(U)==C*mass(U)]
    tight_set=set(tight);T=frozenset().union(*tight)
    assert max((demand(U)/mass(U) for U in sets if U),default=F(0))==C
    pair_checks=0
    for U in tight:
        for V in tight:
            assert U|V in tight_set and U&V in tight_set
            assert not any(d>0 and S<=U|V and not S<=U and not S<=V for S,d in groups.items())
            pair_checks+=1
    signatures={i:tuple(i in U for U in tight) for i in T}
    blocks=[]
    for sig in sorted(set(signatures.values())):
        blocks.append(frozenset(i for i in T if signatures[i]==sig))
    predecessors=[]
    for B in blocks:
        j=next(iter(B))
        predecessors.append({a for a,A in enumerate(blocks) if all(j not in U or A<=U for U in tight)})
    ideals=set()
    for mask in range(1<<len(blocks)):
        chosen={i for i in range(len(blocks)) if mask>>i&1}
        if all(predecessors[i]<=chosen for i in chosen):
            ideals.add(frozenset().union(*(blocks[i] for i in chosen)))
    assert ideals==tight_set
    totals=[F(0)]*len(blocks);edge_checks=0;projections=[]
    for S,d in groups.items():
        if not S<=T:
            assert all(allocation[S].get(i,F(0))==0 for i in T)
            continue
        hit={a for a,B in enumerate(blocks) if B&S}
        top=[a for a in hit if hit<=predecessors[a]]
        assert len(top)==1
        a=top[0];totals[a]+=d
        if not S<=blocks[a]:projections.append(dict(sources=sorted(S),top_block=sorted(blocks[a]),demand=str(d)))
        assert all(v==0 or i in blocks[a] for i,v in allocation[S].items())
        edge_checks+=1
    assert all(totals[a]==C*mass(B) for a,B in enumerate(blocks))
    minimal=[U for U in tight if U and not any(V and V<U for V in tight)]
    return dict(tight_count=len(tight),block_count=len(blocks),pair_checks=pair_checks,
                top_block_edges=edge_checks,proper_projections=projections,minimal_tight_sets=list(map(sorted,minimal)),
                maximal_tight_set=sorted(T),minimal_cuts_cover_maximal=frozenset().union(*minimal)==T)

def main():
    flow=module(ROOT/'rounds/047/flow_probe.py','flow58')
    record=module(ROOT/'rounds/046/threshold_probe.py','record58')
    prefix=module(ROOT/'rounds/046/prefix_probe.py','prefix58')
    lag=module(ROOT/'rounds/054/lag_probe.py','lag58')
    obs=module(ROOT/'rounds/053/obstruction_probe.py','obs58')
    alloc=module(ROOT/'rounds/056/allocation_probe.py','alloc58')
    with tempfile.TemporaryDirectory(prefix='euclidean58_') as folder:
        old=record.load(ROOT,Path(folder)/'scratch.json')
        inputs=[c for c in flow.cases(old,prefix) if c[1] in lag.ORIGINAL|lag.EXTRA]
        inputs.append((0,'round53_four_atom_theta1',sorted(obs.ATOMS),F(1),{}))
        for L in (4,6,10,12):
            inputs.append((2,f'chain_L{L}_alpha3/2',old.old.old.chain(L,den=4,position=F(9,8),weight=F(27,16)),F(3,2),{}))
        reference={c['label']:c for c in json.loads((ROOT/'rounds/057/verification.json').read_text())['poset']['cases']}
        rows=[]
        for batch,label,atoms,alpha,_ in inputs:
            saved,radii,groups,tasks=capture_data(flow,record,old,atoms,alpha)
            weights=[w for y,w in atoms];C,U,_,allocation,_=alloc.solve(groups,weights)
            assert C==F(reference[label]['C_alloc'])
            candidates=set(groups)|{U,frozenset(range(len(atoms)))}
            if len(atoms)<=5:
                candidates|={frozenset(i for i in range(len(atoms)) if mask>>i&1) for mask in range(1,1<<len(atoms))}
            checks=changed_J=0
            for x,J,trace,k,lo,hi,S in tasks:
                mp=maximal(atoms,x);beta=(lo+hi)/2
                kp=next(i for i,g in enumerate(trace,1) if g>alpha/2)
                assert alpha<mp<=2*alpha
                for V in candidates:
                    if not S<=V:continue
                    Q=[atoms[i] for i in sorted(V)]
                    tq=[sum((w for y,w in Q if abs(y-x)<radii[j]),F(0))/(2*radii[j]) for j in range(1,len(trace)+1)]
                    jq=next(i for i,g in enumerate(tq,1) if g>alpha/8)
                    assert J<=jq<=k
                    assert next(i for i,g in enumerate(tq,1) if g>beta)==k
                    assert next(i for i,g in enumerate(tq,1) if g>alpha/2)==kp
                    assert tq[k-1:]==trace[k-1:] and maximal(Q,x)==mp
                    # Remote ballast is checked at this same representative point.
                    p=sum((w for y,w in Q),F(0))
                    if p<1:
                        distant=max(abs(y) for y,w in atoms)+abs(x)+radii[1]+100
                        padded=Q+[(distant,1-p)]
                        assert maximal(padded,x)==mp
                    checks+=1;changed_J+=jq!=J
            cut=cut_audit(groups,weights,allocation,C) if len(atoms)<=10 else None
            rows.append(dict(batch=batch,label=label,N=len(atoms),D=saved['D'],C=str(C),
                             task_cells=len(tasks),restriction_checks=checks,changed_J=changed_J,cut_structure=cut))
        degeneration=[];replication=[]
        for batch,powers in enumerate(((6,8,10),(12,16,20),(24,32,48))):
            for power in powers:
                eps=F(1,2**power);delta=F(5,16)-eps
                atoms=[(F(0),F(1,8)),(delta,F(1,8)),(F(10),F(3,4))]
                saved,radii,groups,tasks=capture_data(flow,record,old,atoms,F(1))
                expected={frozenset({0}):F(1,32)-eps/2,
                          frozenset({1}):F(1,32)-eps/2,frozenset({0,1}):eps}
                assert dict(groups)==expected and F(saved['eligible_volume'])==F(1,8)
                C,U,_,allocation,_=alloc.solve(groups,[w for y,w in atoms])
                assert C==F(1,4) and U==frozenset({0,1})
                cert=cut_audit(groups,[w for y,w in atoms],allocation,C)
                assert cert['tight_count']==2 and cert['minimal_tight_sets']==[[0,1]]
                phi=eps/(C*F(1,8));assert phi==32*eps
                assert {k for x,J,tr,k,lo,hi,S in tasks}=={4,5}
                assert all(next(i for i,g in enumerate(tr,1) if g>F(1,2))==6 for x,J,tr,k,lo,hi,S in tasks)
                degeneration.append(dict(batch=batch,power=power,epsilon=str(eps),source_gap=str(delta),
                                         singleton_demand=str(F(1,32)-eps/2),crossing_demand=str(eps),
                                         C=str(C),normalized_crossing=str(phi),task_cells=len(tasks)))
        # Canonical normalization counterexample and separated-copy checks.
        q=F(1,2);x=F(3,16)
        assert maximal([(F(0),q)],x)==F(4,3)
        def trace_at(atoms,x,alpha,D):
            return [sum((w for y,w in atoms if abs(y-x)<F(4,2**j)/alpha),F(0))/(F(8,2**j)/alpha) for j in range(1,D+1)]
        tq=trace_at([(F(0),q)],x,F(1),8)
        tn=trace_at([(F(0),F(1))],x,F(2),8)
        assert next(i for i,g in enumerate(tq,1) if g>F(1,8))==2
        assert next(i for i,g in enumerate(tn,1) if g>F(2,8))==1
        for batch,ms in enumerate(((0,1,2),(3,4,5),(6,7,8))):
            for m in ms:
                q=F(3,10);lam=F(1,2**m);N=int(1/(lam*q));gap=F(100)/lam
                atoms=[(i*gap,lam*q) for i in range(N)]
                remainder=1-N*lam*q
                if remainder:atoms.append(((N+2)*gap,remainder))
                x=3*q/8;alpha=lam
                tr=trace_at(atoms,x,alpha,10+m)
                original=trace_at([(F(0),q)],x,F(1),10)
                assert tr[m:]==[lam*g for g in original]
                assert all(g<=alpha/8 for g in tr[:m])
                assert maximal(atoms,x)==lam*maximal([(F(0),q)],x)
                assert next(i for i,g in enumerate(tr,1) if g>alpha/8)==m+next(i for i,g in enumerate(original,1) if g>F(1,8))
                assert sum(w for y,w in atoms)==1
                replication.append(dict(batch=batch,m=m,copies=N,filled_mass=str(N*lam*q),remainder=str(remainder)))
        # Abstract nested-tight example: NOT a Euclidean realization.
        gs={frozenset({0}):F(1,2),frozenset({0,1}):F(1,2)};ws=[F(1,2),F(1,2)]
        C,U,_,a,_=alloc.solve(gs,ws);nested=cut_audit(gs,ws,a,C)
        assert not nested['minimal_cuts_cover_maximal'] and nested['tight_count']==3
    assert len(rows)==20
    deps=sorted((ROOT/'runtime').rglob('*.py'))+[ROOT/'rounds'/s for s in
        ('041/square_probe.py','042/row_probe.py','043/global_probe.py','044/tail_probe.py',
         '046/prefix_probe.py','046/threshold_probe.py','047/flow_probe.py',
         '053/obstruction_probe.py','054/lag_probe.py','056/allocation_probe.py','057/verification.json')]+[Path(__file__).resolve()]
    result=dict(status='passed',round=58,date='2026-10-08',python=sys.version.split()[0],
                restriction=dict(batch_counts=[sum(c['batch']==b for c in rows) for b in range(3)],
                                 checks=sum(c['restriction_checks'] for c in rows),
                                 changed_J=sum(c['changed_J'] for c in rows),
                                 exhaustive_cut_cases=sum(c['cut_structure'] is not None for c in rows),cases=rows),
                actual_expansion_counterexample=dict(batch_counts=[3,3,3],cases=degeneration),
                abstract_nested_tight_example=nested,replication_checks=replication,normalization_clock_counterexample=dict(before_J=2,after_J=1),
                scope='Exact n=1 full-E capture demands; restriction checks at rational cell representatives. General restriction and all-epsilon counterexample proved analytically in report; no high-dimensional integration or continuous-input optimal-cut claim.',
                main_weak_type_theorem_proved=False,
                sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in deps})
    print(json.dumps(result,ensure_ascii=False,indent=2))
if __name__=='__main__':main()
