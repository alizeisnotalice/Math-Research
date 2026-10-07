#!/usr/bin/env python3
"""Exact finite capture-poset width certificates; no all-input conclusion."""
from pathlib import Path
from fractions import Fraction as F
from collections import defaultdict, deque
import argparse, importlib.util, json, sys, time, hashlib
sys.dont_write_bytecode = True

def module(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m); return m

def width_certificate(sets):
    """Dilworth via exact integer matching, with antichain and chain cover."""
    sets = sorted(sets, key=lambda S: (len(S), tuple(sorted(S))))
    N = len(sets); adj = [[b for b,T in enumerate(sets) if S < T] for S in sets]
    right = [-1]*N
    def augment(a, seen):
        for b in adj[a]:
            if b in seen: continue
            seen.add(b)
            if right[b] < 0 or augment(right[b], seen):
                right[b] = a; return True
        return False
    for a in range(N): augment(a, set())
    left = [-1]*N
    for b,a in enumerate(right):
        if a >= 0: left[a] = b
    zl = {a for a in range(N) if left[a] < 0}; zr = set(); queue = deque(zl)
    while queue:
        a = queue.popleft()
        for b in adj[a]:
            if left[a] == b or b in zr: continue
            zr.add(b)
            if right[b] >= 0 and right[b] not in zl:
                zl.add(right[b]); queue.append(right[b])
    # Minimum vertex cover supplies an independently checked matching dual.
    cover_l = set(range(N))-zl; cover_r = zr
    matching = [(a,b) for a,b in enumerate(left) if b >= 0]
    assert len(cover_l)+len(cover_r) == len(matching)
    assert all(a in cover_l or b in cover_r for a in range(N) for b in adj[a])
    antichain = [sets[a] for a in range(N) if a in zl and a not in zr]
    chains = []
    for a in range(N):
        if right[a] >= 0: continue
        chain=[]
        while a >= 0: chain.append(sets[a]); a=left[a]
        chains.append(chain)
    width=N-len(matching)
    assert len(antichain)==len(chains)==width
    assert all(not (S <= T or T <= S) for a,S in enumerate(antichain) for T in antichain[a+1:])
    assert all(S < T for chain in chains for S,T in zip(chain,chain[1:]))
    assert {S for chain in chains for S in chain} == set(sets)
    return dict(width=width, antichain=list(map(sorted,antichain)),
                minimum_chain_cover=[list(map(sorted,c)) for c in chains],
                matching_edges=[dict(smaller=sorted(sets[a]),larger=sorted(sets[b])) for a,b in matching],
                vertex_cover_left=[sorted(sets[a]) for a in sorted(cover_l)],
                vertex_cover_right=[sorted(sets[b]) for b in sorted(cover_r)])

def evaluate(flow, record, alloc, old, batch, label, atoms, alpha, params, reference):
    started=time.monotonic(); saved,radii,cells,_=flow.context(old,atoms,alpha)
    I=alpha/4; groups=defaultdict(F); natural=[F(0) for _ in atoms]; witnesses={}
    for a,b,J,trace in cells:
        rec=record.intervals(trace,alpha/4,alpha/2)
        assert sum((hi-lo for lo,hi in rec.values()),F(0))==I
        x=(a+b)/2
        for k,(lo,hi) in rec.items():
            S=frozenset(i for i,(y,w) in enumerate(atoms) if abs(y-x)<radii[k])
            q=sum((atoms[i][1] for i in S),F(0)); g=q/(2*radii[k]);d=(b-a)*(hi-lo)/I*g
            assert d>0;groups[S]+=d
            witnesses.setdefault(S,dict(x_interval=[str(a),str(b)],k=k,radius=str(radii[k]),record_interval=[str(lo),str(hi)]))
            for i in S:natural[i]+=d/q
    weights=[w for y,w in atoms];M=sum(groups.values(),F(0));X=alpha*F(saved['eligible_volume'])
    assert F(3,8)*X<=M<=F(3,4)*X
    C,U,history,allocation,loads=alloc.solve(groups,weights)
    assert C<=max(natural)
    assert sum((weights[i]*natural[i] for i in range(len(weights))),F(0))==M
    per_source=[]
    for i,(position,mass) in enumerate(atoms):
        cert=width_certificate([S for S in groups if i in S])
        assert all(i in S for S in cert['antichain'])
        per_source.append(dict(source=i,position=str(position),mass=str(mass),**cert))
    width=max(row['width'] for row in per_source); strongest=next(row for row in per_source if row['width']==width)
    chains=[]
    for S in sorted(groups,key=lambda S:(-len(S),tuple(sorted(S)))):
        for chain in chains:
            if S<=chain[-1]:chain.append(S);break
        else:chains.append([S])
    root_load=[sum(i in chain[0] for chain in chains) for i in range(len(weights))]
    kappa=max(root_load); assert width<=kappa
    chain_ratios=[]
    for chain in chains:
        for S in chain:
            tail=sum((groups[T] for T in chain if T<=S),F(0));mass=sum((weights[i] for i in S),F(0))
            assert tail<=F(3,2)*mass;chain_ratios.append(tail/mass)
    assert C<=F(3,2)*kappa
    if reference:
        assert C==F(reference['optimal_C']) and kappa==reference['chain_root_congestion']
        assert M==F(reference['M']) and X==F(reference['X'])
    crossing=next(((S,T) for S in groups for T in groups if S&T and not(S<=T or T<=S)),None)
    subtree_rows=[];subtree_C=None
    if crossing is None:
        subtree_demand={}
        for S in sorted(groups,key=lambda T:len(T)):
            descendants=[T for T in groups if T<S]
            children=[T for T in descendants if not any(T<U<S for U in groups)]
            direct=sum((groups[T] for T in groups if T<=S),F(0))
            recursive=groups[S]+sum((subtree_demand[T] for T in children),F(0))
            assert direct==recursive
            subtree_demand[S]=recursive
            mass=sum((weights[i] for i in S),F(0))
            subtree_rows.append(dict(sources=sorted(S),demand=str(groups[S]),children=list(map(sorted,children)),subtree_demand=str(recursive),source_mass=str(mass),ratio=str(recursive/mass)))
        subtree_C=max(F(row['ratio']) for row in subtree_rows)
        assert subtree_C==C
    return dict(laminar=(crossing is None),subtree_C=str(subtree_C) if subtree_C is not None else None,
                subtree_rows=subtree_rows,crossing_witness=list(map(sorted,crossing)) if crossing else None,
                batch=batch,label=label,params=params,N=len(atoms),D=saved['D'],capture_classes=len(groups),
                X=str(X),M=str(M),C_alloc=str(C),natural_max=str(max(natural)),
                atom_width_lower_bound=width,greedy_kappa=kappa,width_over_C=str(F(width)/C),
                greedy_over_width=str(F(kappa,width)),width_equals_greedy=(width==kappa),
                strongest_source=strongest['source'],explicit_antichain=strongest['antichain'],
                explicit_antichain_capture_witnesses=[dict(sources=S,**witnesses[frozenset(S)]) for S in strongest['antichain']],
                per_source_width_certificates=per_source,greedy_chain_cover=[list(map(sorted,c)) for c in chains],
                greedy_root_load=root_load,chain_max_tail_ratio=str(max(chain_ratios)),
                exact_primal_rows=[dict(sources=sorted(S),demand=str(groups[S]),allocation={str(i):str(v) for i,v in row.items()}) for S,row in allocation.items()],
                tight_source_set=sorted(U),tight_source_mass=str(sum((weights[i] for i in U),F(0))),
                tight_cut_demand=str(sum((d for S,d in groups.items() if S<=U),F(0))),
                primal_source_loads=list(map(str,loads)),ratio_iteration_count=len(history),ratio_history=history,
                round56_receipt_check=bool(reference),elapsed_seconds=time.monotonic()-started)

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--repo-root',type=Path,default=Path(__file__).resolve().parents[2])
    ap.add_argument('--output',type=Path,required=True);args=ap.parse_args();root=args.repo_root
    flow=module(root/'rounds/047/flow_probe.py','flow57');record=module(root/'rounds/046/threshold_probe.py','rec57')
    prefix=module(root/'rounds/046/prefix_probe.py','pref57');lag=module(root/'rounds/054/lag_probe.py','lag57')
    obs=module(root/'rounds/053/obstruction_probe.py','obs57');alloc=module(root/'rounds/056/allocation_probe.py','alloc57')
    old=record.load(root,args.output)
    refs={c['label']:c for c in json.loads((root/'rounds/056/verification.json').read_text())['allocation']['cases']}
    inputs=[c for c in flow.cases(old,prefix) if c[1] in lag.ORIGINAL|lag.EXTRA]
    inputs.append((0,'round53_four_atom_theta1',sorted(obs.ATOMS),F(1),{}))
    # L8 is already among the original sixteen; add exactly four new cases.
    for L in (4,6,10,12):
        inputs.append((2,f'chain_L{L}_alpha3/2',old.old.old.chain(L,den=4,position=F(9,8),weight=F(27,16)),F(3,2),dict(L=L,family='chain_depth_pressure')))
    started=time.monotonic();results=[]
    # Standalone posets exercise branch, chain, isolated and empty certificates.
    for sets,expected in (([],0),([frozenset({0})],1),([frozenset({0}),frozenset({0,1})],1),
                          ([frozenset({0,1}),frozenset({0,2}),frozenset({0,3}),frozenset({0,1,2,3})],3)):
        assert width_certificate(sets)['width']==expected
    for batch,label,atoms,alpha,params in inputs:
        c=evaluate(flow,record,alloc,old,batch,label,atoms,alpha,params,refs.get(label));results.append(c)
        print(label,'width',c['atom_width_lower_bound'],'greedy',c['greedy_kappa'],'C',round(float(F(c['C_alloc'])),9),flush=True)
    assert len(results)==20 and sum(c['round56_receipt_check'] for c in results)==16
    out=dict(status='passed',dimension=1,arithmetic='Fraction exact and integer matching',cases=results,
             batch_counts=[sum(c['batch']==b for c in results) for b in range(3)],
             original_receipts_checked=16,standalone_poset_checks=4,
             source_width_certificates_checked=sum(c['N'] for c in results),
             scope='Entire eligible E, actual common-beta first-exit capture classes; per-source antichain lower bound on unweighted fractional chain-root congestion.',
             general_width_growth_proved=False,general_C_alloc_bound_proved=False,
             elapsed_seconds=time.monotonic()-started,script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    args.output.write_text(json.dumps(out,indent=2)+'\n')
    print('PASSED',len(results),out['batch_counts'],'seconds',round(out['elapsed_seconds'],3),flush=True)

if __name__=='__main__':main()
