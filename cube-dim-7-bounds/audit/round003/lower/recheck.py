#!/usr/bin/env python3
"""Read-only exact R003 winner audit; stdout only, no author code imported.

Requires --candidate-root and --reference paths. The latter is the complete R002 results.json.
Reuses the independent R002 full labelled-box fixed-coordinate sweep algorithm.
"""
from pathlib import Path
from fractions import Fraction as Q
from itertools import product
import argparse, hashlib, json, math, platform, sys, time

HERE=Path(__file__).resolve().parent
D=10**35; S=4*D; LIMIT=2**20
def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def need(condition,message):
    if not condition: raise ValueError(message)
def root(n):
    x=1 << ((n.bit_length()+6)//7)
    while True:
        y=(6*x+n//x**6)//7
        if y>=x: break
        x=y
    while x**7>n: x-=1
    while (x+1)**7<=n: x+=1
    need(x**7<=n<(x+1)**7,'Integer seventh-root isolation failed.')
    return x
def add(a,b): return (a[0]+b[0],a[1]+b[1])
def mul(a,c): return tuple(sorted((a[0]*c,a[1]*c)))
def maximum(a,b): return (max(a[0],b[0]),max(a[1],b[1]))
def minimum(a,b): return (min(a[0],b[0]),min(a[1],b[1]))
def units(x):
    y=x*S
    need(y.denominator==1,'Nonintegral common-scale endpoint.')
    return y.numerator

class RawUnion:
    def __init__(self,budget): self.memo={}; self.budget=budget
    def volume(self,boxes):
        if not boxes: return 0
        if boxes in self.memo: return self.memo[boxes]
        self.budget[0]+=1
        need(self.budget[0]<=LIMIT,'Registered state limit reached; incomplete.')
        dims=len(boxes[0])
        need(all(len(b)==dims for b in boxes),'Mixed-dimensional state.')
        if dims==1:
            total=0; right=None
            for ((a,b),) in sorted(boxes):
                if right is None or a>right: total+=b-a; right=b
                elif b>right: total+=b-right; right=b
        else:
            events={}
            for box in boxes:
                a,b=box[0]
                events.setdefault(a,[set(),set()])[0].add(box)
                events.setdefault(b,[set(),set()])[1].add(box)
            active=set(); last=None; total=0
            for x,(starts,ends) in sorted(events.items()):
                if last is not None and active:
                    section=tuple(sorted({b[1:] for b in active}))
                    total+=(x-last)*self.volume(section)
                active.difference_update(ends); active.update(starts); last=x
            need(not active,'Event sweep did not close.')
        self.memo[boxes]=total
        return total

def construct(claim):
    need(claim['event']['key']==[[4,1],[432,1]],'Outside the single registered winning instance.')
    lengths=list(product((1,2),*[(1,2,3,4)]*6))
    roots={m:root(m*D**7) for m in {math.prod(e) for e in lengths}}
    need({str(m):str(v) for m,v in roots.items()}==claim['root_lower_numerators'],'Saved mass-root mismatch.')
    event_roots={m:root(m*D**7) for m in (4,432)}
    need({str(m):str(v) for m,v in event_roots.items()}==claim['event_radical_root_lower_numerators'],'Saved event-root mismatch.')
    lo=sum(event_roots.values()); hi=lo+2
    event=claim['event']
    need(int(event['a1_lower_numerator'])==lo and int(event['a1_upper_numerator'])==hi and int(event['a1_denominator'])==2*D,'Saved event spacing enclosure mismatch.')
    first=(Q(lo,2*D),Q(hi,2*D))
    tail=(Q(D+roots[2],2*D),Q(D+roots[2]+1,2*D))
    inward=[]; outward=[]; ikeys=[]; okeys=[]; checks=0; minwidth=None
    for ell in lengths:
        mass=math.prod(ell); radius=(Q(roots[mass],2*D),Q(roots[mass]+1,2*D))
        ib=[]; ob=[]
        for j,count in enumerate(ell):
            spacing=first if j==0 else tail
            start=-(count//2) if count%2 else 1-count//2
            end=start+count-1
            left=maximum((Q(0),Q(0)),add(mul(spacing,end),mul(radius,-1)))
            right=minimum(mul(spacing,Q(1,2)),add(mul(spacing,start),radius))
            ib.append((units(left[1]),units(right[0]))); ob.append((units(left[0]),units(right[1])))
        if all(a<b for a,b in ob): outward.append(tuple(ob)); okeys.append(ell)
        if all(a<b for a,b in ib):
            inward.append(tuple(ib)); ikeys.append(ell)
            for j,count in enumerate(ell):
                spacing=first if j==0 else tail
                start=-(count//2) if count%2 else 1-count//2
                l,h=(Q(x,S) for x in ib[j])
                for index in range(start,start+count):
                    atom=mul(spacing,index)
                    need(l>=atom[1]-radius[0] and h<=atom[0]+radius[0],'Actual source index containment failed.')
                    checks+=2
                minwidth=h-l if minwidth is None else min(minwidth,h-l)
    need(len(ikeys)==claim['positive_inward_raw_boxes'] and len(okeys)==claim['positive_outward_raw_boxes'],'Positive raw-box count mismatch.')
    for boxes in (inward,outward):
        known=set(boxes)
        for b in boxes:
            for j in range(1,6):
                c=list(b); c[j],c[j+1]=c[j+1],c[j]
                need(tuple(c) in known,'Full labelled tail permutation closure failed.')
    return tuple(sorted(set(inward))),tuple(sorted(set(outward))),len(ikeys),len(okeys),checks,minwidth,max(first[1],tail[1])

def recheck(candidate,reference=None):
    t=time.time(); cpu=time.process_time(); candidate=Path(candidate).expanduser().resolve()
    lock=json.loads((HERE/'input_lock.json').read_text())
    before={n:sha(candidate/n) for n in lock['candidate_files']}
    need(before==lock['candidate_files'],'Candidate differs from final audited source lock.')
    frozen=json.loads((candidate/'frozen_inputs.json').read_text())
    interval=tuple(map(Q,frozen['complete_R002_interval']))
    need(reference is not None,'R002 reference required')
    reference_path=Path(reference).expanduser().resolve(); reference_hash=sha(reference_path)
    need(reference_hash==lock['reference_sha256'],'R002 reference hash mismatch')
    ref=json.loads(reference_path.read_text())['full_construction_comparison_interval']
    need(tuple(map(Q,ref))==interval,'Frozen and actual R002 intervals differ')
    comparison={'reference_sha256':reference_hash,'complete_R002_interval':[str(x) for x in interval]}
    claim=json.loads((candidate/'results.json').read_text())
    need(Q(claim['complete_round002_upper'])==interval[1],'Saved comparison did not use the complete R002 upper.')
    inner,outer,ni,no,checks,minwidth,maxspacing=construct(claim)
    budget=[0]; values=[]; counts=[]
    for boxes in (inner,outer):
        engine=RawUnion(budget); values.append(Q(128*engine.volume(boxes),S**7)); counts.append(len(engine.memo))
    need(values==list(map(Q,claim['full_construction_comparison_interval'])),'Independent exact union mismatch.')
    need(values[0]<=values[1] and values[0]>interval[1],'Union direction or strict complete-R002 comparison failed.')
    need(Q(claim['strict_margin_vs_complete_round002_upper'])==values[0]-interval[1],'Saved strict ideal margin mismatch.')
    need(Q(claim['minimum_inward_width'])==minwidth,'Saved minimum inward width mismatch.')
    need(claim['source_extreme_checks']==14*ni,'Saved source extreme-check count mismatch.')
    eps=Q(1,10**25); eta=eps/2; theta=1-Q(1,10**18); K=10**19
    need(minwidth>2*eps and maxspacing/2<2,'Finite L1 geometric budget failed.')
    loss=128*ni*14*eps*2**6; finite=theta*(values[0]-loss)*Q(2*K+1,2*K+5)**7
    fc=claim['finite_L1']
    need((Q(fc['epsilon']),Q(fc['eta']),Q(fc['theta']),int(fc['K']),int(fc['source_collar']))==(eps,eta,theta,K,2),'Finite L1 parameter mismatch.')
    need(loss==Q(fc['union_shrink_loss_upper']) and finite==Q(fc['exact_ratio_lower']),'Finite L1 exact fraction mismatch.')
    need(finite>interval[1] and finite>Q(2136429741978,10**12),'Strict finite L1 thresholds failed.')
    need(Q(fc['strict_margin_vs_complete_round002_upper'])==finite-interval[1],'Saved strict finite margin mismatch.')
    need({n:sha(candidate/n) for n in before}==before,'Candidate changed during run.')
    if reference_path is not None: need(sha(reference_path)==reference_hash,'Reference changed during run.')
    return {'status':'PASS independent finite winner audit','candidate_root':str(candidate),'comparison':comparison,
        'recheck_sha256':sha(__file__),'input_lock_sha256':sha(HERE/'input_lock.json'),
        'scope':'one winning instance, all 8192 lengths; no old baseline evaluations or new search',
        'algorithm':'full labelled raw boxes, fixed-coordinate event sweep; independent integer Newton/Fraction construction',
        'positive_inward_raw_boxes':ni,'positive_outward_raw_boxes':no,'distinct_inward_boxes':len(inner),'distinct_outward_boxes':len(outer),
        'all_source_coordinate_index_checks':checks,'memo_states_each_direction':counts,'memo_states_both_directions':budget[0],
        'union_lower':str(values[0]),'union_upper':str(values[1]),'finite_L1_exact_lower':str(finite),
        'strict_ideal_margin_vs_complete_R002_upper':str(values[0]-interval[1]),'strict_finite_margin_vs_complete_R002_upper':str(finite-interval[1]),
        'proved_simple_strict_bound':'C7 > 2136429741978/10^12','all_source_hashes_unchanged_after_run':True,
        'wall_seconds':time.time()-t,'cpu_seconds':time.process_time()-cpu,
        'runtime':{'executable':sys.executable,'python':sys.version,'platform':platform.platform()}}

class Parser(argparse.ArgumentParser):
    def error(self,message): raise ValueError(message)
def main():
    p=Parser(description=__doc__)
    p.add_argument('--candidate-root',type=Path,required=True)
    p.add_argument('--reference',type=Path,required=True,help='complete original R002 results.json')
    try:
        args=p.parse_args(); print(json.dumps(recheck(args.candidate_root,args.reference),indent=2),flush=True); return 0
    except Exception as exc:
        print(json.dumps({'status':'FAIL independent finite winner audit','error':str(exc)}),flush=True); return 1
if __name__=='__main__': sys.exit(main())
