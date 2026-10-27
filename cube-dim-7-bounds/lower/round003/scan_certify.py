#!/usr/bin/env python3
"""R003: registered continuation event scan and full rational lattice certification.
Union and tail geometry reused from accepted R002; no old baseline execution.
"""
from pathlib import Path
from itertools import product
from math import prod, factorial
from fractions import Fraction as F
import json,hashlib,time,os,sys
HERE=Path(__file__).resolve().parent
D=10**35; S=4*D; LIMIT=2**20
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def root(m):
    a,b=0,1
    while b**7<=m:b*=2
    while b-a>1:
        c=(a+b)//2
        if c**7<=m:a=c
        else:b=c
    assert a**7<=m<(a+1)**7
    return a
def write(name,obj):
    with (HERE/name).open('x') as f:json.dump(obj,f,indent=2);f.write('\n')
class Exhausted(Exception): pass

class Union:
    def __init__(self, intervals, budget):
        self.iv=intervals; self.cache={}; self.budget=budget; self.stack=[]
    def call(self,state):
        if not state: return 0
        if state in self.cache: return self.cache[state]
        if self.budget['states']>=LIMIT: raise Exhausted()
        self.budget['states']+=1
        self.stack.append(state)
        k=len(state[0]); assert all(len(b)==k for b in state)
        if k==1:
            ints=sorted(self.iv[b[0]] for b in state)
            a,b=ints[0]; total=0
            for c,d in ints[1:]:
                if b<c: total+=b-a; a,b=c,d
                else: b=max(b,d)
            total+=b-a
        else:
            labels={t for box in state for t in box}
            cuts=sorted({z for t in labels for z in self.iv[t]})
            total=0
            for a,b in zip(cuts,cuts[1:]):
                child=set()
                for box in state:
                    for t in set(box):
                        l,h=self.iv[t]
                        if l<=a and b<=h:
                            z=list(box); z.remove(t); child.add(tuple(z))
                if child: total+=(b-a)*self.call(tuple(sorted(child)))
        self.cache[state]=total
        self.stack.pop()
        return total

def tail_build():
    masses={prod(e) for e in product([1,2],*[range(1,5)]*6)}
    rr={m:root(m*D**7) for m in masses}; u=rr[2]
    aa=(2*(D+u),2*(D+u+1)); records=[]
    for e in product([1,2],*[range(1,5)]*6):
        m=prod(e); rad=(2*rr[m],2*(rr[m]+1)); inner=[];outer=[]
        for ell in e[1:]:
            q=-(ell//2) if ell%2 else 1-ell//2;c=q+ell-1
            l=(max(0,c*aa[0]-rad[1]),max(0,c*aa[1]-rad[0]))
            h=(min(aa[0]//2,q*aa[1]+rad[0]),min(aa[1]//2,q*aa[0]+rad[1]))
            inner.append((l[1],h[0]));outer.append((l[0],h[1]))
        if all(l<h for l,h in outer):
            records.append({'ell':e,'mass':m,'radius':rad,'inner':tuple(inner),'outer':tuple(outer)})
    orbit=[r for r in records if tuple(sorted(r['ell'][1:]))==r['ell'][1:]]
    assert len(orbit)==34
    assert len(records)==sum(factorial(6)//prod(factorial(r['ell'][1:].count(v)) for v in set(r['ell'][1:])) for r in orbit)
    return rr,aa,records,orbit

def first_interval(ell,aa,rr,mode):
    low=(0,0) if ell==1 else (max(0,aa[0]-rr[1]),max(0,aa[1]-rr[0]))
    high=(min(aa[0]//2,rr[0]),min(aa[1]//2,rr[1]))
    return (low[1],high[0]) if mode=='inner' else (low[0],high[1])

def event_compute(event,orbit,engines,ivs,ids):
    aa=(2*int(event['a1_lower_numerator']),2*int(event['a1_upper_numerator']))
    answer={}
    for mode in ('inner','outer'):
        boxes=[]
        for r in orbit:
            first=first_interval(r['ell'][0],aa,r['radius'],mode)
            if first[0]<first[1] and all(l<h for l,h in r[mode]):
                boxes.append((first,tuple(sorted(ids[mode][z] for z in r[mode]))))
        cuts=sorted({z for pair,_ in boxes for z in pair});total=0
        for lo,hi in zip(cuts,cuts[1:]):
            state=tuple(sorted({b for (l,h),b in boxes if l<=lo and hi<=h}))
            total+=(hi-lo)*engines[mode].call(state)
        answer[mode]=F(128*total,S**7)
    assert answer['inner']<=answer['outer']
    return aa,answer

def certify_winner(event,aa,answer,roots,tail,records,budget,t0,c0):
    positive=[];outer_n=0;checks=0;width=[]
    for r in records:
        inside=(first_interval(r['ell'][0],aa,r['radius'],'inner'),)+r['inner']
        outside=(first_interval(r['ell'][0],aa,r['radius'],'outer'),)+r['outer']
        if all(l<h for l,h in outside):outer_n+=1
        if not all(l<h for l,h in inside):continue
        positive.append(r)
        for i,(ell,(l,h)) in enumerate(zip(r['ell'],inside)):
            a=aa if i==0 else tail;q=-(ell//2) if ell%2 else 1-ell//2
            assert l>=(q+ell-1)*a[1]-r['radius'][0]
            assert h<=q*a[1]+r['radius'][0]
            width.append(h-l);checks+=2
    eps=F(1,10**25);eta=eps/2;theta=1-F(1,10**18);K=10**19
    assert F(min(width),S)>2*eps
    assert max(aa[1],tail[1])<4*S
    loss=128*len(positive)*14*eps*2**6
    finite=theta*(answer['inner']-loss)*F(2*K+1,2*K+5)**7
    ref=json.loads((HERE/'frozen_inputs.json').read_text())
    baseline=F(ref['complete_R002_interval'][1])
    assert answer['inner']>baseline and finite>baseline
    result={'status':'complete new construction, search stopped at first strict winner','dimension':7,'event':event,'first_spacing_definition':'2*a1=sum(c*m^(1/7)) for (m,c) in event.key; positive roots','coordinate_scale':str(S),'root_denominator':str(D),'root_lower_numerators':{str(k):str(v) for k,v in roots.items()},'event_radical_root_lower_numerators':{str(m):str(root(int(m)*D**7)) for m,c in event['key']},'full_construction_comparison_interval':[str(answer['inner']),str(answer['outer'])],'full_display':[float(answer['inner']),float(answer['outer'])],'complete_round002_upper':str(baseline),'strict_margin_vs_complete_round002_upper':str(answer['inner']-baseline),'positive_inward_raw_boxes':len(positive),'positive_outward_raw_boxes':outer_n,'source_extreme_checks':checks,'minimum_inward_width':str(F(min(width),S)),'global_memo_states':budget['states'],'finite_L1':{'epsilon':str(eps),'eta':str(eta),'theta':str(theta),'K':K,'source_collar':2,'union_shrink_loss_upper':str(loss),'exact_ratio_lower':str(finite),'display':float(finite),'strict_margin_vs_complete_round002_upper':str(finite-baseline)},'input_hashes':{p:sha(HERE/p) for p in ['registration.json','frozen_inputs.json']},'self_sha256':sha(__file__),'pid':os.getpid(),'cpu_seconds':time.process_time()-c0,'wall_seconds':time.time()-t0,'python':sys.executable,'old_seed_evaluations':0,'old_baseline_evaluations':0,'family_exhausted':False,'next_search_authorized':False}
    write('results.json',result)

def main():
    t0=time.time();c0=time.process_time();budget={'states':0}
    assert sha(HERE/'registration.json')=='04fffe13d07393758265ad3223cb25811689730f9e02f476f078338146e96f84'
    frozen=json.loads((HERE/'frozen_inputs.json').read_text())
    assert len(frozen['selected_events'])<=10
    oldkeys={tuple(map(tuple,k)) for k in frozen['prior_evaluated_event_keys']+frozen['excluded_old_seed_keys']}
    assert len({tuple(map(tuple,e['key'])) for e in frozen['selected_events']})==len(frozen['selected_events'])
    assert all(tuple(map(tuple,e['key'])) not in oldkeys for e in frozen['selected_events'])
    roots,tail,records,orbit=tail_build();ivs={};ids={};engines={}
    for mode in ('inner','outer'):
        ivs[mode]=sorted({z for r in orbit for z in r[mode]})
        ids[mode]={z:i for i,z in enumerate(ivs[mode])}
        engines[mode]=Union(ivs[mode],budget)
    ref=F(frozen['complete_R002_interval'][1]);trials=[]
    try:
        for n,event in zip(frozen['selected_event_indices'],frozen['selected_events']):
            assert not event['old_seed_do_not_evaluate']
            low=sum(c*root(m*D**7) for m,c in event['key'])
            high=sum(c*(root(m*D**7)+1) for m,c in event['key'])
            assert low==int(event['a1_lower_numerator']) and high==int(event['a1_upper_numerator'])
            aa,answer=event_compute(event,orbit,engines,ivs,ids)
            trial={'index':n,'event_key':event['key'],'a1_display':event['display'],'full_union_interval':[str(answer['inner']),str(answer['outer'])],'states':budget['states'],'strict_above_complete_round002_upper':answer['inner']>ref}
            trials.append(trial);print(json.dumps(trial),flush=True)
            if answer['inner']>ref:
                certify_winner(event,aa,answer,roots,tail,records,budget,t0,c0)
                break
    except Exhausted:
        write('frontier.json',{'status':'incomplete; not family exclusion','current_event':event,'trials':trials,'global_states':budget['states'],'stacks':{m:e.stack for m,e in engines.items()}})
    write('trials.json',{'trials':trials,'new_event_evaluations':len(trials),'old_seed_evaluations':0,'old_baseline_evaluations':0,'global_states':budget['states'],'stop':'first winner or finite registered limit; no further search'})
    print('registered R003 batch stopped',flush=True)

if __name__=='__main__':main()
