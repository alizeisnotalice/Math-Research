#!/usr/bin/env python3
"""Registered n=7 finite lattice witness. Integer outward arithmetic only.

All coordinates have denominator 4D. No old baseline is evaluated.
Tail boxes are complete permutation orbits; recursive sweep respects orbits.
"""
from pathlib import Path
from fractions import Fraction as F
from itertools import product
from math import prod, factorial
import sys, json, hashlib, time, os, platform

HERE=Path(__file__).resolve().parent
D=10**35
S=4*D
LIMIT=2**20

def sha(path): return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def floor_root(n,k=7):
    lo,hi=0,1
    while hi**k<=n: hi*=2
    while hi-lo>1:
        mid=(lo+hi)//2
        if mid**k<=n: lo=mid
        else: hi=mid
    assert lo**k<=n<(lo+1)**k
    return lo

def enc(x): return str(x)
def dec(x): return format(float(x),'.17g')
def jwrite(p,obj):
    with Path(p).open('x') as f: json.dump(obj,f,indent=2); f.write('\n')

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

def build(power):
    roots={m:floor_root(m*D**7) for m in {prod(e) for e in product([1,2],*[range(1,5)]*6)}}
    ru=roots[2]; rp=floor_root(2**power*D**7)
    # Each pair is a rational enclosure in common units 1/(4D).
    tail=(2*(D+ru),2*(D+ru+1))
    first=(4*D+2*rp,4*D+2*(rp+1))
    widths=[]; records=[]; contain_checks=0; positive_raw=0
    for e in product([1,2],*[range(1,5)]*6):
        m=prod(e); rr=(2*roots[m],2*(roots[m]+1)); boxes=[]
        for i,ell in enumerate(e):
            aa=first if i==0 else tail
            q=-(ell//2) if ell%2 else 1-ell//2
            # (q+ell-1)>=0, q<=0; explicit outward interval arithmetic.
            c=q+ell-1
            ll=(max(0,c*aa[0]-rr[1]),max(0,c*aa[1]-rr[0]))
            hh=(min(aa[0]//2,q*aa[1]+rr[0]),min(aa[1]//2,q*aa[0]+rr[1]))
            assert aa[0]%2==0 and aa[1]%2==0
            boxes.append((ll,hh))
        inner=tuple((l[1],h[0]) for l,h in boxes)
        outer=tuple((l[0],h[1]) for l,h in boxes)
        if all(l<h for l,h in outer):
            positive_raw+=1
            assert all(l<h for l,h in inner), 'undecided positive width'
            for i,(ell,(l,h)) in enumerate(zip(e,inner)):
                aa=first if i==0 else tail
                q=-(ell//2) if ell%2 else 1-ell//2
                # Every actual source atom lies in [q*a,(q+ell-1)*a].
                assert l >= (q+ell-1)*aa[1]-rr[0]
                assert h <= q*aa[1]+rr[0]
                contain_checks+=2
                widths.append(h-l)
            records.append({'ell':e,'mass':m,'inner':inner,'outer':outer})
    # Full list has all 8192 vectors; quotient only after complete checks.
    orbit=[]
    for e in product([1,2],*[range(1,5)]*6):
        if list(e[1:])!=sorted(e[1:]): continue
        r=next((r for r in records if r['ell']==e),None)
        if r: orbit.append(r)
    assert sum(factorial(6)//prod(factorial(e['ell'][1:].count(x)) for x in set(e['ell'][1:])) for e in orbit)==positive_raw
    return roots,first,tail,records,orbit,contain_checks,min(widths)

def compute(power):
    t0=time.time(); c0=time.process_time()
    roots,first,tail,records,orbit,checks,minwidth=build(power)
    budget={'states':0}; results={}; engines=[]
    for direction in ('inner','outer'):
        iv=sorted({i for r in orbit for i in r[direction][1:]})
        ids={x:i for i,x in enumerate(iv)}
        boxes=[(r[direction][0],tuple(sorted(ids[z] for z in r[direction][1:]))) for r in orbit]
        cuts=sorted({x for b,_ in boxes for x in b})
        engine=Union(iv,budget); engines.append(engine)
        total=0
        for a,b in zip(cuts,cuts[1:]):
            state=tuple(sorted({tailbox for (lo,hi),tailbox in boxes if lo<=a and b<=hi}))
            try: total+=(b-a)*engine.call(state)
            except Exhausted:
                jwrite(HERE/f'seed{power}_frontier.json',{'direction':direction,'completed_first_volume':str(total),'first_cuts':list(map(str,cuts)),'current_slab':[str(a),str(b)],'pending_state':state,'recursive_stack':engine.stack,'intervals':[[str(l),str(h)] for l,h in iv],'states':budget['states'],'status':'incomplete; not an exclusion'})
                raise
        results[direction]=F(128*total,S**7)
        print(json.dumps({'seed':power,'direction':direction,'states':budget['states'],'value_display':dec(results[direction]),'wall':time.time()-t0}),flush=True)
    ref=json.loads((HERE/'frozen_internal_reference.json').read_text())
    reflo,refhi=map(F,ref['interval'])
    assert results['inner']<=results['outer']
    # Shrink each of the <=8192 raw boxes by epsilon at every face.
    eps=F(1,10**25); eta=eps/2; theta=1-F(1,10**18); K=10**19
    assert F(minwidth,S)>2*eps
    # Every side is <= max(a1/2,a/2)<2. Thus face-slab bound 14 eps 2^6.
    assert max(first[1],tail[1])<4*S
    loss=128*len(records)*14*eps*2**6
    shrinklb=results['inner']-loss
    finite=theta*shrinklb*F(2*K+1,2*K+5)**7
    out={'status':'complete strict rational full union certificate','dimension':7,'seed_first_power':power,'raw_length_vectors':8192,'positive_raw_boxes':len(records),'symmetric_box_orbits':len(orbit),'source_extreme_containment_checks':checks,'deduplicated_states':budget['states'],'state_limit':LIMIT,'root_denominator':str(D),'root_lower_numerators':{str(m):str(v) for m,v in roots.items()},'first_root_lower_numerator':str(floor_root(2**power*D**7)),'coordinate_scale':str(S),'union_lower':str(results['inner']),'union_upper':str(results['outer']),'union_display':[dec(results['inner']),dec(results['outer'])],'complete_R075_upper':str(refhi),'strict_above_complete_R075_upper':results['inner']>refhi,'exact_margin_vs_R075_upper':str(results['inner']-refhi),'finite_L1':{'epsilon':str(eps),'eta':str(eta),'theta':str(theta),'K':K,'source_collar':2,'raw_box_count_for_loss':len(records),'union_shrink_loss_upper':str(loss),'finite_ratio_lower':str(finite),'finite_ratio_display':dec(finite),'finite_ratio_above_R075_upper':finite>refhi},'min_inward_width':str(F(minwidth,S)),'wall_seconds':time.time()-t0,'cpu_seconds':time.process_time()-c0,'pid':os.getpid(),'python':sys.executable,'platform':platform.platform(),'registration_sha256':sha(HERE/'registration.json'),'reference_sha256':sha(HERE/'frozen_internal_reference.json'),'self_sha256':sha(__file__)}
    jwrite(HERE/f'seed{power}_results.json',out)
    return out

def main():
    assert sha(HERE/'registration.json')=='7e799455eab9787fbddbaeca131a69f0c8649de49d9a3ecfd18ccfaff7dc9a16'
    for p in (4,5):
        try: compute(p)
        except Exhausted:
            print(json.dumps({'seed':p,'status':'state limit exhausted; frontier saved; not excluded'}),flush=True)
    print('registered finite batch stopped',flush=True)

if __name__=='__main__': main()
