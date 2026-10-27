#!/usr/bin/env python3
from pathlib import Path
from itertools import product
from math import prod
import json,hashlib
HERE=Path(__file__).resolve().parent
D=10**35
def root(n):
    a,b=0,1
    while b**7<=n: b*=2
    while b-a>1:
        c=(a+b)//2
        if c**7<=n:a=c
        else:b=c
    return a
def radical(m):
    a=b=0
    while m%2==0:a+=1;m//=2
    while m%3==0:b+=1;m//=3
    assert m==1
    return 2**(a//7)*3**(b//7),2**(a%7)*3**(b%7)
def key(terms):
    out={}
    for m,c in terms:
        d,t=radical(m);out[t]=out.get(t,0)+c*d
    return tuple(sorted(out.items()))
def main():
    masses={prod(e) for e in product([1,2],*[range(1,5)]*6)}
    rr={m:root(m*D**7) for m in masses}
    u=rr[2]; a=(2*(D+u),2*(D+u+1))
    active=[]
    for e in product([1,2],*[range(1,5)]*6):
        if tuple(sorted(e[1:]))!=e[1:]:continue
        r=(2*rr[prod(e)],2*(rr[prod(e)]+1))
        ok=True
        for ell in e[1:]:
            q=-(ell//2) if ell%2 else 1-ell//2
            lo=max(0,(q+ell-1)*a[0]-r[1]);hi=min(a[1]//2,q*a[0]+r[1])
            if hi<=lo:ok=False;break
        if ok:active.append(e)
    ms=sorted({prod(e) for e in active}); ev={}
    for m in ms:
        ev[key([(m,1)])]=['r',m]
        ev[key([(m,2)])]=['2r',m]
    for i,m in enumerate(ms):
        for n in ms[i:]:ev[key([(m,1),(n,1)])]=['r+r',m,n]
    ev[((1,3),)]=['endpoint','3/2'];ev[((1,4),)]=['endpoint','2']
    kept=[]
    for k,form in ev.items():
        low=sum(c*root(m*D**7) for m,c in k)
        high=sum(c*(root(m*D**7)+1) for m,c in k)
        # key encodes 2*a1; retain any event enclosure intersecting domain.
        if high<3*D or low>4*D:continue
        old=k in (key([(1,2),(16,1)]),key([(1,2),(32,1)]))
        kept.append({'key':k,'formula':form,'a1_lower_numerator':str(low),'a1_upper_numerator':str(high),'a1_denominator':str(2*D),'display':float(low/(2*D)),'old_seed_do_not_evaluate':old})
    kept.sort(key=lambda x:int(x['a1_lower_numerator']))
    old5=next(i for i,e in enumerate(kept) if tuple(map(tuple,e['key']))==key([(1,2),(32,1)]))
    order=sorted((i for i,e in enumerate(kept) if not e['old_seed_do_not_evaluate']),key=lambda i:(abs(i-old5),-i))
    out={'status':'finite event census only; no union evaluations','tail_positive_orbits':len(active),'tail_positive_masses':ms,'events_in_domain':len(kept),'new_events':len(order),'probe_priority_event_indices':order[:512],'all_new_events_within_registered_limit':len(order)<=512,'events':kept,'self_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    with (HERE/'event_census.json').open('x') as f:json.dump(out,f,indent=2);f.write('\n')
    print(json.dumps({k:v for k,v in out.items() if k!='events'},indent=2))
if __name__=='__main__':main()
