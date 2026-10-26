#!/usr/bin/env python3
"""Only equal-base raw endpoint speed comparisons; no union evaluations."""
from pathlib import Path
from fractions import Fraction as Q
from functools import lru_cache,cmp_to_key
from itertools import product,combinations
from math import prod
from collections import defaultdict,Counter
import sys,os,time,json,hashlib
sys.dont_write_bytecode=True
P=Path(__file__).resolve().parent;ROOT=P.parent
sys.path[:0]=[str(ROOT/'lower'),str(P)]
from certify_rectangular import field
from profile_variation import ENV,sha,load,save,rt,enc

def main():
    start=time.monotonic();reg=load(P/'event_registration.json')
    assert sha(Path(__file__))==reg['source_sha256']
    for row in reg['inputs']:assert sha(Path(row['path']))==row['sha256']
    F,u,v,bounds,roots,den=field();zero=F(0);a=(1+u**3/2,)+((1+u)/2,)*4
    @lru_cache(None)
    def radius(n):
        t=n;e2=e3=0
        while t%2==0:t//=2;e2+=1
        while t%3==0:t//=3;e3+=1
        assert t==1
        return u**e2*v**e3/2
    groups=[defaultdict(set) for _ in range(5)];counts=Counter();records=[]
    pattern=((Q(1,2),Q(-1,2)),(Q(0),Q(1)),(Q(-1,2),Q(-1,2)))
    for ls in product(range(1,5),repeat=5):
        n=prod(ls);r=radius(n);tail=prod(ls[1:])
        for j in range(-1,5):
            counts['raw']+=1;starts=[j]+[-(ell-1)//2 if ell%2 else 1-ell//2 for ell in ls[1:]]
            dmA=sum(pattern[(j+t)%3][0] for t in range(ls[0]))*tail
            dmB=sum(pattern[(j+t)%3][1] for t in range(ls[0]))*tail
            dA=r*dmA/(5*n);dB=r*dmB/(5*n);axis=[]
            for k,(ell,st) in enumerate(zip(ls,starts)):
                domain=3*a[k] if k==0 else a[k]/2
                rawlo=(st+ell-1)*a[k]-r;rawhi=st*a[k]+r
                lo=max(zero,rawlo);hi=min(domain,rawhi)
                los={(zero,zero)} if rawlo<zero else {(-dA,-dB)} if zero<rawlo else {(zero,zero),(-dA,-dB)}
                his={(dA,dB)} if rawhi<domain else {(zero,zero)} if domain<rawhi else {(zero,zero),(dA,dB)}
                axis.append((lo,hi,los,his))
            if any(hi<lo for lo,hi,_,_ in axis):counts['strict_negative_base']+=1;continue
            if any(lo==hi for lo,hi,_,_ in axis):counts['zero_base_raw']+=1
            else:counts['positive_base_raw']+=1
            for k,(lo,hi,los,his) in enumerate(axis):
                groups[k][lo].update(los);groups[k][hi].update(his)
            records.append({'lengths':ls,'start':j,'mass_slope_affine':[str(dmA),str(dmB)]})
    assert counts['raw']==6144 and counts['positive_base_raw']==179 and counts['zero_base_raw']==0
    paircount=0;events=[];groupcounts=[]
    for k,gg in enumerate(groups):
        groupcounts.append({'coordinate':k,'base_positions':len(gg),'velocity_form_counts':sorted(len(v) for v in gg.values())})
        for base,forms in gg.items():
            for (A,B),(C,D) in combinations(sorted(forms),2):
                paircount+=1
                if paircount>200000:
                    save(P/'event_checkpoint.json',{'status':'LIMIT_paircomparison','completed_comparisons':paircount-1,'counts':dict(counts),'runtime':rt(start)});return
                if B==D:continue
                numerator=C-A;denominator=B-D
                if denominator<zero:numerator=-numerator;denominator=-denominator
                if numerator<denominator*Q(-1,3) or denominator*Q(1,3)<numerator:continue
                events.append((numerator,denominator,{'coordinate':k,'base':enc(base),'left_velocity':[enc(A),enc(B)],'right_velocity':[enc(C),enc(D)]}))
    def cmp(x,y):
        diff=x[0]*y[1]-y[0]*x[1]
        return 0 if diff==zero else (-1 if diff<zero else 1)
    events.sort(key=cmp_to_key(cmp));unique=[]
    for x in events:
        if unique and cmp(x,unique[-1])==0:continue
        unique.append(x)
    def event_row(x):
        num,dn,evidence=x;nl,nh=bounds(num.a);dl,dh=bounds(dn.a);assert dl>0
        options=[nl/dl,nl/dh,nh/dl,nh/dh];lo,hi=min(options),max(options)
        rational=None
        for v,w in zip(num.a,dn.a):
            if w:
                q=v/w
                if num==dn*q:rational=str(q)
                break
        return {'numerator':enc(num),'positive_denominator':enc(dn),'rational_value_if_exact':rational,'isolating_interval':[str(lo),str(hi)],'comparison_witness':evidence}
    out={'status':'COMPLETE_all_raw_equalbase_velocity_events','counts':dict(counts),'raw_endpoint_pair_comparisons':paircount,'event_occurrences':len(events),'distinct_breakpoint_candidates_in_closed_window':len(unique),'events':[event_row(x) for x in unique],'groups':groupcounts,'full_active_raw_affine_mass_inventory':records,'mechanism':'all6144raw beforefirstvariationpruning;allclippingalternativesatequalbase;onlyequalbasevelocitycomparisonsmatter;candidateeventsetmayincludeinactiveorunionirrelevantcrossingsbutcannotomittruefirstorderchanges','runtime':rt(start),'source_sha256':reg['source_sha256'],'registration_sha256':sha(P/'event_registration.json')}
    save(P/'events.json',out);print(json.dumps({k:v for k,v in out.items() if k not in ('events','groups','full_active_raw_affine_mass_inventory')},indent=2))
    print('event rational display only:',[x['rational_value_if_exact'] for x in out['events']],flush=True)
if __name__=='__main__':
    assert all(os.environ.get(k)=='1' for k in ENV) and os.getpriority(os.PRIO_PROCESS,0)>=10
    main()
