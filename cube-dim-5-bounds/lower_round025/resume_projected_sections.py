#!/usr/bin/env python3
"""Project a successful prefix; integrate only its unfinished region on safe coarse cuts."""
from pathlib import Path
from fractions import Fraction as Q
from itertools import product
from math import prod
from collections import defaultdict
import sys,os,time
sys.dont_write_bytecode=True
P=Path(__file__).resolve().parent;ROOT=P.parent
sys.path[:0]=[str(ROOT/'lower'),str(P)]
from certify_rectangular import field
from position_sections import ENV,load,save,sha,rt,enc
def main():
    start=time.monotonic();reg=load(P/'repair_registration.json')
    assert sha(Path(__file__))==reg['source_sha256']
    for x in reg['inputs']:assert sha(Path(x['path']))==x['sha256']
    assert not (P/'tail_histogram.json').exists(),'completehistogram alreadyexists'
    F,u,v,bounds,_,_=field();zero=F(0);a=1+u**3/2;s=(1+u)/2
    def dec(d):
        z=[Q(0)]*25
        for ij,c in d.items():i,j=map(int,ij.split(','));z[5*i+j]=Q(c)
        return F(a=z)
    geom=load(P/'tail_geometry.json');old=load(P/'tail_checkpoint.json');completed=old['completed_first_tail_slabs'];cut=dec(geom['axes'][0][completed])
    kept=[];excluded=[];min_gap_lower=None
    for row in geom['families']:
        r=dec(row['radius']);span=(row['lengths'][0]-1)*a
        if 2*r<span:
            gap=span-2*r;lower=bounds(gap.a)[0];assert lower>0
            min_gap_lower=lower if min_gap_lower is None else min(min_gap_lower,lower)
            excluded.append({'family_id':row['family_id'],'strict_negative_first_axis_width':enc(-gap)})
        else:kept.append(dict(row,box=tuple(dec(z) for z in row['tail_box'])))
    keepmask=sum(1<<r['family_id'] for r in kept)
    # Coefficient-wise aggregation avoids caching tens of thousands of useless
    # fine masks. Every already-computed geometric cell is reused exactly once.
    arrays=defaultdict(lambda:[Q(0)]*25)
    for m,coeff in old['histogram'].items():
        target=arrays[int(m)&keepmask]
        for ij,c in coeff.items():i,j=map(int,ij.split(','));target[5*i+j]+=Q(c)
    hist={m:F(a=z) for m,z in arrays.items()};prefix_mass=sum(hist.values(),F(0))
    assert prefix_mass==cut*(s/2)**3
    axes=[]
    for k in range(4):
        left=cut if k==0 else zero
        axes.append(sorted({left,s/2}|{x for r in kept for x in r['box'][2*k:2*k+2] if left<x<s/2}))
    segments=[]
    for k,axis in enumerate(axes):
        segments.append([(hi-lo,sum(1<<r['family_id'] for r in kept if r['box'][2*k]<=lo and hi<=r['box'][2*k+1])) for lo,hi in zip(axis,axis[1:])])
    save(P/'projection_certificate.json',{'status':'EXACT_reuse_of_successful_fine_prefix','original_completed_slabs':completed,'original_completed_cells':old['nodes_this_process'],'prefix_upper_first_tail_coordinate':enc(cut),'projected_prefix_section_count':len(hist),'kept_first_span_families':[r['family_id'] for r in kept],'excluded_strict_span_families':excluded,'valid_epsilon_upper':str(min(Q(1,8),min_gap_lower/4)),'exclusion_proof':'firstaxisrawunclippedwidth=2r-(ell1-1)a+epsilon*(zfirst-zlast);abs(zfirst-zlast)<=2;sostrictnegativewidthcannotopenwhen epsilon<minimum_gap/4;zero_widthfamiliesnotexcluded','remaining_axis_segments':list(map(len,segments)),'remaining_cells':prod(map(len,segments)),'old_checkpoint_sha256':sha(P/'tail_checkpoint.json'),'source_sha256':reg['source_sha256'],'runtime':rt(start)})
    cp=P/'projected_checkpoint.json';prefix=0;nodes=0
    if cp.exists():
        ck=load(cp);assert ck['source_sha256']==reg['source_sha256']
        prefix=ck['completed_remaining_slabs'];hist={int(m):dec(z) for m,z in ck['histogram'].items()}
    def checkpoint(status):
        save(cp,{'status':status,'completed_remaining_slabs':prefix,'total_remaining_slabs':len(segments[0]),'new_cells_this_process':nodes,'histogram':{str(m):enc(z) for m,z in sorted(hist.items())},'source_sha256':reg['source_sha256'],'runtime':rt(start)})
    print({'stage':'projected_unfinished_region_only','pid':os.getpid(),'kept':len(kept),'excluded':len(excluded),'remaining_cells':prod(map(len,segments)),'axis_segments':list(map(len,segments))},flush=True)
    for i in range(prefix,len(segments[0])):
        if nodes+prod(map(len,segments[1:]))>200000:checkpoint('FINITE_NODE_LIMIT');return
        dx,mx=segments[0][i]
        for (dy,my),(dz,mz),(dw,mw) in product(*segments[1:]):
            nodes+=1;m=mx&my&mz&mw;hist[m]=hist.get(m,F(0))+dx*dy*dz*dw
        prefix=i+1;checkpoint('EXACT_remaining_region_prefix_saved')
    assert sum(hist.values(),F(0))==(s/2)**4
    result={'status':'COMPLETE_static_four_tail_family_sections','histogram':{str(m):enc(z) for m,z in sorted(hist.items())},'section_count':len(hist),'full_domain_volume':enc((s/2)**4),'positive_tail_families':len(kept),'successfully_reused_fine_cells':old['nodes_this_process'],'new_remaining_coarse_cells':nodes,'old_cells_recomputed':0,'projection_certificate_sha256':sha(P/'projection_certificate.json'),'geometry_sha256':sha(P/'tail_geometry.json'),'source_sha256':reg['source_sha256'],'runtime':rt(start)}
    save(P/'tail_histogram.json',result);checkpoint('COMPLETE_static_projected_sections')
    print({k:v for k,v in result.items() if k not in ('histogram','full_domain_volume')},flush=True)
if __name__=='__main__':
    assert all(os.environ.get(k)=='1' for k in ENV) and os.getpriority(os.PRIO_PROCESS,0)>=10
    main()
