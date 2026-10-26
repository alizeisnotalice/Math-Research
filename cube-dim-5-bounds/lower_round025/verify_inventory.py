#!/usr/bin/env python3
"""Independent full1024 scalar tail/first-span qualification; no integration."""
from pathlib import Path
from itertools import product
from math import prod
from fractions import Fraction as Q
from functools import lru_cache
import sys,os,time
sys.dont_write_bytecode=True
P=Path(__file__).resolve().parent;ROOT=P.parent
sys.path[:0]=[str(ROOT/'lower'),str(P)]
from certify_rectangular import field
from position_sections import ENV,load,save,sha,rt,enc
def main():
    start=time.monotonic();F,u,v,_,_,_=field();a=1+u**3/2;s=(1+u)/2;zero=F(0)
    @lru_cache(None)
    def radius(n):
        x=n;i=j=0
        while x%2==0:x//=2;i+=1
        while x%3==0:x//=3;j+=1
        assert x==1
        r=u**i*v**j/2;assert 32*r**5==F(n) and r<=F(2);return r
    g=load(P/'tail_geometry.json');fams={x['family_id']:x for x in g['families']};projection=load(P/'projection_certificate.json')
    assert len(g['raw_inventory'])==1024
    positive=[];guards=0
    for fid,(ls,row) in enumerate(zip(product(range(1,5),repeat=5),g['raw_inventory'])):
        assert row['family_id']==fid and list(ls)==row['lengths'] and row['N']==prod(ls)
        r=radius(prod(ls));box=[]
        for ell in ls[1:]:
            j=-(ell-1)//2 if ell%2 else 1-ell//2
            box.extend((max(zero,(j+ell-1)*s-r),min(s/2,j*s+r)));guards+=2
        valid=all(box[2*k]<box[2*k+1] for k in range(4))
        assert valid==(row['status']=='positive_tail_box')==(fid in fams)
        if valid:
            positive.append(fid);assert fams[fid]['radius']==enc(r) and fams[fid]['tail_box']==list(map(enc,box))
    keep=set(projection['kept_first_span_families']);excluded={x['family_id'] for x in projection['excluded_strict_span_families']}
    assert keep.isdisjoint(excluded) and keep|excluded==set(positive)
    for fid in positive:
        row=g['raw_inventory'][fid];span=(row['lengths'][0]-1)*a;diameter=2*radius(row['N'])
        assert (span<=diameter)==(fid in keep)
    out={'status':'PASS_full1024_independent_scalar_qualification','raw_length_vectors':1024,'tail_endpoint_guards':guards,'positive_tail_boxes':len(positive),'local_firstspan_eligible':len(keep),'strict_firstspan_exclusions':len(excluded),'zero_width_firstspan_retained_by_nonstrict_test':True,'no_union_or_section_recomputation':True,'runtime':rt(start),'source_sha256':sha(Path(__file__)),'geometry_sha256':sha(P/'tail_geometry.json')}
    save(P/'inventory_verification.json',out);print(out)
if __name__=='__main__':
    assert all(os.environ.get(k)=='1' for k in ENV) and os.getpriority(os.PRIO_PROCESS,0)>=10
    main()
