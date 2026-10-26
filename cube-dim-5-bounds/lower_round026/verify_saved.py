#!/usr/bin/env python3
"""Audit full raw endpoint labels, affine responses and signed chamber, no unions."""
from pathlib import Path
from fractions import Fraction as Q
from functools import lru_cache
from itertools import permutations,combinations
from collections import Counter
import sys,os,time
sys.dont_write_bytecode=True
P=Path(__file__).resolve().parent;ROOT=P.parent;OLD=ROOT/'lower_round025'
sys.path[:0]=[str(ROOT/'lower'),str(P)]
from certify_rectangular import field
from position_fan import ENV,sha,load,save,enc,rt
def main():
    start=time.monotonic();reg=load(P/'registration.json');assert sha(P/'position_fan.py')==reg['source_sha256']
    for x in reg['inputs']:assert sha(Path(x['path']))==x['sha256']
    for x in load(OLD/'archive_whitelist.json')['files']:assert sha(OLD/x['path'])==x['sha256']
    F,u,v,bounds,_,_=field();a=1+u**3/2;zero=F(0)
    def dec(d):
        z=[Q(0)]*25
        for ij,c in d.items():i,j=map(int,ij.split(','));z[5*i+j]=Q(c)
        return F(a=z)
    @lru_cache(None)
    def radius(n):
        i=j=0
        while n%2==0:n//=2;i+=1
        while n%3==0:n//=3;j+=1
        assert n==1
        return u**i*v**j/2
    geo=load(OLD/'tail_geometry.json');ev=load(P/'events.json');raw=ev['all_raw_inventory'];assert len(raw)==6144
    labels_checked=0
    for index,row in enumerate(geo['raw_inventory']):
        ell=row['lengths'][0];r=radius(row['N'])
        for offset,j in enumerate(range(-1,5)):
            x=raw[6*index+offset];assert x['family_id']==row['family_id'] and x['start']==j
            low=(j+ell-1)*a-r;high=j*a+r
            lo={3} if low<zero else {(j+ell-1)%3} if zero<low else {3,(j+ell-1)%3}
            hi={j%3} if high<3*a else {3} if 3*a<high else {3,j%3}
            assert x['lower_velocity_labels']==sorted(lo) and x['upper_velocity_labels']==sorted(hi);labels_checked+=2
    forms=((Q(1,2),Q(-1,2)),(Q(0),Q(1)),(Q(-1,2),Q(-1,2)),(Q(0),Q(0)))
    global_possible=set()
    for (A,B),(C,D) in combinations(forms,2):
        if B!=D:
            t=(C-A)/(B-D)
            if 0<=t<=Q(1,3):global_possible.add(t)
    assert global_possible=={Q(0),Q(1,3)}
    assert all(Q(x['t']) in global_possible for x in ev['events']) and ev['internal_events']==[]
    hist={int(m):dec(c) for m,c in load(OLD/'tail_histogram.json')['histogram'].items()};checks=[]
    for name in ('plus','minus'):
        r=load(P/f'{name}_affine.json');ck=load(P/f'{name}_checkpoint.json');assert ck['status']=='COMPLETE_all_affine_sections' and ck['completed_sections']==ck['total_sections']==2224
        A=B=F(0);dist=Counter();assert len(ck['section_affine_derivatives'])==len(hist)
        for x in ck['section_affine_derivatives']:
            m=int(x['family_mask']);aa,bb=Q(x['constant']),Q(x['slope']);A+=hist[m]*aa;B+=hist[m]*bb;dist[(str(aa),str(bb))]+=1
        assert A==dec(ck['partial_A']) and B==dec(ck['partial_B'])
        A=A*Q(16,3);B=B*Q(16,3)
        assert A==dec(r['D_constant']) and B==dec(r['D_slope'])==zero
        oldD=dec(load(OLD/f'{name}_derivative.json')['directional_derivative_coefficients']);assert A==oldD/2
        assert list(map(str,bounds(A.a)))==r['sample_derivative_interval'] and bounds(A.a)[1]<0
        checks.append({'profile':name,'section_A_B_counts':[{'constant':k[0],'slope':k[1],'count':n} for k,n in sorted(dist.items())],'constant_coefficients':enc(A),'slope_exactly_zero':True})
    # Full range1 hexagon boundary consists of ordinary permutations of v(t).
    # Even permutations are G*v; odd permutations are G*(-v). Each edge can
    # be halved because -P02 maps v(t) to v(-t), also for -v.
    vf=forms[:3]
    def parity(p):return -1 if sum(p[i]>p[j] for i in range(3) for j in range(i+1,3))%2 else 1
    for p in permutations(range(3)):
        ordinary=tuple(vf[p[i]] for i in range(3));profile=vf if parity(p)==1 else tuple((-A,-B) for A,B in vf)
        signed=tuple(tuple(parity(p)*x for x in profile[p[i]]) for i in range(3));assert ordinary==signed
    flipped=tuple(tuple(-x for x in vf[i]) for i in (2,1,0));assert flipped==tuple((A,-B) for A,B in vf)
    endpoint=tuple(A+B*Q(1,3) for A,B in vf);assert endpoint[0]==endpoint[1]
    negative=tuple(-x for x in endpoint);assert tuple(-endpoint[i] for i in (1,0,2))==negative
    final=load(P/'fan_results.json');assert final['status']=='STRICT_NEGATIVE_full_position_fan' and final['bothhalfedgesconstant']
    assert bounds(dec(final['full_fan_maximum_coefficients']).a)[1]<0
    out={'status':'PASS_full_signedD3_position_fan_saved_certificate','raw_endpoint_label_formulas_checked':labels_checked,'raw_objects_checked':6144,'global_velocity_equalities_have_no_open_internal_roots':True,'actual_tie_event_count':len(ev['events']),'section_profiles':checks,'signedG_hexagon_coverage_and_halfedge_reflection_checked':True,'closed_t1_3_signed_equivalence_checked':True,'old025_archive_files_unchanged':True,'new_union_or_histogram_jobs':0,'runtime':rt(start),'source_sha256':sha(Path(__file__)),'fan_result_sha256':sha(P/'fan_results.json')}
    save(P/'verification.json',out);print({k:v for k,v in out.items() if k!='section_profiles'})
if __name__=='__main__':
    assert all(os.environ.get(k)=='1' for k in ENV) and os.getpriority(os.PRIO_PROCESS,0)>=10
    main()
