#!/usr/bin/env python3
"""Check saved sections and affine certificates; no new joint-mask pairs."""
from pathlib import Path
from fractions import Fraction as Q
import json,hashlib,os,time,signal
from sections import G,ENV,rt,dec
P=Path(__file__).resolve().parent
def load(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def area(b):return max(0,b[1]-b[0])*max(0,b[3]-b[2])
def direct_length(desc,a):
    # Generic clipped interval union, independent of the anchored-length formula.
    intervals=[]
    for ell,r in desc:
        ii=range(-(ell-1)//2,(ell-1)//2+1) if ell%2 else range(1-ell//2,1+ell//2)
        lo=max(Q(0),max(ii)*a-r);hi=min(a/2,min(ii)*a+r)
        if lo<hi:intervals.append((lo,hi))
    intervals.sort();answer=Q(0);last=None
    for lo,hi in intervals:
        if last is None:last=[lo,hi]
        elif lo<=last[1]:last[1]=max(last[1],hi)
        else:answer+=last[1]-last[0];last=[lo,hi]
    return answer+(last[1]-last[0] if last else Q(0))
def main():
    start=time.monotonic();signal.alarm(180)
    reg=load(P/'registration.json');r=load(P/'results.json')
    assert sha(P/'sections.py')==reg['source_sha256']
    for x in reg['inputs']:assert sha(P.parent/x['path'])==x['sha256']
    # Every previously archived021 file remains byte-identical.
    old=P.parent/'lower_round021';manifest=load(old/'archive_whitelist.json')
    for x in manifest['files']:assert sha(old/x['path'])==x['sha256']
    tail=load(P/'tail_histogram.json');joint=load(P/'section_histogram.json');ck=load(P/'checkpoint.json')
    th={int(k):v for k,v in tail['histogram'].items()};jh={int(k):v for k,v in joint['histogram'].items()}
    assert len(tail['shapes'])==54 and sum(th.values())==G**2//4 and sum(jh.values())==G**4//4
    assert ck['status']=='COMPLETE_fourD_shape_sections' and ck['completed_pairs']==ck['total_pairs']==110302
    assert ck['joint_histogram']==joint['histogram'] and ck['next_pair_indices']==[262,0]
    plane=load(old/'plane_histogram.json');ph={int(k):v for k,v in plane['histogram'].items()}
    marginal_checks=0
    for s in tail['shapes']:
        bit=1<<s['shape_id'];tb=area(s['tail_inner_box'])
        assert sum(v for m,v in th.items() if m&bit)==tb
        planar=sum(v for m,v in ph.items() if m&(1<<s['family_id']))
        assert sum(v for m,v in jh.items() if m&bit)==tb*planar
        marginal_checks+=1
    descriptor_weights={}
    for mask,w in jh.items():
        best={}
        for s in tail['shapes']:
            if mask&(1<<s['shape_id']):
                ell=s['shape'][0];best[ell]=max(best.get(ell,Q(0)),Q(s['r_lower']))
        key=tuple(sorted(best.items()));descriptor_weights[key]=descriptor_weights.get(key,0)+w
    raw_bound=Q(0);pieces_checked=0;point_checks=0
    for row in r['descriptors']:
        desc=tuple((ell,Q(v)) for ell,v in row['descriptor']);weight=Q(row['weight'])
        assert weight==Q(descriptor_weights.pop(desc),G**4)
        segments=row['segments'];assert Q(segments[0]['left'])==Q(3,2) and Q(segments[-1]['right'])==2
        maximum=Q(0);previous=None
        for seg in segments:
            a,b=Q(seg['left']),Q(seg['right']);s,t=Q(seg['slope']),Q(seg['intercept'])
            assert a<b and (previous is None or previous==a);previous=b
            for x in (a,(2*a+b)/3,(a+b)/2,(a+2*b)/3,b):
                assert direct_length(desc,x)==s*x+t;point_checks+=1
            maximum=max(maximum,s*a+t,s*b+t);pieces_checked+=1
        assert maximum==Q(row['independent_max_length']);raw_bound+=weight*maximum
    assert not descriptor_weights
    geom=load(old/'candidate_geometry.json');sh=Q(geom['s_interval'][1])
    assert 8*sh**4*raw_bound==Q(r['frozen_sectionwise_independent_max_upper'])
    err=sum((Q(x['uniform_product_error_upper']) for x in r['shell_rows']),Q(0))
    assert err==Q(r['uniform_nominal_error_before_8s4'])
    upper=8*sh**4*(raw_bound+err);assert upper==Q(r['uniform_nominal_finitefunctional_upper'])
    C=load(P.parent/'lower_round003/a2_equals_a3_results.json');cl=Q(C['ratio_interval'][0]);assert cl-upper>0
    out={'status':'PASS_saved_section_and_affine_interface_audit','all021_whitelist_files_unchanged':True,'new_joint_mask_pairs':0,'tail_and_joint_single_shape_marginals_verified':marginal_checks,'descriptor_weights_verified':len(r['descriptors']),'affine_segments_verified_against_generic_interval_union':pieces_checked,'generic_interval_union_checks':point_checks,'section_weight_total':'1/4','strict_C_gap_positive':str(cl-upper),'window_nominal_upper_decimal_display':dec(upper),'runtime':rt(start),'source_sha256':sha(Path(__file__)),'result_sha256':sha(P/'results.json')}
    (P/'verification.json').write_text(json.dumps(out,separators=(',',':'))+'\n');signal.alarm(0);print(json.dumps(out,indent=2))
if __name__=='__main__':
    assert all(os.environ.get(k)=='1' for k in ENV) and os.getpriority(os.PRIO_PROCESS,0)>=10
    main()
