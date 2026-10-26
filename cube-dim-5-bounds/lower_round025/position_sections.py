#!/usr/bin/env python3
"""Period3 position shifts: one static four-tail histogram, two 1D union jets."""
from pathlib import Path
from fractions import Fraction as Q
from functools import lru_cache
from itertools import product
from math import prod
from collections import defaultdict,Counter
import sys,os,time,threading,resource,hashlib,json
sys.dont_write_bytecode=True
P=Path(__file__).resolve().parent;ROOT=P.parent
sys.path.insert(0,str(ROOT/'lower'))
from certify_rectangular import field
ENV=('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','VECLIB_MAXIMUM_THREADS','BLIS_NUM_THREADS')
def load(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def save(p,d):
    tmp=p.with_suffix(p.suffix+'.tmp');tmp.write_text(json.dumps(d,separators=(',',':'))+'\n');tmp.replace(p)
def enc(x):return {f'{i//5},{i%5}':str(c) for i,c in enumerate(x.a) if c}
def rt(start):
    rss=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return {'pid':os.getpid(),'nice':os.getpriority(os.PRIO_PROCESS,0),'threads':threading.active_count(),'peak_rss_bytes':rss*(1 if sys.platform=='darwin' else 1024),'wall_seconds':time.monotonic()-start,'thread_environment':{k:os.environ.get(k) for k in ENV}}

def main(stage):
    start=time.monotonic();reg=load(P/'registration.json')
    assert sha(Path(__file__))==reg['source_sha256']
    for x in reg['inputs']:assert sha(Path(x['path']))==x['sha256']
    F,u,v,bounds,roots,den=field();zero=F(0);a=1+u**3/2;s=(1+u)/2
    def dec(d):
        x=[Q(0)]*25
        for ij,c in d.items():i,j=map(int,ij.split(','));x[5*i+j]=Q(c)
        return F(a=x)
    @lru_cache(None)
    def rad(n):
        t=n;i=j=0
        while t%2==0:t//=2;i+=1
        while t%3==0:t//=3;j+=1
        assert t==1
        r=u**i*v**j/2;assert 32*r**5==F(n);return r
    def build_families():
        families=[];inventory=[]
        for family_id,ls in enumerate(product(range(1,5),repeat=5)):
            n=prod(ls);r=rad(n);box=[]
            for ell in ls[1:]:
                st=-(ell-1)//2 if ell%2 else 1-ell//2
                box.extend((max(zero,(st+ell-1)*s-r),min(s/2,st*s+r)))
            status='positive_tail_box' if all(box[2*k]<box[2*k+1] for k in range(4)) else 'tail_empty_or_fixed_zero_width'
            inventory.append({'family_id':family_id,'lengths':ls,'N':n,'status':status})
            if status=='positive_tail_box':families.append({'family_id':family_id,'lengths':ls,'N':n,'radius':r,'tail_box':tuple(box)})
        assert len(inventory)==1024
        return families,inventory
    if stage=='tail':
        assert not (P/'tail_histogram.json').exists(),'complete static histogram already exists;do not rerun'
        families,inventory=build_families();axes=[sorted({zero,s/2}|{f['tail_box'][2*k+j] for f in families for j in (0,1)}) for k in range(4)]
        segments=[]
        for k,axis in enumerate(axes):
            row=[]
            for lo,hi in zip(axis,axis[1:]):
                mask=sum(1<<f['family_id'] for f in families if f['tail_box'][2*k]<=lo and hi<=f['tail_box'][2*k+1])
                row.append((hi-lo,mask))
            segments.append(row)
        save(P/'tail_geometry.json',{'status':'COMPLETE_full1024_tail_geometry','raw_inventory':inventory,'families':[{'family_id':f['family_id'],'lengths':f['lengths'],'N':f['N'],'radius':enc(f['radius']),'tail_box':[enc(x) for x in f['tail_box']]} for f in families],'axes':[[enc(x) for x in axis] for axis in axes],'source_sha256':reg['source_sha256'],'runtime':rt(start)})
        hist=defaultdict(lambda:F(0));prefix=0;nodes=0;cp=P/'tail_checkpoint.json'
        if cp.exists():
            ck=load(cp);assert ck['source_sha256']==reg['source_sha256']
            prefix=ck['completed_first_tail_slabs'];hist.update({int(k):dec(x) for k,x in ck['histogram'].items()})
        def checkpoint(status):
            save(cp,{'status':status,'completed_first_tail_slabs':prefix,'total_first_tail_slabs':len(segments[0]),'nodes_this_process':nodes,'histogram':{str(k):enc(x) for k,x in sorted(hist.items())},'source_sha256':reg['source_sha256'],'runtime':rt(start)})
        print(json.dumps({'stage':'static_four_tail_sections','positive_tail_families':len(families),'axis_segments':list(map(len,segments)),'total_cells':prod(map(len,segments)),'runtime':rt(start)}),flush=True)
        for ix in range(prefix,len(segments[0])):
            if nodes+prod(map(len,segments[1:]))>200000:checkpoint('FINITE_NODE_LIMIT_no_complete_histogram');return
            dx,mx=segments[0][ix]
            for (dy,my),(dz,mz),(dw,mw) in product(*segments[1:]):
                nodes+=1;hist[mx&my&mz&mw]+=dx*dy*dz*dw
            prefix=ix+1;checkpoint('EXACT_first_tail_slab_prefix_saved')
        total=sum(hist.values(),F(0));assert total==(s/2)**4
        result={'status':'COMPLETE_static_four_tail_family_sections','histogram':{str(k):enc(x) for k,x in sorted(hist.items())},'section_count':len(hist),'full_domain_volume':enc(total),'positive_tail_families':len(families),'cells':nodes,'geometry_sha256':sha(P/'tail_geometry.json'),'source_sha256':reg['source_sha256'],'runtime':rt(start)}
        save(P/'tail_histogram.json',result);checkpoint('COMPLETE_static_four_tail_sections')
        print(json.dumps({k:v for k,v in result.items() if k not in ('histogram','full_domain_volume')},indent=2),flush=True);return
    assert stage in reg['profiles']
    assert not (P/f'{stage}_derivative.json').exists(),'complete derivative exists;do not rerun'
    pattern=tuple(map(Q,reg['profiles'][stage]));geom=load(P/'tail_geometry.json');hh=load(P/'tail_histogram.json');assert hh['geometry_sha256']==sha(P/'tail_geometry.json')
    families={f['family_id']:f for f in geom['families']};hist=sorted((int(k),dec(x)) for k,x in hh['histogram'].items())
    # Jet is (exact algebraic base coordinate, rational speed), ordered epsilon->0+.
    interval_by_family={};counts=Counter();zero_openings=[];ties=defaultdict(set)
    for f in geom['raw_inventory']:
        fid=f['family_id'];ell=f['lengths'][0];r=rad(f['N']);intervals=[]
        for j in range(-1,5):
            counts['raw_candidates']+=1
            if fid not in families:counts['fixed_tail_empty_or_zero']+=1;continue
            lo=max((zero,Q(0)),((j+ell-1)*a-r,pattern[(j+ell-1)%3]))
            hi=min((3*a,Q(0)),(j*a+r,pattern[j%3]))
            if hi<lo or hi==lo:counts['first_axis_jet_empty']+=1;continue
            counts['active_raw']+=1
            if lo[0]==hi[0]:zero_openings.append({'family_id':fid,'start':j,'width_derivative':str(hi[1]-lo[1])})
            for c,d in (lo,hi):ties[c].add(d)
            intervals.append((lo,hi))
        interval_by_family[fid]=intervals
    assert counts['raw_candidates']==6144
    def union_jet(items):
        intervals=sorted(items);c=F(0);d=Q(0);left=right=None
        for lo,hi in intervals:
            if left is None:left,right=lo,hi
            elif right<lo:c+=right[0]-left[0];d+=right[1]-left[1];left,right=lo,hi
            else:right=max(right,hi)
        if left is not None:c+=right[0]-left[0];d+=right[1]-left[1]
        return c,d
    def mask_jet(items):
        cuts=sorted({x for interval in items for x in interval});answer=Q(0)
        for lo,hi in zip(cuts,cuts[1:]):
            if any(l<=lo and hi<=h for l,h in items):answer+=hi[1]-lo[1]
        return answer
    D=F(0);prefix=0;cp=P/f'{stage}_checkpoint.json';section_rows=[]
    if cp.exists():
        ck=load(cp);assert ck['source_sha256']==reg['source_sha256'] and ck['tail_histogram_sha256']==sha(P/'tail_histogram.json')
        prefix=ck['completed_sections'];D=dec(ck['partial_derivative_volume']);section_rows=ck['sections']
    def checkpoint(status):
        save(cp,{'status':status,'completed_sections':prefix,'total_sections':len(hist),'partial_derivative_volume':enc(D),'sections':section_rows,'tail_histogram_sha256':sha(P/'tail_histogram.json'),'source_sha256':reg['source_sha256'],'runtime':rt(start)})
    for index in range(prefix,len(hist)):
        mask,weight=hist[index];items=[]
        for fid,intervals in interval_by_family.items():
            if mask&(1<<fid):items.extend(intervals)
        _,speed=union_jet(items);assert speed==mask_jet(items)
        D+=weight*speed;section_rows.append({'family_mask':str(mask),'oneD_union_length_derivative':str(speed)})
        prefix=index+1
        if prefix%128==0:checkpoint('EXACT_section_prefix_saved')
    ratioD=D*Q(16,3);lo,hi=bounds(ratioD.a)
    status='STRICT_POSITIVE_DIRECTIONAL_DERIVATIVE' if lo>0 else 'STRICT_NEGATIVE_DIRECTIONAL_DERIVATIVE' if hi<0 else 'EXACT_ZERO_DIRECTIONAL_DERIVATIVE' if ratioD==F(0) else 'UNRESOLVED_SIGN'
    result={'status':status,'profile':stage,'pattern':list(map(str,pattern)),'position':'p_j=j*a1+epsilon*z_(jmod3);tailC unchanged','period_mass':'3','normalization':'16*physicalperiodunionvolume/3','common_radius':'r=N^(1/5)/2 independentofepsilon,same5coordinates','counts':dict(counts),'new_zero_width_openings':zero_openings,'same_base_different_speed_groups':sum(len(x)>1 for x in ties.values()),'directional_derivative_coefficients':enc(ratioD),'derivative_interval':[str(lo),str(hi)],'display_only':[float(lo),float(hi)],'both_oneD_union_engines_agree_all_sections':True,'baseline_not_evaluated':'C inherited from original full source plusperiod/reflection identity;no newbaselineunionjob','tail_histogram_sha256':sha(P/'tail_histogram.json'),'source_sha256':reg['source_sha256'],'registration_sha256':sha(P/'registration.json'),'root_denominator':str(den),'root_lower_numerators':roots,'runtime':rt(start),'scope':'one of two finiteperiod3position directions only;no fullfan or Cnclaim;main ownsproof,L1andcollar'}
    save(P/f'{stage}_derivative.json',result);checkpoint('COMPLETE_all_section_derivatives')
    print(json.dumps({k:v for k,v in result.items() if k not in ('directional_derivative_coefficients','derivative_interval')},indent=2),flush=True)
if __name__=='__main__':
    assert all(os.environ.get(k)=='1' for k in ENV) and os.getpriority(os.PRIO_PROCESS,0)>=10
    if len(sys.argv)!=2:raise SystemExit('usage: position_sections.py tail|plus|minus')
    main(sys.argv[1])
