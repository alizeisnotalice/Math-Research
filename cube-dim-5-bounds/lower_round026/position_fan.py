#!/usr/bin/env python3
"""Full signed-D3 position fan: raw velocity events, two affine section unions."""
from pathlib import Path
from fractions import Fraction as Q
from functools import lru_cache
from itertools import product,combinations,permutations
from math import prod
from collections import defaultdict,Counter
import sys,os,time,threading,resource,hashlib,json
sys.dont_write_bytecode=True
P=Path(__file__).resolve().parent;ROOT=P.parent;OLD=ROOT/'lower_round025'
sys.path.insert(0,str(ROOT/'lower'))
from certify_rectangular import field
ENV=('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','VECLIB_MAXIMUM_THREADS','BLIS_NUM_THREADS')
def load(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def save(p,d):
    t=p.with_suffix(p.suffix+'.tmp');t.write_text(json.dumps(d,separators=(',',':'))+'\n');t.replace(p)
def enc(x):return {f'{i//5},{i%5}':str(v) for i,v in enumerate(x.a) if v}
def rt(start):
    rss=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return {'pid':os.getpid(),'nice':os.getpriority(os.PRIO_PROCESS,0),'threads':threading.active_count(),'peak_rss_bytes':rss*(1 if sys.platform=='darwin' else 1024),'wall_seconds':time.monotonic()-start,'thread_environment':{k:os.environ.get(k) for k in ENV}}

def main(stage):
    start=time.monotonic();reg=load(P/'registration.json')
    assert sha(Path(__file__))==reg['source_sha256']
    for row in reg['inputs']:assert sha(Path(row['path']))==row['sha256']
    F,u,v,bounds,_,_=field();zero=F(0);a=1+u**3/2
    def dec(d):
        z=[Q(0)]*25
        for ij,c in d.items():i,j=map(int,ij.split(','));z[5*i+j]=Q(c)
        return F(a=z)
    @lru_cache(None)
    def radius(n):
        x=n;i=j=0
        while x%2==0:x//=2;i+=1
        while x%3==0:x//=3;j+=1
        assert x==1
        r=u**i*v**j/2;assert 32*r**5==F(n);return r
    geo=load(OLD/'tail_geometry.json');hh=load(OLD/'tail_histogram.json');projection=load(OLD/'projection_certificate.json')
    keep=set(projection['kept_first_span_families']);keepmask=sum(1<<x for x in keep)
    hist=sorted((int(m),dec(c)) for m,c in hh['histogram'].items())
    assert len(geo['raw_inventory'])==1024 and len(hist)==2224 and all(m&~keepmask==0 for m,_ in hist)
    affine=((Q(1,2),Q(-1,2)),(Q(0),Q(1)),(Q(-1,2),Q(-1,2)),(Q(0),Q(0)))
    print({'stage':stage,'pid':os.getpid(),'raw_objects':6144,'static_sections':len(hist),'runtime':rt(start)},flush=True)
    if stage=='events':
        assert not (P/'events.json').exists(),'eventcertificatealreadycomplete'
        groups=defaultdict(set);counts=Counter();inventory=[]
        for row in geo['raw_inventory']:
            fid=row['family_id'];ell=row['lengths'][0];r=radius(row['N'])
            for j in range(-1,5):
                counts['all_raw']+=1;rawlo=(j+ell-1)*a-r;rawhi=j*a+r;lo=max(zero,rawlo);hi=min(3*a,rawhi)
                lower_labels={3} if rawlo<zero else {(j+ell-1)%3} if zero<rawlo else {3,(j+ell-1)%3}
                upper_labels={j%3} if rawhi<3*a else {3} if 3*a<rawhi else {3,j%3}
                groups[lo].update(lower_labels);groups[hi].update(upper_labels)
                status='tail_fixed_empty_or_zero' if row['status']!='positive_tail_box' else 'first_axis_strict_negative' if hi<lo else 'first_axis_zero' if hi==lo else 'positive_base_raw'
                counts[status]+=1
                inventory.append({'family_id':fid,'start':j,'first_phase':j%3,'last_phase':(j+ell-1)%3,'lower_velocity_labels':sorted(lower_labels),'upper_velocity_labels':sorted(upper_labels),'status':status})
        assert counts['all_raw']==6144 and counts['positive_base_raw']==179
        event_witness={};pairs=0
        for base,labels in groups.items():
            for l,r in combinations(sorted(labels),2):
                pairs+=1
                if pairs>reg['max_event_comparisons']:save(P/'event_checkpoint.json',{'status':'FINITE_PAIR_LIMIT','completed':pairs-1,'runtime':rt(start)});return
                A,B=affine[l];C,D=affine[r]
                if B==D:continue
                t=(C-A)/(B-D)
                if 0<=t<=Q(1,3):event_witness.setdefault(t,{'base':enc(base),'velocity_labels':[l,r],'affine_left':list(map(str,(A,B))),'affine_right':list(map(str,(C,D)))})
        # Sign reversal on the minus halfedge multiplies both velocities by-1,
        # hence has exactly the same equality events, including clipping0.
        events=sorted(event_witness);internal=[t for t in events if 0<t<Q(1,3)]
        out={'status':'COMPLETE_all6144_raw_and2224section_velocity_events','counts':dict(counts),'all_raw_inventory':inventory,'base_groups':len(groups),'velocity_label_set':['z0','z1','z2','0'],'pair_comparisons':pairs,'events':[{'t':str(t),'witness':event_witness[t]} for t in events],'internal_events':list(map(str,internal)),'minus_halfedge_same_equalities_by_sign_reversal':True,'static_sections_checked':len(hist),'section_raw_interval_memberships':sum(m.bit_count()*6 for m,_ in hist),'all_sections_subset_of_full_raw_inventory':True,'source_sha256':reg['source_sha256'],'runtime':rt(start)}
        save(P/'events.json',out);print({k:v for k,v in out.items() if k!='all_raw_inventory'},flush=True);return
    event=load(P/'events.json');assert event['internal_events']==[], 'unexpectedinternalevents;registeractualsubcellsbeforeunions'
    if stage in ('plus','minus'):
        assert not (P/f'{stage}_affine.json').exists(),'completeprofileexists'
        sign=1 if stage=='plus' else -1;t=Q(reg['sample_t']);forms=tuple((sign*A,sign*B) for A,B in affine)
        def endpoint(c,label):
            A,B=forms[label];return c,A+B*t,A,B
        def key(x):return x[:2]
        def le(x,y):return key(x)<=key(y)
        def lt(x,y):return key(x)<key(y)
        intervals={};counts=Counter();zero_openings=[]
        for row in geo['raw_inventory']:
            fid=row['family_id'];ell=row['lengths'][0];r=radius(row['N']);items=[]
            for j in range(-1,5):
                counts['all_raw']+=1
                if row['status']!='positive_tail_box':counts['tail_fixed_empty_or_zero']+=1;continue
                lo=max((endpoint(zero,3),endpoint((j+ell-1)*a-r,(j+ell-1)%3)),key=key)
                hi=min((endpoint(3*a,3),endpoint(j*a+r,j%3)),key=key)
                if not lt(lo,hi):counts['first_axis_jet_empty']+=1;continue
                counts['active_raw']+=1
                if lo[0]==hi[0]:zero_openings.append({'family_id':fid,'start':j})
                items.append((lo,hi))
            intervals[fid]=items
        def merge(items):
            A=B=d=Q(0);left=right=None
            for lo,hi in sorted(items,key=lambda x:(key(x[0]),key(x[1]))):
                if left is None:left,right=lo,hi
                elif lt(right,lo):A+=right[2]-left[2];B+=right[3]-left[3];d+=right[1]-left[1];left,right=lo,hi
                else:right=max((right,hi),key=key)
            if left is not None:A+=right[2]-left[2];B+=right[3]-left[3];d+=right[1]-left[1]
            assert A+B*t==d;return A,B
        def masks(items):
            unique={}
            for pair in items:
                for x in pair:
                    if key(x) in unique:assert unique[key(x)][2:]==x[2:]
                    unique[key(x)]=x
            cuts=sorted(unique.values(),key=key);A=B=Q(0)
            for lo,hi in zip(cuts,cuts[1:]):
                if any(le(l,lo) and le(hi,h) for l,h in items):A+=hi[2]-lo[2];B+=hi[3]-lo[3]
            return A,B
        outA=F(0);outB=F(0);prefix=0;rows=[];cp=P/f'{stage}_checkpoint.json'
        if cp.exists():
            ck=load(cp);assert ck['source_sha256']==reg['source_sha256']
            prefix=ck['completed_sections'];outA=dec(ck['partial_A']);outB=dec(ck['partial_B']);rows=ck['section_affine_derivatives']
        def checkpoint(status):
            save(cp,{'status':status,'completed_sections':prefix,'total_sections':len(hist),'partial_A':enc(outA),'partial_B':enc(outB),'section_affine_derivatives':rows,'source_sha256':reg['source_sha256'],'runtime':rt(start)})
        for i in range(prefix,len(hist)):
            mask,weight=hist[i];items=[]
            for fid,values in intervals.items():
                if mask&(1<<fid):items.extend(values)
            A,B=merge(items);assert (A,B)==masks(items);outA+=weight*A;outB+=weight*B
            rows.append({'family_mask':str(mask),'constant':str(A),'slope':str(B)});prefix=i+1
            if prefix%128==0:checkpoint('EXACT_affine_section_prefix_saved')
        A=outA*Q(16,3);B=outB*Q(16,3);D=A+B*t;lo,hi=bounds(D.a)
        result={'status':'COMPLETE_exact_affine_halfedge','profile':stage,'sample_t':str(t),'pattern_at_sample':list(map(str,(sign*(1-t)/2,sign*t,-sign*(1+t)/2))),'D_constant':enc(A),'D_slope':enc(B),'sample_derivative_interval':list(map(str,(lo,hi))),'sample_display_only':[float(lo),float(hi)],'counts':dict(counts),'zero_width_births':zero_openings,'all2224_sections_dual_oneD_engines_agree_A_B':True,'no_baseline_or_static_histogram_recomputation':True,'source_sha256':reg['source_sha256'],'events_sha256':sha(P/'events.json'),'runtime':rt(start)}
        save(P/f'{stage}_affine.json',result);checkpoint('COMPLETE_all_affine_sections');print(result,flush=True);return
    assert stage=='finalize'
    plus=load(P/'plus_affine.json');minus=load(P/'minus_affine.json');PA,PB=dec(plus['D_constant']),dec(plus['D_slope']);MA,MB=dec(minus['D_constant']),dec(minus['D_slope'])
    oldplus=dec(load(OLD/'plus_derivative.json')['directional_derivative_coefficients']);oldminus=dec(load(OLD/'minus_derivative.json')['directional_derivative_coefficients'])
    assert PA==oldplus/2 and MA==oldminus/2 and PA+PB/3==MA+MB/3
    values=[PA,MA,PA+PB/3,MA+MB/3];bounds_all=[bounds(z.a) for z in values]
    maximum=max(values);lo,hi=bounds(maximum.a)
    status='STRICT_NEGATIVE_full_position_fan' if hi<0 else 'POSITIVE_direction_exists' if lo>0 else 'NONPOSITIVE_full_position_fan' if all(x<=zero for x in values) else 'UNRESOLVED_sign'
    # Groupqualification is on actual signed permutations, not the massfan S3.
    perms=list(permutations(range(3)))
    def parity(p):return -1 if sum(p[i]>p[j] for i in range(3) for j in range(i+1,3))%2 else 1
    def act(z,p):return tuple(parity(p)*z[p[i]] for i in range(3))
    middle=(Q(1,2),Q(0),Q(-1,2));op={act(middle,p) for p in perms};om={act(tuple(-x for x in middle),p) for p in perms};assert len(op)==len(om)==3 and op.isdisjoint(om)
    endpoint=(Q(1,3),Q(1,3),Q(-2,3));assert tuple(-x for x in endpoint) in {act(endpoint,p) for p in perms}
    for oldname in ('lower_round023','lower_round024','lower_round025'):
        old=ROOT/oldname
        for x in load(old/'archive_whitelist.json')['files']:assert sha(old/x['path'])==x['sha256']
    result={'status':status,'full_fan_domain':'range1,sum0,bothsignedD3halfedges +/-((1-t)/2,t,-(1+t)/2),0<=t<=1/3','halfedge_formulas':{'plus':{'A':enc(PA),'B':enc(PB)},'minus':{'A':enc(MA),'B':enc(MB)}},'closed_endpoint_intervals':[[str(l),str(h)] for l,h in bounds_all],'full_fan_maximum_coefficients':enc(maximum),'full_fan_maximum_interval':[str(lo),str(hi)],'display_only':[float(lo),float(hi)],'zero_direction':'0','positive_homogeneity':'D(lambda*z)=lambda*D(z),lambda>=0','sample_t0_reuses025_after_divide_by_range2':True,'t1_3_plus_minus_equivalence_under_signedG_verified':True,'bothhalfedgesconstant':PB==MB==zero,'actual_internal_events':event['internal_events'],'closed_extension':'finiteendpointLipschitz epsilon*directiondifference impliesD Lipschitz;affineopencellformulaextends toclosedt0,t1/3','old023024025_archives_unchanged':True,'new_union_jobs_in_finalizer':0,'source_sha256':reg['source_sha256'],'runtime':rt(start),'scope':'onlylocalfirstorderpositionsofregistered1024lengthfinitefunctional;noCn/generalinputupper,noquantifiedepsilon;mainownsproofL1collar'}
    save(P/'fan_results.json',result);print(result,flush=True)
if __name__=='__main__':
    assert all(os.environ.get(k)=='1' for k in ENV) and os.getpriority(os.PRIO_PROCESS,0)>=10
    if len(sys.argv)!=2:raise SystemExit('usage: position_fan.py events|plus|minus|finalize')
    main(sys.argv[1])
