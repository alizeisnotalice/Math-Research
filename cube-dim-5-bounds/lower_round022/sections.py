#!/usr/bin/env python3
"""One static 4D shape-section histogram, then exact affine first-axis interfaces."""
from pathlib import Path
from fractions import Fraction as Q
from collections import defaultdict
from decimal import Decimal,localcontext
import hashlib,json,os,resource,sys,threading,time,signal
P=Path(__file__).resolve().parent;ROOT=P.parent;G=5*2**40
ENV=('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','VECLIB_MAXIMUM_THREADS','BLIS_NUM_THREADS')
LEFT=Q(3,2);RIGHT=Q(2)
def load(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def save(p,d):
    t=p.with_suffix(p.suffix+'.tmp');t.write_text(json.dumps(d,separators=(',',':'))+'\n');t.replace(p)
def volume(b):
    z=1
    for i in range(0,len(b),2):z*=max(0,b[i+1]-b[i])
    return z
def rt(start):
    rss=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return {'pid':os.getpid(),'nice':os.getpriority(os.PRIO_PROCESS,0),'threads':threading.active_count(),'peak_rss_bytes':rss*(1 if sys.platform=='darwin' else 1024),'wall_seconds':time.monotonic()-start,'thread_environment':{k:os.environ.get(k) for k in ENV}}
def dec(q):
    with localcontext() as c:c.prec=60;return str(Decimal(q.numerator)/Decimal(q.denominator))
def interval(ell,r,a):
    # Odd blocks are left anchored; even blocks are right anchored after clipping.
    length=max(Q(0),min(a/2,r-Q(ell-1,2)*a))
    return (Q(0),length) if ell%2 else (a/2-length,a/2)
def union_length(desc,a):
    left=max((interval(ell,r,a)[1] for ell,r in desc if ell%2),default=Q(0))
    right=max((a/2-interval(ell,r,a)[0] for ell,r in desc if not ell%2),default=Q(0))
    return min(a/2,left+right)
def pieces(desc):
    # Every raw interval endpoint is affine. Their crossing values include all
    # changes of clipping, emptiness, nesting, overlap and union endpoint order.
    lines={(Q(0),Q(0)),(Q(1,2),Q(0))}
    for ell,r in desc:
        k=ell//2 if ell%2==0 else (ell-1)//2
        lines.add((Q(k),-r))
        lines.add((Q(-((ell-1)//2)),r))
    lines=sorted(lines);events={LEFT,RIGHT}
    for i,(s,t) in enumerate(lines):
        for ss,tt in lines[i+1:]:
            if s!=ss:
                a=(tt-t)/(s-ss)
                if LEFT<a<RIGHT:events.add(a)
    ev=sorted(events);out=[]
    for a,b in zip(ev,ev[1:]):
        va,vb=union_length(desc,a),union_length(desc,b)
        slope=(vb-va)/(b-a);intercept=va-slope*a
        assert union_length(desc,(a+b)/2)==slope*(a+b)/2+intercept
        row={'left':str(a),'right':str(b),'slope':str(slope),'intercept':str(intercept)}
        if out and row['slope']==out[-1]['slope'] and row['intercept']==out[-1]['intercept']:out[-1]['right']=str(b)
        else:out.append(row)
    maximum=max(union_length(desc,a) for a in ev)
    return out,maximum

def main():
    start=time.monotonic();reg=load(P/'registration.json')
    def hard_wall(signum,frame):
        save(P/'limit_result.json',{'status':'WALL_LIMIT_180seconds_outputs_remain_recoverable','no_complete_bound_claim':not (P/'results.json').exists(),'runtime':rt(start)})
        raise SystemExit(0)
    signal.signal(signal.SIGALRM,hard_wall);signal.alarm(180)
    assert sha(Path(__file__))==reg['source_sha256']
    for r in reg['inputs']:assert sha(ROOT/r['path'])==r['sha256']
    geom=load(ROOT/'lower_round021/candidate_geometry.json');plane=load(ROOT/'lower_round021/plane_histogram.json')
    assert geom['raw_case_count']==1024 and len(geom['families'])==36 and len(geom['new18_accounting'])==18
    assert plane['status']=='COMPLETE_joint_planar_stripe_signatures'
    families=geom['families'];shapes=[];family_shapes=defaultdict(int)
    for f in families:
        for shape,ib,ob in zip(f['shapes'],f['inner_boxes'],f['outer_boxes']):
            i=len(shapes);shapes.append({'shape_id':i,'family_id':f['family_id'],'shape':shape,'r_lower':f['strict_common_radius'],'r_upper':f['nominal_radius_interval'][1],'tail_inner_box':ib[2:],'tail_outer_box':ob[2:]})
            family_shapes[f['family_id']]|=1<<i
    assert len(shapes)==54
    tailpath=P/'tail_histogram.json'
    if tailpath.exists():tail=load(tailpath)
    else:
        axes=[sorted({0,G//2}|{s['tail_inner_box'][2*q+k] for s in shapes for k in (0,1)}) for q in range(2)]
        th=defaultdict(int)
        for lo,hi in zip(axes[0],axes[0][1:]):
            for low,high in zip(axes[1],axes[1][1:]):
                mask=sum(1<<s['shape_id'] for s in shapes if s['tail_inner_box'][0]<=lo and hi<=s['tail_inner_box'][1] and s['tail_inner_box'][2]<=low and high<=s['tail_inner_box'][3])
                th[mask]+=(hi-lo)*(high-low)
        assert sum(th.values())==G**2//4
        tail={'status':'COMPLETE_two_tail_axes_shape_signatures','shapes':shapes,'axes':axes,'cells':(len(axes[0])-1)*(len(axes[1])-1),'histogram':{str(k):v for k,v in sorted(th.items())},'runtime':rt(start)};save(tailpath,tail)
    assert tail['shapes']==shapes
    ph=[]
    for k,v in plane['histogram'].items():
        family_mask=int(k);sm=0
        for f in families:
            if family_mask&(1<<f['family_id']):sm|=family_shapes[f['family_id']]
        ph.append((int(k),sm,v))
    ph.sort();th=sorted((int(k),v) for k,v in tail['histogram'].items())
    assert sum(v for _,_,v in ph)==G**2 and sum(v for _,v in th)==G**2//4
    joint=defaultdict(int);prefix=0;ckpath=P/'checkpoint.json'
    if ckpath.exists():
        ck=load(ckpath);assert ck['source_sha256']==reg['source_sha256'] and ck['tail_sha256']==sha(tailpath)
        prefix=ck['completed_pairs'];joint.update({int(k):v for k,v in ck['joint_histogram'].items()})
    initial=prefix;total_pairs=len(ph)*len(th)
    def checkpoint(status):
        save(ckpath,{'status':status,'completed_pairs':prefix,'total_pairs':total_pairs,'next_pair_indices':list(divmod(prefix,len(th))),'pair_order':'ascending numeric plane family mask, then ascending numeric two-tail shape mask','joint_histogram':{str(k):v for k,v in sorted(joint.items())},'new_pairs_this_run':prefix-initial,'source_sha256':reg['source_sha256'],'tail_sha256':sha(tailpath),'runtime':rt(start)})
    print(json.dumps({'stage':'fourD_joint_shape_sections','tail_cells':tail['cells'],'tail_signatures':len(th),'plane_signatures':len(ph),'total_pairs':total_pairs,'initial_prefix':prefix,'runtime':rt(start)}),flush=True)
    while prefix<total_pairs:
        if prefix-initial>=200000 or time.monotonic()-start>=180:
            checkpoint('CHECKPOINT_LIMIT_no_complete_window_bound');return
        i,j=divmod(prefix,len(th));_,pm,pv=ph[i];tm,tv=th[j]
        joint[pm&tm]+=pv*tv;prefix+=1
        if prefix%len(th)==0:checkpoint('FOURD_exact_prefix_saved')
    assert sum(joint.values())==G**4//4
    checkpoint('COMPLETE_fourD_shape_sections')
    save(P/'section_histogram.json',{'status':'COMPLETE_static_fourD_shape_sections','histogram':{str(k):v for k,v in sorted(joint.items())},'weight_denominator':str(G**4),'total_section_weight':'1/4','nonzero_shape_sections':sum(k!=0 for k in joint),'plane_sha256':sha(ROOT/'lower_round021/plane_histogram.json'),'tail_sha256':sha(tailpath),'source_sha256':reg['source_sha256'],'runtime':rt(start)})
    descriptors=defaultdict(int)
    for mask,w in joint.items():
        best={}
        for s in shapes:
            if mask&(1<<s['shape_id']):
                ell=s['shape'][0];best[ell]=max(best.get(ell,Q(0)),Q(s['r_lower']))
        descriptors[tuple(sorted(best.items()))]+=w
    rows=[];bound=Q(0)
    for desc,weight in sorted(descriptors.items()):
        pw,maximum=pieces(desc);w=Q(weight,G**4);bound+=w*maximum
        rows.append({'descriptor':[[ell,str(r)] for ell,r in desc],'weight':str(w),'segments':pw,'independent_max_length':str(maximum)})
    sl,sh=map(Q,geom['s_interval']);error=Q(0);shells=[]
    for f in families:
        rl,rh=map(Q,f['nominal_radius_interval']);du=Q(0);outer=Q(0)
        for s in shapes:
            if s['family_id']!=f['family_id']:continue
            ai=Q(volume(s['tail_inner_box']),G**2);ao=Q(volume(s['tail_outer_box']),G**2)
            assert ao>=ai
            # interval radius enlargement changes union length by at most2*dr.
            du+=2*(rh-rl)*ao+RIGHT/2*(ao-ai);outer+=RIGHT/2*ao
        vu=min(RIGHT/8,outer)
        rho=Q(f['rho_grid'],G);rhohi=Q(f['rho_outer_interval_upper'])
        dp=min(Q(1),36*(rhohi**2-rho**2));err=du+vu*dp;error+=err
        shells.append({'family_id':f['family_id'],'uniform_three_axis_shell_upper':str(du),'uniform_outer_three_axis_volume_upper':str(vu),'planar_stripe_shell_area_upper':str(dp),'uniform_product_error_upper':str(err)})
    frozen_upper=8*sh**4*bound;nominal_upper=8*sh**4*(bound+error)
    C=load(ROOT/'lower_round003/a2_equals_a3_results.json');cl,ch=map(Q,C['ratio_interval'])
    result={'status':'COMPLETE_static_sections_and_uniform_window_envelope','window':[str(LEFT),str(RIGHT)],'raw_inventory_source_count':1024,'shapes':54,'families':36,'new18_accounting':geom['new18_accounting'],'tail_histogram_cells':tail['cells'],'tail_signature_count':len(th),'plane_signature_count':len(ph),'completed_fourD_pairs':prefix,'new_pairs_this_run':prefix-initial,'section_signature_count':len(joint),'descriptor_count':len(rows),'descriptors':rows,'frozen_functional_expression':'R_inner(a)=8*s^4*sum_section weight*union_length_of_physical_first_axis_intervals(a)','per_shape_first_axis_interval':'[max(0,max(I_ell)*a-r_lower),min(a/2,min(I_ell)*a+r_lower)],empty_if_lower>=upper','frozen_sectionwise_independent_max_upper':str(frozen_upper),'uniform_nominal_error_before_8s4':str(error),'uniform_nominal_finitefunctional_upper':str(nominal_upper),'uniform_upper_decimal_display':dec(nominal_upper),'below_complete_C_cert_on_entire_window':nominal_upper<cl,'strict_C_low_minus_window_upper':str(cl-nominal_upper),'C_gap_decimal_display':dec(cl-nominal_upper),'shell_rows':shells,'normalization':'physical firstaxis contributes length without a1 Jacobian;other four axes Jacobian s^4;periodmass25 cancels25 stripes;8 unrotated reflections','upper_relaxation':'each section independently takes max_a unionlength before weighted sum;different sections may maximize at different a;this is an upper bound,not a realizable witness ratio','common_radius':'same stored strict r_lower for all5coordinates per family;firstaxis now exact physical affine,other4axes are frozen inward sections;nominal enlargement separately upper bounded','source_sha256':reg['source_sha256'],'registration_sha256':sha(P/'registration.json'),'section_histogram_sha256':sha(P/'section_histogram.json'),'tail_histogram_sha256':sha(tailpath),'runtime':rt(start),'scope':'one finite family window envelope and exact persection PWA interface;no a1 point scan,no allscale or operator optimum,no repeated021 union;main retains L1 and collar review'}
    save(P/'results.json',result)
    signal.alarm(0)
    print(json.dumps({k:v for k,v in result.items() if k not in ('descriptors','shell_rows','new18_accounting','strict_C_low_minus_window_upper')},indent=2),flush=True)

if __name__=='__main__':
    assert all(os.environ.get(k)=='1' for k in ENV) and os.getpriority(os.PRIO_PROCESS,0)>=10
    main()
