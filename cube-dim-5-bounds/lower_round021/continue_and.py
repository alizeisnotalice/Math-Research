#!/usr/bin/env python3
"""Resume only the saved signature AND; never build geometry or histograms."""
from pathlib import Path
from fractions import Fraction as Q
from collections import defaultdict
import hashlib,json,os,time
from finite_union import G,A,ENV,runtime,volume,dec

HERE=Path(__file__).resolve().parent
PROJECT=HERE.parent
def load(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def save(p,d):
    tmp=p.with_suffix(p.suffix+'.tmp')
    tmp.write_text(json.dumps(d,separators=(',',':'))+'\n')
    tmp.replace(p)

def main():
    start=time.monotonic()
    reg=load(HERE/'repair_registration.json')
    assert sha(Path(__file__))==reg['source_sha256']
    for row in reg['locked_inputs']:
        assert sha(HERE/row['path'])==row['sha256'],row
    orig=load(HERE/'registration.json')
    for row in orig['inputs']:
        assert sha(PROJECT/row['path'])==row['sha256'],row
    geom=load(HERE/'candidate_geometry.json')
    plane=load(HERE/'plane_histogram.json')
    unrot=load(HERE/'unrotated_histogram.json')
    assert plane['status']=='COMPLETE_joint_planar_stripe_signatures'
    assert unrot['status']=='COMPLETE_3D_family_signatures'
    assert geom['a1']==str(A) and geom['grid']==G and geom['raw_case_count']==1024
    assert len(geom['raw_inventory'])==1024 and len(geom['new18_accounting'])==18
    assert all(x['not_ignored'] for x in geom['new18_accounting'])
    families=geom['families']
    assert len(families)==36 and [f['family_id'] for f in families]==list(range(36))
    ph=sorted((int(k),v) for k,v in plane['histogram'].items())
    bh=sorted((int(k),v) for k,v in unrot['histogram'].items())
    assert len(ph)==262 and len(bh)==649
    assert sum(v for _,v in ph)==G**2 and sum(v for _,v in bh)==G**3//8
    valid=(1<<len(families))-1
    assert all(0<=m<=valid and v>0 for m,v in ph+bh)
    formulas={(x['K'],x['n']):x for x in load(PROJECT/'lower_round020/results.json')['all_case_exact_formulas']}
    marginal_checks=[]
    for f in families:
        assert f['N']==f['K']*f['n']
        rl=Q(f['strict_common_radius'])
        assert (2*rl)**5<f['N']
        sl,sh=map(Q,geom['s_interval'])
        assert (2*sl-1)**5<2<(2*sh-1)**5
        assert Q(f['alpha_grid'],G)*A<=rl and Q(f['rho_grid'],G)*sh<=rl
        key=f['K'],f['n']
        marginal=Q(25*sum(v for m,v in ph if m&(1<<f['family_id'])),G**2)
        if key in formulas:
            aa,bb,cc=map(Q,formulas[key]['normalized_exact_area_polynomial'])
            rho=Q(f['rho_grid'],G)
            expected=aa*rho*rho+bb*rho+cc
            assert marginal==expected,(key,marginal,expected)
            marginal_checks.append({'K':key[0],'n':key[1],'rho_inner':str(rho),'area':str(marginal),'equals020_exact_polynomial':True})
        elif key in ((8,1),(8,2)):
            assert marginal==25
            marginal_checks.append({'K':key[0],'n':key[1],'area':'25','equals017018_wholeplane':True})
    required={tuple(x) for x in load(PROJECT/'lower_round020/results.json')['previous_partial_cases_resolved']}
    assert required<={(x['K'],x['n']) for x in marginal_checks}
    total_pairs=len(ph)*len(bh)
    ckpath=HERE/'and_checkpoint.json'
    joint=defaultdict(int);prefix=0;total=0
    if ckpath.exists():
        ck=load(ckpath)
        assert ck['source_sha256']==reg['source_sha256']
        assert ck['repair_registration_sha256']==sha(HERE/'repair_registration.json')
        prefix=ck['completed_pairs'];total=ck['union_grid_volume']
        joint.update({int(k):v for k,v in ck['joint_histogram'].items()})
        assert 0<=prefix<=total_pairs and sum(joint.values())==total
    initial_prefix=prefix
    def checkpoint(status):
        save(ckpath,{'status':status,'completed_pairs':prefix,'total_pairs':total_pairs,'next_pair_indices':list(divmod(prefix,len(bh))),'pair_order':'ascending numeric plane mask, then ascending numeric unrotated mask','union_grid_volume':total,'joint_histogram':{str(k):v for k,v in sorted(joint.items())},'source_sha256':reg['source_sha256'],'repair_registration_sha256':sha(HERE/'repair_registration.json'),'new_pairs_this_run':prefix-initial_prefix,'runtime':runtime(start)})
    print(json.dumps({'stage':'saved_histogram_AND_only','initial_prefix':prefix,'total_pairs':total_pairs,'runtime':runtime(start)}),flush=True)
    while prefix<total_pairs:
        if prefix-initial_prefix>=200000 or time.monotonic()-start>=120:
            checkpoint('CHECKPOINT_LIMIT_no_complete_ratio')
            print(json.dumps({'status':'CHECKPOINT_LIMIT_no_complete_ratio','completed_pairs':prefix,'runtime':runtime(start)}),flush=True)
            return
        i,j=divmod(prefix,len(bh));pm,pa=ph[i];bm,bv=bh[j];mask=pm&bm
        if mask:
            joint[mask]+=pa*bv
            total+=pa*bv
        prefix+=1
        if prefix%len(bh)==0:checkpoint('AND_exact_prefix_saved')
    assert sum(joint.values())==total and 0<=total<=G**5//8
    V=Q(total,G**5);shell=Q(0);shell_rows=[]
    for f in families:
        du=sum((Q(volume(ob)-volume(ib),G**3) for ib,ob in zip(f['inner_boxes'],f['outer_boxes'])),Q(0))
        vu=min(Q(1,8),sum((Q(volume(b),G**3) for b in f['outer_boxes']),Q(0)))
        lowrho=Q(f['rho_grid'],G);hirho=Q(f['rho_outer_interval_upper'])
        dp=min(Q(1),36*(hirho**2-lowrho**2))
        assert du>=0 and dp>=0
        err=du+vu*dp;shell+=err
        shell_rows.append({'family_id':f['family_id'],'unrotated_shell_volume_upper':str(du),'outer_unrotated_union_volume_upper':str(vu),'planar_stripe_threshold_shell_area_upper':str(dp),'product_difference_upper':str(err)})
    low=8*A*sl**4*V;inner_hi=8*A*sh**4*V
    high=8*A*sh**4*min(Q(1,8),V+shell)
    C=load(PROJECT/'lower_round003/a2_equals_a3_results.json');cl,ch=map(Q,C['ratio_interval'])
    old=load(PROJECT/'lower_round016/delivery007/result.json');ol,oh=map(Q,old['ratio_interval'])
    result={'status':'COMPLETE_one_new_finite_5D_union_candidate','a1':str(A),'tail_spacing':'s=(1+u)/2,u^5=2','raw_cases_accounted':1024,'new18_cases_accounted':geom['new18_accounting'],'included_K_n_family_count':len(families),'normalized_inner_stripe_union_volume':str(V),'exact_inner_union_ratio_expression':f'8*(151/100)*s^4*({V})','inner_union_ratio_interval':[str(low),str(inner_hi)],'nominal_finite_family_ratio_interval':[str(low),str(high)],'nominal_ratio_decimal_display':[dec(low),dec(high)],'normalized_nominal_minus_inner_union_upper':str(shell),'below_C_cert_for_this_finite_candidate':high<cl,'above_C_cert_certified':low>ch,'difference_from_C_cert_interval':[str(low-ch),str(high-cl)],'difference_from_old016_interval':[str(low-oh),str(high-ol)],'old016_is_different_a1_and_only26_families':True,'strictly_above_old016_certified':low>oh,'exact_marginal_checks':marginal_checks,'joint_union_rule':'plane family mask AND unrotated family mask nonzero, count each pair once','normalization':'periodmass25 cancels period mass denominator;8 unrotated reflections;Jacobian a1*s^4;ratio=8*a1*s^4*stripevolume','shell_bound_proof':'225 window annuli bound wholeperiod area by900*(rho_hi^2-rho_inner^2);divide by25 gives stripe upper36*(...);perfamily product error<=unrotated_shell+outer_unrotated_volume*planar_shell;sum only error terms','shell_rows':shell_rows,'joint_histogram':{str(k):v for k,v in sorted(joint.items())},'event_component_indices':[geom['event_component_containing_candidate']['left_event']['event_index'],geom['event_component_containing_candidate']['right_event']['event_index']],'geometry_sha256':sha(HERE/'candidate_geometry.json'),'plane_histogram_sha256':sha(HERE/'plane_histogram.json'),'unrotated_histogram_sha256':sha(HERE/'unrotated_histogram.json'),'original_source_sha256':sha(HERE/'finite_union.py'),'source_sha256':reg['source_sha256'],'registration_sha256':sha(HERE/'registration.json'),'repair_registration_sha256':sha(HERE/'repair_registration.json'),'completed_AND_pairs':prefix,'new_AND_pairs_this_run':prefix-initial_prefix,'recomputed_geometry_or_histogram_cells':0,'runtime':runtime(start),'scope':'one fixed finite-family candidate only;no second candidate;no allscale optimum or C5 operator upper;main owns L1 and collar qualification'}
    save(HERE/'results.json',result)
    checkpoint('COMPLETE_saved_histogram_AND')
    print(json.dumps({k:v for k,v in result.items() if k not in ('joint_histogram','shell_rows','new18_cases_accounted','difference_from_C_cert_interval','difference_from_old016_interval')},indent=2),flush=True)

if __name__=='__main__':
    assert all(os.environ.get(k)=='1' for k in ENV)
    assert os.getpriority(os.PRIO_PROCESS,0)>=10
    main()
