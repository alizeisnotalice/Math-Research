#!/usr/bin/env python3
"""Exact fixed-fine-record columns for full actual n1 E and common beta.
Different k at one z are never merged into a column.
"""
from pathlib import Path
from fractions import Fraction as F
from collections import defaultdict
import argparse, importlib.util, json, hashlib, time, sys
sys.dont_write_bytecode=True
KINDS=('all','sigma=1','sigma>=1/2','sigma<=1/4')


def module(path,name):
    s=importlib.util.spec_from_file_location(name,path)
    m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m


def outside_length(z,xa,xb,R):
    return xb-xa-max(F(0),min(xb,z+R)-max(xa,z-R))


def evaluate(flow,rec,glob,old,batch,label,atoms,alpha,params):
    began=time.monotonic()
    saved,radii,cells,_=flow.context(old,atoms,alpha)
    I=alpha/4;lo=alpha/4;hi=alpha/2
    Evolume=F(saved['eligible_volume']);X=alpha*Evolume
    data=[]
    for a,b,J,trace in cells:
        mid=(a+b)/2;records=rec.intervals(trace,lo,hi)
        members={k:frozenset(i for i,(y,w) in enumerate(atoms) if abs(y-mid)<radii[k])
                 for k in range(1,saved['D']+1)}
        if data and data[-1][1]==a and data[-1][2:]==(records,members):
            data[-1]=(data[-1][0],b,records,members)
        else:data.append((a,b,records,members))
    assert sum((b-a for a,b,records,members in data),F(0))==Evolume
    O_direct=O_columns=M=Jaccard_O=tau_O=F(0)
    first_union_checks=0
    costs={kind:F(0) for kind in KINDS}
    max_column={kind:F(0) for kind in KINDS}
    max_conditional={kind:F(0) for kind in KINDS}
    witnesses={kind:None for kind in KINDS}
    columns=[];fine_groups={};breakpoint_checks=0;active_edges=0
    gram_local=None
    for zc,(za,zb,iz,mz) in enumerate(data):
        for k,(ka,kb) in iz.items():
            fine_set=mz[k]
            fine_mass=sum((atoms[i][1] for i in fine_set),F(0))
            Vk=2*radii[k];gk=fine_mass/Vk;fine_probability=(kb-ka)/I
            assert alpha/4<gk<=alpha
            M+=gk*fine_probability*(zb-za)
            entries=[];area_integral={kind:F(0) for kind in KINDS}
            common_groups=defaultdict(F)
            for xc,(xa,xb,ix,mx) in enumerate(data):
                for j,(ja,jb) in ix.items():
                    if j>=k:continue
                    W=min(jb,kb)-max(ja,ka)
                    if W<=0:continue
                    common_set=mx[j]&fine_set
                    common=sum((atoms[i][1] for i in common_set),F(0))
                    if not common:continue
                    R=radii[j];Vj=2*R
                    out=(xb-xa)*(zb-za)-glob.strip_area(xa,xb,za,zb,R)
                    if out==0:continue
                    fine_lag=glob.strip_area(xa,xb,za,zb,R+radii[k])-glob.strip_area(xa,xb,za,zb,R)
                    assert fine_lag==out
                    sigma=common/fine_mass
                    coarse_mass=sum((atoms[i][1] for i in mx[j]),F(0))
                    union_set=mx[j]|fine_set
                    union_mass=sum((atoms[i][1] for i in union_set),F(0))
                    assert union_mass==coarse_mass+fine_mass-common
                    jaccard=common/union_mass;tau=common/coarse_mass
                    # Conditional on first iid hit of the union, its source i
                    # has exact probability w_i/P(union); no sampling is used.
                    first_union_same_color=sum((atoms[i][1]/union_mass for i in common_set),F(0))
                    assert first_union_same_color==jaccard
                    assert 0<jaccard<=tau<=1
                    if sigma==1:assert jaccard==tau
                    first_union_checks+=1
                    coefficient=W/I*sigma/Vj
                    flags=dict(all=True,**{'sigma=1':sigma==1,'sigma>=1/2':sigma>=F(1,2),
                                          'sigma<=1/4':sigma<=F(1,4)})
                    integral=coefficient*out
                    O_direct+=W*common*out/(I*Vj*Vk)
                    Jaccard_O+=jaccard*gk*integral
                    tau_O+=tau*gk*integral
                    active_edges+=1
                    for kind,flag in flags.items():
                        if flag:area_integral[kind]+=integral
                    common_groups[tuple(sorted(common_set))]+=integral
                    entries.append((j,xa,xb,R,coefficient,flags))
            knots={za,zb}
            for j,xa,xb,R,coefficient,flags in entries:
                for x in (xa,xb):
                    for sign in (-1,1):
                        z=x+sign*R
                        if za<z<zb:knots.add(z)
            maxima={kind:F(0) for kind in KINDS}
            witness_z={kind:None for kind in KINDS}
            integrated={kind:F(0) for kind in KINDS}
            previous_z=previous_value=None
            for z in sorted(knots):
                values={kind:F(0) for kind in KINDS}
                for j,xa,xb,R,coefficient,flags in entries:
                    value=coefficient*outside_length(z,xa,xb,R)
                    for kind,flag in flags.items():
                        if flag:values[kind]+=value
                for kind in KINDS:
                    if values[kind]>maxima[kind]:
                        maxima[kind]=values[kind];witness_z[kind]=z
                    if previous_z is not None:
                        integrated[kind]+=(z-previous_z)*(previous_value[kind]+values[kind])/2
                previous_z=z;previous_value=values;breakpoint_checks+=1
            assert integrated==area_integral
            conditional_max={kind:maxima[kind]/fine_probability for kind in KINDS}
            for kind in KINDS:
                costs[kind]+=gk*integrated[kind]
                max_column[kind]=max(max_column[kind],maxima[kind])
                max_conditional[kind]=max(max_conditional[kind],conditional_max[kind])
                if maxima[kind] and max_column[kind]==maxima[kind]:
                    witnesses[kind]=dict(zcell=zc,k=k,z_one_sided_limit=str(witness_z[kind]),
                                         z_interval=[str(za),str(zb)],record_interval=[str(ka),str(kb)],
                                         fine_source_set=sorted(fine_set),fine_mass=str(fine_mass),
                                         raw_column_max=str(maxima[kind]),
                                         conditional_column_max=str(conditional_max[kind]),
                                         column_O=str(gk*integrated[kind]))
            O_columns+=gk*integrated['all']
            if 'gram_k' in params:
                z0=-F(3,5)*params['radius']
                if k==params['gram_k'] and za<z0<zb:
                    value=sum((coefficient*outside_length(z0,xa,xb,R)
                               for j,xa,xb,R,coefficient,flags in entries),F(0))
                    gram_local=dict(k=k,z0=str(z0),c_k_z0=str(value),
                                    conditional_c_k_z0=str(value/fine_probability),
                                    gk=str(gk),fine_record_probability=str(fine_probability))
            group_key=(k,tuple(sorted(fine_set)))
            if group_key not in fine_groups:
                fine_groups[group_key]=dict(k=k,fine_source_set=sorted(fine_set),fine_mass=fine_mass,
                  gk=gk,z_volume=F(0),column_integral=F(0),O=F(0),max=F(0),
                  sigma1_O=F(0),common_subsets=defaultdict(F),record_columns=0)
            group=fine_groups[group_key]
            group['z_volume']+=zb-za;group['column_integral']+=integrated['all']
            group['O']+=gk*integrated['all'];group['sigma1_O']+=gk*integrated['sigma=1']
            group['max']=max(group['max'],maxima['all']);group['record_columns']+=1
            for subset,value in common_groups.items():group['common_subsets'][subset]+=gk*value
            if integrated['all']:
                assert integrated['all']/(zb-za)<=maxima['all']
                columns.append(dict(zcell=zc,k=k,z_interval=[str(za),str(zb)],
                  record_interval=[str(ka),str(kb)],fine_record_probability=str(fine_probability),
                  fine_source_set=sorted(fine_set),fine_mass=str(fine_mass),gk=str(gk),
                  column_integral={kind:str(v) for kind,v in integrated.items()},
                  column_mean={kind:str(v/(zb-za)) for kind,v in integrated.items()},
                  column_max={kind:str(v) for kind,v in maxima.items()},
                  conditional_column_max={kind:str(v) for kind,v in conditional_max.items()},
                  O={kind:str(gk*v) for kind,v in integrated.items()},
                  common_source_subsets=[dict(common_source_set=list(subset),O=str(gk*v))
                                         for subset,v in common_groups.items()]))
    assert O_direct==O_columns==costs['all']
    assert costs['sigma=1']<=costs['sigma>=1/2']<=costs['all']
    assert M>0 and M<=F(3,4)*X
    assert Jaccard_O<=F(3,8)*X and tau_O<=F(3,4)*X
    assert Jaccard_O<=tau_O<=O_columns
    assert O_columns/M<=max_conditional['all']
    groups=[]
    for key,g in fine_groups.items():
        if not g['O']:continue
        assert g['gk']*g['column_integral']==g['O']
        assert sum(g['common_subsets'].values(),F(0))==g['O']
        groups.append(dict(k=g['k'],fine_source_set=g['fine_source_set'],fine_mass=str(g['fine_mass']),
                           z_volume=str(g['z_volume']),column_integral=str(g['column_integral']),
                           raw_column_max=str(g['max']),O=str(g['O']),sigma1_O=str(g['sigma1_O']),
                           common_source_subsets=[dict(common_source_set=list(subset),O=str(v))
                                                  for subset,v in g['common_subsets'].items()]))
    assert sum((F(g['O']) for g in groups),F(0))==O_columns
    out=dict(batch=batch,label=label,D=saved['D'],N=len(atoms),alpha=str(alpha),X=str(X),
       Evolume=str(Evolume),O=str(O_columns),O_over_X=str(O_columns/X),fine_record_mass_M=str(M),
       Jaccard_O=str(Jaccard_O),Jaccard_O_over_X=str(Jaccard_O/X),
       tau_O=str(tau_O),tau_O_over_X=str(tau_O/X),
       first_union_same_color_probability_checks=first_union_checks,
       sigma1_Jaccard_equals_tau_exact=True,
       Jaccard_O_le_3X_over_8=True,tau_O_le_3X_over_4=True,
       M_over_X=str(M/X),O_over_M_weighted_mean_conditional_column=str(O_columns/M),
       maximum_raw_column={kind:str(v) for kind,v in max_column.items()},
       maximum_conditional_column={kind:str(v) for kind,v in max_conditional.items()},
       max_raw_column_witness=witnesses,component_O={kind:str(v) for kind,v in costs.items()},
       positive_columns=len(columns),active_edges=active_edges,breakpoint_checks=breakpoint_checks,
       independent_rectangle_equals_column_trapezoid=True,original_O_equals_weighted_columns=True,
       fine_labels_at_same_z_are_kept_separate=True,positive_column_records=columns,fine_source_groups=groups,
       atoms=[dict(location=str(y),mass=str(w)) for y,w in atoms],elapsed_seconds=time.monotonic()-began)
    if gram_local is not None:
        out.update(gram_k=params['gram_k'],gram_radius=str(params['radius']),gram_local=gram_local,
                   local_point_fee=str(params['local_point_fee']))
    return out


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--repo-root',type=Path,default=Path(__file__).resolve().parents[2])
    p.add_argument('--output',type=Path,required=True);args=p.parse_args()
    flow=module(args.repo_root/'rounds/047/flow_probe.py','flow56')
    rec=module(args.repo_root/'rounds/046/threshold_probe.py','record56')
    prefix=module(args.repo_root/'rounds/046/prefix_probe.py','prefix56')
    glob=module(args.repo_root/'rounds/049/global_probe.py','global56')
    lag=module(args.repo_root/'rounds/054/lag_probe.py','lag56')
    obs=module(args.repo_root/'rounds/053/obstruction_probe.py','obs56')
    gram=module(args.repo_root/'rounds/055/gram_certificate.py','gram56')
    old=rec.load(args.repo_root,args.output)
    refs={c['label']:c for c in json.loads((args.repo_root/'rounds/055/verification.json').read_text())['shared_tests']['cases']}
    inputs=[c for c in flow.cases(old,prefix) if c[1] in lag.ORIGINAL|lag.EXTRA]
    inputs.append((0,'round53_four_atom_theta1',sorted(obs.ATOMS),F(1),{}))
    for k,batch in ((9,0),(12,1),(16,2),(24,2)):
        local=gram.case(k,batch);r=gram.radius(k);mass=F(23,100);tiny=F(3,2)*r
        inputs.append((batch,f'gram_three_atoms_k{k}',[(F(0),tiny),(F(7,20),mass),(F(10),1-mass-tiny)],
                       F(1),dict(gram_k=k,radius=r,local_point_fee=F(local['original_point_fee']))))
    began=time.monotonic();cases=[]
    for batch,label,atoms,alpha,params in inputs:
        c=evaluate(flow,rec,glob,old,batch,label,atoms,alpha,params)
        if label in refs:
            assert F(c['O'])==F(refs[label]['O']);c['round55_original_O_receipt']='passed'
        else:
            check=flow.kernel_probe(old,rec,atoms,alpha)
            assert F(check['positive'])==F(c['O']);c['independent_round47_full_E_kernel_O']='passed'
        cases.append(c)
        print(batch,label,'cmax',float(F(c['maximum_raw_column']['all'])),
              'conditional',float(F(c['maximum_conditional_column']['all'])),
              'O/X',float(F(c['O_over_X'])),'D',c['D'],flush=True)
    assert len(cases)==20
    out=dict(status='passed',cases=cases,batch_counts=[sum(c['batch']==b for c in cases) for b in range(3)],
             scope='full actual atomic E/P/common beta; fixed (z,k) columns never mix fine labels',
             arithmetic='Fraction exact',total_positive_columns=sum(c['positive_columns'] for c in cases),
             total_active_edges=sum(c['active_edges'] for c in cases),
             total_breakpoint_checks=sum(c['breakpoint_checks'] for c in cases),
             general_column_bound_proved=False,general_weak_type_bound_proved=False,
             elapsed_seconds=time.monotonic()-began,
             script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    args.output.write_text(json.dumps(out,indent=2)+'\n')
    print('PASSED',len(cases),out['batch_counts'],'positive columns',out['total_positive_columns'],
          'seconds',round(out['elapsed_seconds'],2),flush=True)


if __name__=='__main__':main()
