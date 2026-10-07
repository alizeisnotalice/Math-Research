#!/usr/bin/env python3
"""Exact n1 true common-record conditional-source sigma diagnostics.
Global cases use full atomic E; the C-infinity extension is a local box certificate.
Fraction geometry/records/costs; Decimal entropy diagnostics only.
"""
from pathlib import Path
from fractions import Fraction as F
from decimal import Decimal as D, localcontext
from collections import defaultdict
import argparse, importlib.util, json, hashlib, time, sys
sys.dont_write_bytecode = True


def module(path,name):
    spec=importlib.util.spec_from_file_location(name,path)
    m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
    return m


def dec(q): return D(q.numerator)/D(q.denominator)


def band(q):
    if q <= F(1,4): return '0<q<=1/4'
    if q <= F(1,2): return '1/4<q<=1/2'
    if q <= F(3,4): return '1/2<q<=3/4'
    if q < 1: return '3/4<q<1'
    return 'q=1'


def outside_length_integral(xa,xb,za,zb,R):
    knots={xa,xb}
    for z in (za,zb):
        for sign in (-1,1):
            x=z+sign*R
            if xa < x < xb: knots.add(x)
    def value(x):
        return zb-za-max(F(0),min(zb,x+R)-max(za,x-R))
    knots=sorted(knots)
    return sum(((b-a)*(value(a)+value(b))/2 for a,b in zip(knots,knots[1:])),F(0))


def evaluate(flow,record,glob,shear,old,batch,label,atoms,alpha):
    began=time.monotonic()
    saved,radii,cells,_=flow.context(old,atoms,alpha)
    I=alpha/4;lo=alpha/4;hi=alpha/2
    Evolume=F(saved['eligible_volume']);X=alpha*Evolume
    data=[]
    for a,b,J,trace in cells:
        mid=(a+b)/2;rec=record.intervals(trace,lo,hi)
        members={k:frozenset(i for i,(y,w) in enumerate(atoms) if abs(y-mid)<radii[k])
                 for k in range(1,saved['D']+1)}
        if data and data[-1][1]==a and data[-1][2:]==(rec,members):
            data[-1]=(data[-1][0],b,rec,members)
        else:data.append((a,b,rec,members))
    assert sum((b-a for a,b,rec,members in data),F(0))==Evolume
    O=O_polygon=O_length=O_two_source=T=Ss=F(0)
    Es=D(0);rows=[];area_checks=0
    sigma_bins=defaultdict(F);theta_joint=defaultdict(lambda:defaultdict(F))
    groups={key:F(0) for key in ('sigma=1','sigma>=3/4','sigma>=1/2','sigma<=1/4')}
    group_gap={key:defaultdict(F) for key in groups}
    sigma_gap=defaultdict(lambda:defaultdict(F))
    O_gap=defaultdict(F);joint_sigma1_theta1=F(0)
    for xc,(xa,xb,ix,mx) in enumerate(data):
        for zc,(za,zb,iz,mz) in enumerate(data):
            rect_area=(xb-xa)*(zb-za)
            for j,(ja,jb) in ix.items():
                R=radii[j];Vj=2*R
                out=rect_area-glob.strip_area(xa,xb,za,zb,R)
                if not out:continue
                polygon_area=length_area=None
                for k,(ka,kb) in iz.items():
                    if k<=j:continue
                    W=min(jb,kb)-max(ja,ka)
                    if W<=0:continue
                    common_set=mx[j]&mz[k]
                    common=sum((atoms[i][1] for i in common_set),F(0))
                    if not common:continue
                    r=radii[k];Vk=2*r
                    fine_mass=sum((atoms[i][1] for i in mz[k]),F(0))
                    coarse_mass=sum((atoms[i][1] for i in mx[j]),F(0))
                    lost=sum((atoms[i][1] for i in mx[j]-mz[j]),F(0))
                    assert lost>0 and W*Vj<=lost
                    sigma=common/fine_mass;tau=common/coarse_mass;theta=W*Vj/lost
                    assert 0<sigma<=1 and 0<tau<=1 and 0<theta<=1
                    gk=fine_mass/Vk
                    assert alpha/4<gk<=alpha
                    fine_geom=glob.strip_area(xa,xb,za,zb,R+r)-glob.strip_area(xa,xb,za,zb,R)
                    assert fine_geom==out
                    if polygon_area is None:
                        rectangle=[(xa,za),(xb,za),(xb,zb),(xa,zb)]
                        polygon_area=sum((shear.area(shear.clip(rectangle,A,B,-R))
                                          for A,B in ((F(-1),F(1)),(F(1),F(-1)))),F(0))
                        length_area=outside_length_integral(xa,xb,za,zb,R)
                        assert polygon_area==length_area==out
                        area_checks+=1
                    factor=W*common/(I*Vj*Vk)
                    cost=factor*out
                    t=W*sigma*fine_geom/(I*Vj)
                    assert cost==gk*t and alpha*t/4<=cost<=alpha*t
                    O+=cost;O_polygon+=factor*polygon_area;O_length+=factor*length_area
                    O_two_source+=theta*common*lost*out/(I*Vj*Vj*Vk)
                    T+=t;Ss+=sigma*cost
                    if sigma!=1: Es+=-dec(sigma).ln()*dec(cost)
                    sb,tb=band(sigma),band(theta)
                    sigma_bins[sb]+=cost;theta_joint[sb][tb]+=cost
                    sigma_gap[sb][k-j]+=cost;O_gap[k-j]+=cost
                    conditions={'sigma=1':sigma==1,'sigma>=3/4':sigma>=F(3,4),
                                'sigma>=1/2':sigma>=F(1,2),'sigma<=1/4':sigma<=F(1,4)}
                    for key,condition in conditions.items():
                        if condition:groups[key]+=cost;group_gap[key][k-j]+=cost
                    if sigma==theta==1:joint_sigma1_theta1+=cost
                    rows.append(dict(xcell=xc,zcell=zc,x_interval=[str(xa),str(xb)],
                                     z_interval=[str(za),str(zb)],j=j,k=k,gap=k-j,
                                     W=str(W),lost=str(lost),fine_mass=str(fine_mass),
                                     coarse_mass=str(coarse_mass),common_mass=str(common),
                                     sigma=str(sigma),theta=str(theta),tau=str(tau),gk=str(gk),
                                     outside_area=str(out),O=str(cost),T_sigma=str(t)))
    assert O==O_polygon==O_length==O_two_source
    assert alpha*T/4<=O<=alpha*T and 0<=Ss<=O
    assert sum(sigma_bins.values(),F(0))==sum(O_gap.values(),F(0))==O
    return dict(batch=batch,label=label,alpha=str(alpha),D=saved['D'],Evolume=str(Evolume),X=str(X),
                O=str(O),O_over_X=str(O/X),T_sigma=str(T),T_sigma_over_E=str(T/Evolume),
                O_over_alpha_Tsigma=str(O/(alpha*T)) if T else None,
                sigma_min=str(min(F(r['sigma']) for r in rows)) if rows else None,
                sigma_max=str(max(F(r['sigma']) for r in rows)) if rows else None,
                tau_min=str(min(F(r['tau']) for r in rows)) if rows else None,
                tau_max=str(max(F(r['tau']) for r in rows)) if rows else None,
                sigma_cost={key:str(v) for key,v in sigma_bins.items()},
                sigma_requested_groups={key:str(v) for key,v in groups.items()},
                sigma_requested_groups_overlap=True,
                sigma_group_gap={key:{str(gap):str(v) for gap,v in value.items()} for key,value in group_gap.items()},
                sigma_gap={key:{str(gap):str(v) for gap,v in value.items()} for key,value in sigma_gap.items()},
                sigma_theta_joint={key:{tb:str(v) for tb,v in value.items()} for key,value in theta_joint.items()},
                O_gap={str(gap):str(v) for gap,v in O_gap.items()},
                sigma1_theta1_cost=str(joint_sigma1_theta1),
                E_sigma_endpoint_entropy=str(Es),S_sigma_square=str(Ss),
                positive_cells=len(rows),three_independent_area_checks=area_checks,
                fine_lag_O_equals_full_outer_O=True,two_source_O_identity=True,
                alpha_quarter_Tsigma_le_O_le_alpha_Tsigma=True,
                positive_cell_records=rows,elapsed_seconds=time.monotonic()-began)


def smooth_local_certificate(obs):
    baseline=obs.main();assert baseline['status']=='passed'
    tx,tz=obs.trace(obs.X),obs.trace(obs.Z)
    alpha=F(1);j,k=4,5;R,r=obs.radius(j),obs.radius(k)
    coarse={i for i,(y,w) in enumerate(obs.ATOMS) if abs(obs.X-y)<R}
    fine={i for i,(y,w) in enumerate(obs.ATOMS) if abs(obs.Z-y)<r}
    common=sum((obs.ATOMS[i][1] for i in coarse&fine),F(0))
    fine_mass=sum((obs.ATOMS[i][1] for i in fine),F(0))
    coarse_mass=sum((obs.ATOMS[i][1] for i in coarse),F(0))
    lost=sum((w for y,w in obs.ATOMS if abs(obs.X-y)<R and abs(obs.Z-y)>=R),F(0))
    W=F(baseline['W']);sigma=common/fine_mass;tau=common/coarse_mass
    theta=W*2*R/lost
    assert sigma==theta==1 and tau==F(20,31)
    rows=[]
    for c in baseline['smooth_cases']:
        eps=F(c['epsilon']);h=F(c['smoothing_support'])
        assert obs.X-obs.Z-2*eps>R and obs.X-obs.Z+2*eps<R+r
        # Baseline interval certificate preserves every source/radius member,
        # every relevant clock, and the full MP band on these positive boxes.
        area=4*eps*eps;O=W/(alpha/4)*common/(4*R*r)*area
        T=W/(alpha/4)*sigma/(2*R)*area
        assert O==F(c['positive_fee']) and O==F(3,5)*T
        rows.append(dict(batch=c['batch'],epsilon=str(eps),h=str(h),sigma='1',theta='1',tau='20/31',
                         local_box_O=str(O),local_box_Tsigma=str(T),E_sigma='0',S_sigma=str(O),
                         actual_smooth_band_and_record_interval_certificate='passed',
                         global_smooth_E_integral_computed=False))
    return rows


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--repo-root',type=Path,default=Path(__file__).resolve().parents[2])
    parser.add_argument('--output',type=Path,required=True);args=parser.parse_args()
    flow=module(args.repo_root/'rounds/047/flow_probe.py','flow55')
    record=module(args.repo_root/'rounds/046/threshold_probe.py','record55')
    prefix=module(args.repo_root/'rounds/046/prefix_probe.py','prefix55')
    glob=module(args.repo_root/'rounds/049/global_probe.py','global55')
    shear=module(args.repo_root/'rounds/048/shear_probe.py','shear55')
    lag=module(args.repo_root/'rounds/054/lag_probe.py','lag55')
    obs=module(args.repo_root/'rounds/053/obstruction_probe.py','obstruction55')
    old=record.load(args.repo_root,args.output)
    refs={c['label']:c for c in json.loads((args.repo_root/'rounds/054/verification.json').read_text())['lag_tests']['cases']}
    inputs=[c for c in flow.cases(old,prefix) if c[1] in lag.ORIGINAL|lag.EXTRA]
    inputs.append((0,'round53_four_atom_theta1',sorted(obs.ATOMS),F(1),dict(family='round53_atomic_ancestor_of_smooth_box')))
    started=time.monotonic();cases=[]
    with localcontext() as ctx:
        ctx.prec=80
        for batch,label,atoms,alpha,params in inputs:
            c=evaluate(flow,record,glob,shear,old,batch,label,atoms,alpha)
            if label in refs:
                assert F(c['O'])==F(refs[label]['O'])
                c['round54_actual_O_receipt']='passed'
            else:
                check=flow.kernel_probe(old,record,atoms,alpha)
                assert F(check['positive'])==F(c['O'])
                c['independent_round47_full_E_kernel_O']='passed'
            cases.append(c)
            O=F(c['O']);group=F(c['sigma_requested_groups']['sigma=1'])
            print(batch,label,'positive',c['positive_cells'],'sigma1/O',float(group/O) if O else None,
                  'O/(alphaT)',c['O_over_alpha_Tsigma'],flush=True)
    assert len(cases)==16
    out=dict(status='passed',cases=cases,batch_counts=[sum(c['batch']==b for c in cases) for b in range(3)],
             local_smooth_sigma1_theta1_certificates=smooth_local_certificate(obs),
             total_positive_cells=sum(c['positive_cells'] for c in cases),
             total_independent_area_checks=sum(c['three_independent_area_checks'] for c in cases),
             scope='global atomic full E and actual shared beta; local smooth boxes separately identified',
             arithmetic='Fraction; 80-digit Decimal entropy diagnostics',
             entropy_is_interval_certificate=False,general_sigma_budget_proved=False,
             general_weak_type_bound_proved=False,elapsed_seconds=time.monotonic()-started,
             script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    args.output.write_text(json.dumps(out,indent=2)+'\n')
    print('PASSED',len(cases),out['batch_counts'],'positive',out['total_positive_cells'],
          'seconds',round(out['elapsed_seconds'],2),flush=True)


if __name__=='__main__':main()
