#!/usr/bin/env python3
"""Exact real-band high-occupancy transport truncations; finite n=1 samples."""
from pathlib import Path
from fractions import Fraction as F
from collections import defaultdict
import argparse, importlib.util, sys, json, time, hashlib
sys.dont_write_bytecode=True
ROOT=Path(__file__).resolve().parents[2]
THRESHOLDS=(F(1,2),F(1),F(2),F(4),F(8),F(16))

def archived43(output):
    """Import constructors with argv side effects contained and no old main run."""
    previous=list(sys.argv)
    try:
        sys.argv=['global_probe','--output',str(output)]
        spec=importlib.util.spec_from_file_location('archived_stage43',ROOT/'rounds/043/global_probe.py')
        mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
        return mod
    finally:sys.argv=previous

def below_length(a,b,va,vb,t):
    """Length where the affine interpolation is STRICTLY below t."""
    if min(va,vb)>=t:return F(0)
    if max(va,vb)<=t:return b-a
    root=a+(b-a)*(t-va)/(vb-va)
    assert a<root<b
    return root-a if va<t else b-root

def cases(old):
    labels={'separated_equal_4','separated_binary','two_near_half','old41_micro_k10','old41_micro_k12',
       'old41_chain_L12','old41_chain_L24','old41_chain_L48',
       'big_contact_L6_b8_s-1_side1','band_edge_sliver_L6_b10_side1',
       'random_two_clouds_s420701_L6','random_two_clouds_s420720_L10','random_two_clouds_s420739_L14',
       'asymmetric_tree_d3_s420924','asymmetric_tree_d4_s420947','asymmetric_tree_d5_s421154'}
    for batch,label,atoms,alpha,params in old.cases():
        if label in labels:yield batch,label,atoms,alpha,params
    for L in (12,24):
        yield 1,f'enhanced_chain_L{L}',old.old.old.chain(L,den=16,position=F(9,8),weight=F(7,4)),F(1),dict(L=L)
    yield 2,'old41_chain_L96',old.old.old.chain(96,den=4,position=F(9,8),weight=F(27,16)),F(1),dict(L=96)

def evaluate(old,batch,label,atoms,alpha,params):
    began=time.monotonic();saved,_,radii,_=old.observer(atoms,alpha,'band',4)
    cells=[dict(c,radius=str(radii[c['K']])) for c in saved['observer_cells']]
    hs=old.build_h(cells,alpha/8,'eligible');A=old.add(hs.values());av={y:A.value(y) for y,w in atoms}
    all_h_knots={z for h in hs.values() for z in h.knots}
    prefixes={k:old.add(hs[j] for j in hs if j<k) for k in hs}
    prefix_source={y:{k:s.value(y) for k,s in prefixes.items()} for y,w in atoms}
    prefix_difference_max=F(0);prefix_witness=None;edge_endpoint_count=0
    X=alpha*F(saved['eligible_volume']);M=sum((w*av[y] for y,w in atoms),F(0));Q=sum((w*av[y]**2 for y,w in atoms),F(0))
    assert X>0 and X/2<=M<=X
    raw=defaultdict(F);low=defaultdict(F);obsmean=F(0);obssquare=F(0)
    stop_levels=(F(0),F(1),F(2),F(4));obs_stop={u:F(0) for u in stop_levels}
    tab={t:dict(rho_max=F(0),rho_mass=F(0),rho_square_mass=F(0),rho_eq_one_mass=F(0),fixed_low_mass=F(0),rho_witness=None) for t in THRESHOLDS}
    for c in cells:
        k=c['K'];l=F(c['lo']);r=F(c['hi']);rad=radii[k];vk=2*rad
        cuts=sorted({l,r}|{z for z in all_h_knots if l<z<r}|{z for y,w in atoms for z in (y-rad,y+rad) if l<z<r})
        for a,b in zip(cuts,cuts[1:]):
            mid=(a+b)/2;present=[(y,w) for y,w in atoms if abs(y-mid)<rad];mass=sum((w for y,w in present),F(0));g=mass/vk
            assert alpha/2<g<=alpha
            va,vb=A.value(a),A.value(b);width=b-a;obsmean+=width*(va+vb)*g/2
            obssquare+=width*(va*va+va*vb+vb*vb)*g/3
            for u in stop_levels:obs_stop[u]+=g*below_length(a,b,-va,-vb,-u)
            for y,w in present:
                sy=prefix_source[y][k]
                for endpoint in (a,b):
                    difference=abs(sy-prefixes[k].value(endpoint));edge_endpoint_count+=1
                    assert difference<=F(1,2)
                    if difference>prefix_difference_max:
                        prefix_difference_max=difference;prefix_witness=dict(y=str(y),z=str(endpoint),K=k,source_prefix=str(sy),observer_prefix=str(prefixes[k].value(endpoint)))
                raw[y]+=width/vk
                low[y]+=below_length(a,b,va,vb,av[y]/2)/vk
            for t in THRESHOLDS:
                gm=sum((w for y,w in present if av[y]>t),F(0))/vk;rho=gm/g;z=tab[t]
                z['rho_mass']+=width*gm;z['rho_square_mass']+=width*g*rho*rho
                z['fixed_low_mass']+=below_length(a,b,va,vb,t/2)*gm
                if rho==1:z['rho_eq_one_mass']+=width*g
                if rho>z['rho_max']:z['rho_max']=rho;z['rho_witness']=dict(lo=str(a),hi=str(b),K=k,rho=str(rho),g=str(g),high_ball_mass=str(gm*vk))
    for y,w in atoms:
        assert raw[y]==av[y] and 0<=low[y]<=av[y]
        if av[y]:assert low[y]/av[y]<=min(F(1),F(1,2)+F(3,2)/av[y])
        direct=sum((max(F(0),min(F(c['hi']),y+F(c['radius']))-max(F(c['lo']),y-F(c['radius'])))/(2*F(c['radius'])) for c in cells),F(0))
        assert direct==av[y]
    assert obsmean<=X
    thresholds=[]
    for t in THRESHOLDS:
        high=[(y,w) for y,w in atoms if av[y]>t];qt=sum((w*av[y]**2 for y,w in high),F(0));mt=sum((w*av[y] for y,w in high),F(0))
        H=sum((w*av[y]*low[y] for y,w in high),F(0));mlow=sum((w*low[y] for y,w in high),F(0))
        z=tab[t];assert z['rho_mass']==mt and 0<=H<=qt and 0<=mlow<=mt
        assert H<=(F(1,2)+F(3,2)/t)*qt
        # This is a proved truncation budget, checked only as an implementation identity.
        assert Q<=(t+2)*X+H
        ps=[low[y]/av[y] for y,w in high]
        thresholds.append(dict(t=str(t),high_source_count=len(high),Qtail=str(qt),mu_tail_mass=str(mt),H=str(H),low_mu_tail_mass=str(mlow),
          H_over_Qtail=str(H/qt) if qt else None,H_over_X=str(H/X),Qtail_over_X=str(qt/X),t_Qtail_over_X=str(t*qt/X),mu_tail_over_X=str(mt/X),n1_contraction_factor=str(F(1,2)+F(3,2)/t),
          high_source_low_observer_probability=str(mlow/mt) if mt else None,
          min_source_low_probability=str(min(ps)) if ps else None,max_source_low_probability=str(max(ps)) if ps else None,
          rho_sup=str(z['rho_max']),rho_G_mean=str(mt/M),rho_G_second_moment=str(z['rho_square_mass']/M),rho_eq_one_G_mass=str(z['rho_eq_one_mass']),
          fixed_low_observer_G_high_mass=str(z['fixed_low_mass']),fixed_low_given_high=str(z['fixed_low_mass']/mt) if mt else None,rho_supremum_witness=z['rho_witness']))
    shifted={q:sum((w*max(F(0),av[y]-F(3,2))**q for y,w in atoms),F(0)) for q in (2,3)}
    assert shifted[2]<=2*obsmean and shifted[3]<=3*obssquare
    sharp={q:sum((w*max(F(0),av[y]-F(1,2))**q for y,w in atoms),F(0)) for q in (2,3)}
    diagonal=sum((w*sum((h.value(y)**2 for h in hs.values()),F(0)) for y,w in atoms),F(0))
    assert sharp[2]<=2*obsmean+diagonal
    assert sharp[3]<=3*obssquare+3*obsmean+M
    stop_loss=[]
    for u in stop_levels:
        value=sum((w*max(F(0),av[y]-F(3,2)-u) for y,w in atoms),F(0))
        assert value<=obs_stop[u]
        stop_loss.append(dict(u=str(u),source_loss=str(value),observer_tail_mass=str(obs_stop[u])))
    return dict(batch=batch,label=label,parameters=params,N=len(atoms),D=saved['D'],alpha=str(alpha),X=str(X),M=str(M),Q=str(Q),Q_over_X=str(Q/X),observer_A_mean_pairing=str(obsmean),
      max_source_A=str(max(av.values())),rank_budget=dict(c='3/2',moment_shift='1/2',sharp_shifted_square=str(sharp[2]),sharp_shifted_cube=str(sharp[3]),shifted_square=str(shifted[2]),shifted_cube=str(shifted[3]),observer_first=str(obsmean),observer_second=str(obssquare),stop_loss=stop_loss),thresholds=thresholds,prefix_Lipschitz_audit=dict(status='passed',bound='1/2',max_difference=str(prefix_difference_max),edge_endpoint_count=edge_endpoint_count,witness=prefix_witness),source=[dict(location=str(y),mass=str(w),A=str(av[y]),low_kernel_mass=str(low[y]),conditional_low_probability=str(low[y]/av[y]) if av[y] else None) for y,w in atoms],
      observer_cells=cells,original_intervals=saved['original_intervals'],original_volume=saved['original_volume'],J1_volume=saved['J1_volume'],exact_error='0',elapsed_seconds=time.monotonic()-began)

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--output',type=Path,required=True);args=parser.parse_args();old=archived43(args.output)
    start=time.monotonic();out=dict(status='running',arithmetic='Fraction exact',scope='n=1 full original band/J/K; true mu; strict A>t and A(y)>2A(z)',cases=[])
    for batch,label,atoms,alpha,params in cases(old):
        c=evaluate(old,batch,label,atoms,alpha,params);out['cases'].append(c);args.output.write_text(json.dumps(out)+'\n')
        print(batch,label,'maxA',float(F(c['max_source_A'])),'H/Qtail',[(r['t'],round(float(F(r['H_over_Qtail'])),9) if r['H_over_Qtail'] else None) for r in c['thresholds']],'sec',round(c['elapsed_seconds'],2),flush=True)
    out['batch_maxima']=[]
    for batch in sorted({c['batch'] for c in out['cases']}):
        group=[c for c in out['cases'] if c['batch']==batch];maxima=[]
        for t in THRESHOLDS:
            active=[(c,r) for c in group for r in c['thresholds'] if F(r['t'])==t and r['H_over_Qtail'] is not None]
            maximum={}
            for key in ('H_over_Qtail','H_over_X','Qtail_over_X','t_Qtail_over_X','high_source_low_observer_probability','max_source_low_probability','rho_sup'):
                if active:
                    c,r=max(active,key=lambda v:F(v[1][key]));maximum[key]=dict(label=c['label'],exact=r[key],display=float(F(r[key])))
                else:maximum[key]=None
            maxima.append(dict(t=str(t),active_count=len(active),maximum=maximum))
        out['batch_maxima'].append(dict(batch=batch,count=len(group),thresholds=maxima))
    out.update(status='passed',case_count=len(out['cases']),elapsed_seconds=time.monotonic()-start,exact_error='0',sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),finite_tests_prove_uniform_contraction=False)
    most=max(out['cases'],key=lambda c:F(c['prefix_Lipschitz_audit']['max_difference']))
    out['prefix_audit_summary']=dict(status='passed',bound='1/2',edge_endpoint_count=sum(c['prefix_Lipschitz_audit']['edge_endpoint_count'] for c in out['cases']),maximum_label=most['label'],maximum=most['prefix_Lipschitz_audit']['max_difference'],all_source_low_probability_bounds_checked=True,all_H_contraction_bounds_checked=True)
    args.output.write_text(json.dumps(out,indent=2)+'\n');print('PASSED',out['case_count'],'seconds',out['elapsed_seconds'],flush=True)
if __name__=='__main__':main()
