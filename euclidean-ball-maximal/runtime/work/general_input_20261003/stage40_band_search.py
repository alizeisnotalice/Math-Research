"""Bounded new actual band searches, with exact spatial/source and E1 ledgers.
Original band={alpha<M<=2alpha}, common-c stops, and eligible J>=2 are
implemented by the unchanged stage24 kernel. Eta is substituted for Pi only
after the whole original h support saturation certificate has passed.
"""
from fractions import Fraction as F
from pathlib import Path
from collections import defaultdict
import argparse,json,random,time,sys
sys.dont_write_bytecode=True
from stage39_excess_search import evaluate,support_saturation
from stage30_transition_weight_probe import build_h,Linear

ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/'output/general_input_20261003/stage40_band_search.json'
SEEDS={1:400711,2:400727,3:400743}

def chain(L,den=16,mirror=False,alpha=F(1),q=8,weight=F(9,8),position=F(3,2),jitter=False,seed=0,central_index=None):
    rng=random.Random(seed);scale=F(1,8) if mirror else F(1,4)
    rs=[scale/F(q**i) for i in range(L)];atoms=[]
    for i,r in enumerate(rs):
        w=weight*r*(1+F(rng.randint(-2,2),128) if jitter else 1)
        y=position*r/alpha
        atoms.append((y,w))
        if mirror:atoms.append((-y,w))
    p=rs[-1 if central_index is None else central_index]/den
    atoms.append((F(0),p));total=sum(w for x,w in atoms)
    assert total<1
    atoms.append((F(10)/alpha,1-total));atoms.sort()
    assert len({x for x,w in atoms})==len(atoms)
    return atoms

def cases(rnd):
    if rnd==1:
        specs=[dict(L=L,den=1000,pilot=True) for L in (4,16,40,64)]
        specs += [dict(L=L,den=16) for L in (4,16,24,32,40)]
    elif rnd==2:
        specs=[dict(L=L,den=den,mirror=True,alpha=alpha) for L,den,alpha in
               [(8,16,F(1)),(12,16,F(3,2)),(16,16,F(3,4)),(24,16,F(1)),
                (24,4,F(3,2)),(24,16,F(1)),(32,16,F(1)),(24,16,F(5,4))]]
        specs[5]['central_index']=8;specs[7]['jitter']=True
    else:
        specs=[dict(L=L,den=den,mirror=mirror,alpha=alpha,weight=F(27,16),position=F(9,8))
            for L,den,mirror,alpha in [(4,16,False,F(1)),(6,16,False,F(1)),
              (8,16,False,F(3,2)),(12,16,False,F(3,4)),(16,16,False,F(1)),
              (6,16,True,F(1)),(8,16,True,F(5,4)),(12,4,True,F(1)),
              (12,16,True,F(3,2)),(8,16,False,F(5,4)),(12,16,False,F(1))]]
        specs[-2]['jitter']=True;specs[-1]['q']=6
    for i,spec in enumerate(specs):
        pilot=spec.pop('pilot',False);spec['seed']=SEEDS[rnd]+i
        atoms=chain(**spec)
        yield i,spec,pilot,atoms

def augment(summary,complete):
    source=complete['complete_source_excess_ledger'];knots=complete['A_hfold']
    source_max=max(source,key=lambda v:F(v['A_exact']))
    space_max=max(knots,key=lambda v:F(v['value'])) if knots else None
    summary.update(max_source_A_exact=source_max['A_exact'],max_source_A_location=source_max['location'],
        center_source_A_exact=next(v['A_exact'] for v in source if v['location']=='0'),
        spatial_A_gt_one=F(summary['max_A_exact'])>1,
        source_A_gt_one=F(source_max['A_exact'])>1,
        maximum_spatial_A_location=space_max['x'] if space_max else None,
        E1_sign=None if summary['E1_exact'] is None else (F(summary['E1_exact'])>0)-(F(summary['E1_exact'])<0),
        E1_over_X_display=float(F(summary['E1_over_X_exact'])) if summary['E1_over_X_exact'] else None)
    stop=complete['complete_actual_matched_height_audit']['complete_actual_observer_stop']
    summary['actual_J_values']=sorted({v['J'] for v in stop['observer_cells']})
    summary['actual_K_values']=sorted({v['K'] for v in stop['observer_cells']})
    summary['original_band_interval_count']=len(stop['original_intervals'])
    complete['summary']=summary

def statistics(cases):
    ss=[c['summary'] for c in cases]
    best=max(ss,key=lambda s:F(s['E1_over_X_exact']) if s['E1_over_X_exact'] else -1)
    high=max(ss,key=lambda s:F(s['max_A_exact']))
    return dict(case_count=len(ss),spatial_A_gt_one_count=sum(s['spatial_A_gt_one'] for s in ss),
        source_A_gt_one_count=sum(s['source_A_gt_one'] for s in ss),
        saturated_count=sum(s['all_h_supports_saturated'] for s in ss),
        certified_E1_positive_count=sum(s['E1_sign']==1 for s in ss),
        certified_E1_negative_count=sum(s['E1_sign']==-1 for s in ss),
        certified_E1_zero_count=sum(s['E1_sign']==0 for s in ss),
        source_screen_only_count=sum(s['source_screen_only'] for s in ss),
        threshold_reached_count=sum(s['threshold_one_sixteenth_reached'] for s in ss),
        maximum_E1_over_X_case=best,maximum_spatial_A_case=high)

def verify(out):
    count=0
    for rnd in out['rounds']:
        for case in rnd['cases']:
            atoms=[(F(v['location']),F(v['mass'])) for v in case['source']]
            assert sum(w for x,w in atoms)==1
            s=case['summary'];assert s['scope']=='band'
            assert s['E1_certified'] and s['all_h_supports_saturated']
            assert F(s['minimum_u_over_t_exact'])>=1
            E=F(s['E1_exact']);X=F(s['X_exact'])
            assert E==F(s['T1_source_screen_exact'])-F(s['returned_excess_exact'])
            assert E/X==F(s['E1_over_X_exact']) and F(s['S_exact'])-E==F(s['base_exact'])<=1
            assert s['threshold_one_sixteenth_reached']==(E/X>=F(1,16));count+=1
        assert rnd['statistics']==statistics(rnd['cases'])
    for c in out['retained_complete_scopes'].values():
        s=c['summary'];row=c['complete_actual_matched_height_audit'];stop=row['complete_actual_observer_stop']
        atoms=[(F(v['location']),F(v['mass'])) for v in c['original_source']]
        hs=build_h(row['complete_original_cells'],F(stop['entrance_a']),'eligible');ev=defaultdict(F)
        for h in hs.values():
            for x,d in h.events.items():ev[x]+=d
        A=Linear(ev);assert A.save()==c['A_hfold']
        eta,saturation=support_saturation(atoms,stop,hs)
        assert saturation==c['complete_original_h_support_saturation']
        T=sum(w*max(F(0),A.value(x)-1) for x,w in atoms);ret=F(0)
        for v in c['complete_positive_clipping_pieces']:
            lo,hi,sl,bi=(F(v[k]) for k in ('lo','hi','B_slope','B_intercept'))
            pl,ph=eta.prefix(lo),eta.prefix(hi)
            assert A.value(lo)-1==sl*lo+bi and A.value(hi)-1==sl*hi+bi
            val=sl*(ph[1]-pl[1])+bi*(ph[0]-pl[0]);assert val==F(v['returned_excess_exact']);ret+=val
        assert T==F(s['T1_source_screen_exact']) and ret==F(s['returned_excess_exact']) and T-ret==F(s['E1_exact'])
    return dict(status='passed',all_saved_summary_count=count,complete_scope_count=len(out['retained_complete_scopes']),
        original_cells_hfold_support_minima_and_exact_clipping_recomputed=True,stopping_reexecuted=False)

def main():
    p=argparse.ArgumentParser();p.add_argument('--verify-saved',action='store_true');p.add_argument('--round',type=int,choices=(1,2,3));args=p.parse_args()
    if args.verify_saved:
        out=json.loads(OUT.read_text());out['verification']=verify(out);OUT.write_text(json.dumps(out)+'\n');print(out['verification']);return
    out=json.loads(OUT.read_text()) if OUT.exists() else dict(status='running',dimension=1,scope='band',rounds=[],retained_complete_scopes={},
        method='Exact original actual n1 band; all source inputs and summaries retained. Full original stop/cells/eta/hfold/saturation/clip ledgers retained for first counterexample and each round maxima. Screening is not an upper-bound proof.',
        seeds=SEEDS)
    done={r['round'] for r in out['rounds']}
    for rnd in ([args.round] if args.round else (1,2,3)):
        if rnd in done:continue
        t=time.monotonic();rd=dict(round=rnd,seed=SEEDS[rnd],cases=[]);out['rounds'].append(rd);best=None;high=None
        for i,spec,pilot,atoms in cases(rnd):
            label=f'stage40_r{rnd}_i{i}_L{spec["L"]}_band'
            summary,complete=evaluate(atoms,spec.get('alpha',F(1)),'band',label)
            complete['complete_actual_matched_height_audit']['label']=label
            augment(summary,complete)
            case=dict(index=i,parameters={k:str(v) if isinstance(v,F) else v for k,v in spec.items()},
                source=[dict(location=str(x),mass=str(w)) for x,w in atoms],summary=summary)
            rd['cases'].append(case)
            if rnd==3 and i==0:
                out['simplest_optimized_positive_case']=label;out['retained_complete_scopes'][label]=complete
            if summary['spatial_A_gt_one'] and 'first_spatial_counterexample' not in out:
                out['first_spatial_counterexample']=label;out['retained_complete_scopes'][label]=complete
            if summary['E1_sign']==1 and 'first_positive_E1' not in out:
                out['first_positive_E1']=label;out['retained_complete_scopes'][label]=complete
            value=F(summary['E1_over_X_exact']) if summary['E1_over_X_exact'] else -1
            if best is None or value>best[0]:best=(value,complete)
            if high is None or F(summary['max_A_exact'])>high[0]:high=(F(summary['max_A_exact']),complete)
            rd['statistics']=statistics(rd['cases']);OUT.write_text(json.dumps(out)+'\n')
            print(label,'D',summary['actual_D'],'A0',float(F(summary['center_source_A_exact'])),'maxA',float(F(summary['max_A_exact'])),
                'E/X',summary['E1_over_X_display'],'sat',summary['all_h_supports_saturated'],flush=True)
        for _,c in (best,high):out['retained_complete_scopes'][c['summary']['label']]=c
        rd['elapsed_seconds']=time.monotonic()-t;rd['status']='passed';OUT.write_text(json.dumps(out)+'\n')
    if len(out['rounds'])==3:
        out['status']='passed';out['statistics']=statistics([c for r in out['rounds'] for c in r['cases']]);out['verification']=verify(out)
    OUT.write_text(json.dumps(out)+'\n');print('saved',OUT,flush=True)

if __name__=='__main__':main()
