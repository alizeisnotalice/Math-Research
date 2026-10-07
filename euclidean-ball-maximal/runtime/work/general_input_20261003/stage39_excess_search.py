"""New actual n1 inputs searching E1/X >= 1/16, not frozen-scope reruns.
T1=< (A-1)+,P> is an upper-bound screen. A true E1 is reported only when
all original hk supports are certified u>=tk, or when A<=1 so E1=0.
All actual common-c stopping and original full/band masks are unchanged.
"""
from fractions import Fraction as F
from collections import defaultdict
from bisect import bisect_right
from pathlib import Path
import json,time,hashlib,random,sys,argparse
sys.dont_write_bytecode=True
from stage24_clock_replacement_probe import run_scope
from stage29_matched_height_probe import audit
from stage25_stopped_energy_probe import MomentMeasure
from stage30_transition_weight_probe import build_h,Linear

ROOT=Path(__file__).resolve().parents[2];DATA=ROOT/'output/general_input_20261003'
OUT=DATA/'stage39_excess_search.json';SEEDS={1:390711,2:390727,3:390743}

def support_saturation(atoms,stop,hs):
    P=MomentMeasure(atoms=atoms)
    eta=MomentMeasure(cells=[tuple(F(v[k]) for k in ('lo','hi','density')) for v in stop['complete_eta_density']])
    u=lambda x:eta.left_distance(x)-P.left_distance(x)
    du=lambda x:eta.prefix(x)[0]-P.prefix(x)[0]
    measures={y for y,w in atoms}|{x for l,r,d in eta.cells for x in (l,r)}
    alpha=F(stop['alpha']);a=F(stop['entrance_a']);layers=[];minimum=None;valid=True
    for k,h in sorted(hs.items()):
        t=alpha*(F(1,2**k)/(2*a))**2/16;cuts=sorted(measures|set(h.knots));pieces=[]
        for lo,hi in zip(cuts,cuts[1:]):
            mid=(lo+hi)/2;hc=h.value(lo);slope=h.slope(mid)
            if hc==slope==0:continue
            idx=bisect_right(eta.points,mid)-1
            density=eta.cells[idx][2] if idx>=0 and mid<eta.cells[idx][1] else F(0)
            q=[u(lo),du(lo),density/2];length=hi-lo
            assert q[0]+q[1]*length+q[2]*length*length==u(hi)
            candidates=[(q[0],F(0)),(u(hi),length)]
            if q[2] and 0<-q[1]/(2*q[2])<length:
                z=-q[1]/(2*q[2]);candidates.append((q[0]-q[1]*q[1]/(4*q[2]),z))
            value,pos=min(candidates);assert value>=0
            ratio=value/t;minimum=ratio if minimum is None else min(minimum,ratio)
            valid &= value>=t
            pieces.append(dict(lo=str(lo),hi=str(hi),u_coefficients=list(map(str,q)),
                minimum_u_exact=str(value),minimum_position=str(lo+pos),height=str(t),
                minimum_u_over_t_exact=str(ratio),all_minimum_candidates=[dict(value=str(v),position=str(lo+p)) for v,p in candidates]))
        layers.append(dict(K=k,height=str(t),hfold=h.save(),complete_support_pieces=pieces))
    return eta,dict(all_h_supports_saturated=valid,minimum_u_over_t_exact=str(minimum) if minimum is not None else None,layers=layers)

def evaluate(atoms,alpha,kind,label):
    stop=run_scope(atoms,alpha,kind,4);row=audit(atoms,stop,label)
    hs=build_h(row['complete_original_cells'],F(stop['entrance_a']),'eligible');events=defaultdict(F)
    for h in hs.values():
        for x,v in h.events.items():events[x]+=v
    A=Linear(events);eta,saturation=support_saturation(atoms,stop,hs)
    M=sum((w*A.value(y) for y,w in atoms),F(0));T1=F(0);sources=[]
    for y,w in atoms:
        value=A.value(y);contribution=w*max(F(0),value-1);T1+=contribution
        sources.append(dict(location=str(y),mass=str(w),A_exact=str(value),source_excess=str(contribution)))
    returned=eta_A=F(0);pieces=[]
    for lo,hi in zip(A.knots,A.knots[1:]):
        slope=A.slope((lo+hi)/2);intercept=A.value(lo)-slope*lo
        pl,ph=eta.prefix(lo),eta.prefix(hi);eta_A+=slope*(ph[1]-pl[1])+intercept*(ph[0]-pl[0])
        cuts=[lo,hi]
        if (A.value(lo)-1)*(A.value(hi)-1)<0:cuts=[lo,(1-intercept)/slope,hi]
        for l,r in zip(cuts,cuts[1:]):
            if A.value((l+r)/2)<=1:continue
            pl,pr=eta.prefix(l),eta.prefix(r)
            value=slope*(pr[1]-pl[1])+(intercept-1)*(pr[0]-pl[0]);assert value>=0
            returned+=value
            pieces.append(dict(lo=str(l),hi=str(r),B_slope=str(slope),B_intercept=str(intercept-1),
                eta_mass=str(pr[0]-pl[0]),eta_first_moment=str(pr[1]-pl[1]),returned_excess_exact=str(value)))
    X=alpha*F(row['eligible_volume']);S=F(row['separate_J1_and_eligible_integrals']['eligible_J_at_least_2']['matched_signed_residual'])
    maxA=max([F(0)]+[A.value(x) for x in A.knots])
    certified=saturation['all_h_supports_saturated'] or maxA<=1
    E=T1-returned if certified else None
    if saturation['all_h_supports_saturated']:assert M-eta_A==S
    if maxA<=1:assert T1==returned==0
    if E is not None:
        base=S-E;assert base<=1
    result=dict(label=label,scope=kind,N=len(atoms),actual_D=stop['D'],alpha=str(alpha),
        max_A_exact=str(maxA),M_exact=str(M),X_exact=str(X),S_exact=str(S),
        T1_source_screen_exact=str(T1),T1_over_X_exact=str(T1/X) if X else None,
        all_h_supports_saturated=saturation['all_h_supports_saturated'],minimum_u_over_t_exact=saturation['minimum_u_over_t_exact'],
        E1_certified=certified,E1_exact=str(E) if E is not None else None,
        E1_over_X_exact=str(E/X) if E is not None and X else None,
        returned_excess_exact=str(returned) if certified else None,
        provisional_eta_pairing_not_true_return_if_unsaturated=str(returned) if not certified else None,
        base_exact=str(S-E) if E is not None else None,
        threshold_one_sixteenth_reached=E is not None and X>0 and E/X>=F(1,16),
        source_screen_only=not certified)
    complete=dict(summary=result,original_source=[dict(location=str(x),mass=str(w)) for x,w in atoms],
        complete_actual_matched_height_audit=row,complete_original_h_support_saturation=saturation,
        A_hfold=A.save(),complete_source_excess_ledger=sources,complete_positive_clipping_pieces=pieces)
    return result,complete

def generate_round(rnd,best=None):
    rng=random.Random(SEEDS[rnd]);cases=[]
    count={1:384,2:320,3:384}[rnd]
    for i in range(count):
        alpha=rng.choice([F(3,4),F(1),F(5,4),F(3,2),F(2)])
        if rnd==1:
            m0=rng.choice([F(1,32),F(1,16),F(3,32),F(1,8)])*rng.choice([F(31,32),F(1),F(33,32)])
            total=rng.choice([F(3,8),F(7,16),F(15,32),F(49,100)])
            m1=(total-m0)*F(rng.randint(2,8),10);m2=total-m0-m1
            z1=F(rng.randint(75,300),1000);z2=z1+F(rng.randint(75,375),1000)
            atoms=[(F(0),m0),(z1/F(3,2),m1),(z2/F(3,2),m2),(F(10),1-total)]
        elif rnd==2:
            depth=3+i%6;m0=rng.choice([F(1,64),F(1,32),F(1,16),F(3,32)])
            total=rng.choice([F(3,8),F(7,16),F(15,32),F(49,100)])
            q=rng.choice([F(1,2),F(1,3),F(1,4),F(2,3)])
            raw=[q**k*rng.choice([F(7,8),F(1),F(9,8)]) for k in range(depth)]
            masses=[(total-m0)*w/sum(raw) for w in raw]
            scale=rng.choice([F(3,4),F(1),F(5,4),F(3,2)])
            atoms=[(F(0),m0)]+[(scale*w/alpha*(1+F((-1)**k,2**(7+i%6))),w) for k,w in enumerate(masses)]+[(F(10),1-total)]
            # Break accidental duplicates using a tiny exact offset.
            atoms=[(x+F(k,2**(16+i%6)) if k and k<len(atoms)-1 else x,w) for k,(x,w) in enumerate(atoms)]
        else:
            assert best is not None
            # Five exact pilot locations are reused from already saved pilot
            # full scopes; subsequent mutations target their clean dyadic best.
            anchor=[(F(0),F(1,8)),(F(1,8),F(1,8)),(F(1,4),F(1,8)),(F(10),F(5,8))]
            atoms=list(anchor);alpha=F(3,2)
            if i<5:
                pilot=[(F(3,20),F(1,4)),(F(1,6),F(1,4)),(F(1,8),F(1,4)),(F(7,48),F(1,4)),(F(3,20),F(13,50))]
                x,y=pilot[i];atoms=[(F(0),F(1,8)),(x,F(1,8)),(y,F(1,8)),(F(10),F(5,8))]
                cases.append(dict(round=rnd,index=i,N=4,alpha=str(alpha),pilot_saved_full_scope=f'/private/tmp/stage39_analytic_trial_{i}.json',
                    source=[dict(location=str(x),mass=str(w)) for x,w in atoms]));continue
            j=rng.randrange(len(atoms)-1);k=rng.randrange(len(atoms)-1)
            if i%3:
                delta=F(rng.randint(-12,12),2**(9+i%6))
                if j!=k and atoms[j][1]+delta>0 and atoms[k][1]-delta>0:
                    atoms[j]=(atoms[j][0],atoms[j][1]+delta);atoms[k]=(atoms[k][0],atoms[k][1]-delta)
            if j>0:
                atoms[j]=(atoms[j][0]+F(rng.randint(-16,16),2**(8+i%7)),atoms[j][1])
            if i%11==0:alpha*=rng.choice([F(31,32),F(33,32),F(15,16),F(17,16)])
            if i%13==0 and len(atoms)<10:
                j=rng.randrange(len(atoms)-1);x,w=atoms[j];part=w/F(16)
                atoms[j]=(x,w-part);atoms.append((x+F((-1)**i,2**(7+i%5)),part))
        atoms=sorted(atoms)
        if len({x for x,w in atoms})<len(atoms):continue
        assert sum((w for x,w in atoms),F(0))==1 and all(w>0 for x,w in atoms)
        cases.append(dict(round=rnd,index=i,N=len(atoms),alpha=str(alpha),
            source=[dict(location=str(x),mass=str(w)) for x,w in atoms]))
    return cases

def verify_saved(out):
    keys=[];scopes=[]
    for r in out['rounds']:
        local=[]
        for case in r['cases']:
            assert sum(F(v['mass']) for v in case['source'])==1
            key=json.dumps(dict(alpha=case['alpha'],source=case['source']),sort_keys=True);keys.append(key);local.append(key)
            for s in case['scopes']:
                assert s['E1_certified'] and s['all_h_supports_saturated'] and F(s['minimum_u_over_t_exact'])>=1
                X=F(s['X_exact']);E=F(s['E1_exact'])
                assert E==F(s['T1_source_screen_exact'])-F(s['returned_excess_exact'])
                assert F(s['E1_over_X_exact'])==E/X and F(s['base_exact'])==F(s['S_exact'])-E<=1
                assert s['threshold_one_sixteenth_reached']==(E/X>=F(1,16));scopes.append(s)
        r['distinct_source_alpha_input_count']=len(set(local))
        r['threshold_reached_scope_count']=sum(s['threshold_one_sixteenth_reached'] for c in r['cases'] for s in c['scopes'])
    chosen={c['summary']['label']:c for c in [out['first_certified_threshold_scope'],out['maximum_certified_E1_scope']]+[r['extremal_complete_scope'] for r in out['rounds']]}
    for complete in chosen.values():
        s=complete['summary'];row=complete['complete_actual_matched_height_audit'];stop=row['complete_actual_observer_stop']
        atoms=[(F(v['location']),F(v['mass'])) for v in complete['original_source']]
        hs=build_h(row['complete_original_cells'],F(stop['entrance_a']),'eligible');events=defaultdict(F)
        for h in hs.values():
            for x,j in h.events.items():events[x]+=j
        A=Linear(events);assert A.save()==complete['A_hfold']
        eta,saturation=support_saturation(atoms,stop,hs)
        assert saturation==complete['complete_original_h_support_saturation'] and saturation['all_h_supports_saturated']
        T=sum(w*max(F(0),A.value(x)-1) for x,w in atoms);returned=F(0)
        for piece in complete['complete_positive_clipping_pieces']:
            lo,hi,sl,bi=(F(piece[k]) for k in ('lo','hi','B_slope','B_intercept'));pl,ph=eta.prefix(lo),eta.prefix(hi)
            assert A.value(lo)-1==sl*lo+bi and A.value(hi)-1==sl*hi+bi
            value=sl*(ph[1]-pl[1])+bi*(ph[0]-pl[0]);assert value==F(piece['returned_excess_exact']);returned+=value
        assert T==F(s['T1_source_screen_exact']) and returned==F(s['returned_excess_exact']) and T-returned==F(s['E1_exact'])
    return dict(status='passed',executed_case_count=len(keys),distinct_new_source_alpha_input_count=len(set(keys)),
        neutral_mutation_duplicate_case_count=len(keys)-len(set(keys)),scope_count=len(scopes),
        all_2176_saved_fraction_summaries_checked=True,complete_extremal_scope_count_recomputed_without_stopping=len(chosen),
        complete_support_minima_and_clipping_pairings_recomputed=True,stopping_reexecuted=False)

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--round',type=int,choices=(1,2,3));parser.add_argument('--verify-saved',action='store_true');args=parser.parse_args()
    start=time.monotonic();out=json.loads(OUT.read_text()) if OUT.exists() else dict(status='running',rounds=[],
        original_common_c_stopping_source_and_full_band_masks_unchanged=True,frozen_42_scopes_reexecuted=False,
        T1_screen_not_identified_with_E1=True,bridge_eta_substitution_requires_exact_support_saturation=True,
        atomic_pressure_only=True,L1_counterexample_not_claimed=True,additive_C_target_not_refuted=True)
    if args.verify_saved:
        out['verification']=verify_saved(out)
        out['distinct_new_input_count']=out['verification']['distinct_new_source_alpha_input_count']
        out['new_input_count_semantics']='executed new cases; neutral mutations include repeated own-round input parameters, explicitly counted in verification'
        out['source_sha256']=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
        OUT.write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out['verification'],indent=2));return
    best=out.get('best_case_parameters');bestscore=F(out.get('maximum_certified_E1_over_X_exact','-1'))
    for rnd in ([args.round] if args.round else [1,2,3]):
        assert not any(r['round']==rnd for r in out['rounds']),'refusing to rerun completed new-input round'
        r=dict(round=rnd,seed=SEEDS[rnd],cases=[]);out['rounds'].append(r);rs=time.monotonic();roundbest=F(-1)
        for case in generate_round(rnd,best):
            atoms=[(F(v['location']),F(v['mass'])) for v in case['source']];case['scopes']=[]
            for kind in ('full','band'):
                label=f'stage39_r{rnd}_i{case["index"]}_N{case["N"]}_{kind}'
                if kind=='full' and case.get('pilot_saved_full_scope') and Path(case['pilot_saved_full_scope']).exists():
                    complete=json.loads(Path(case['pilot_saved_full_scope']).read_text());summary=complete['summary']
                    assert complete['original_source']==case['source'] and summary['alpha']==case['alpha']
                    summary['new_pilot_execution_reused_without_stopping_rerun']=True
                else:summary,complete=evaluate(atoms,F(case['alpha']),kind,label)
                case['scopes'].append(summary)
                score=F(summary['E1_over_X_exact']) if summary['E1_over_X_exact'] is not None else F(-1)
                if score>roundbest:r['extremal_complete_scope']=complete;roundbest=score
                if score>bestscore:
                    bestscore=score;best=dict(source=case['source'],alpha=case['alpha']);out['best_case_parameters']=best
                    out['maximum_certified_E1_over_X_exact']=str(score);out['maximum_certified_E1_scope']=complete
                if summary['threshold_one_sixteenth_reached'] and 'first_certified_threshold_scope' not in out:
                    out['first_certified_threshold_scope']=complete
                    print('FIRST_CERTIFIED_THRESHOLD',label,'E1/X',summary['E1_over_X_exact'],'minu/t',summary['minimum_u_over_t_exact'],flush=True)
                print(label,'A',float(F(summary['max_A_exact'])),'E1/X',float(score) if summary['E1_certified'] else 'T1-screen-only',
                    'sat',summary['all_h_supports_saturated'],'D',summary['actual_D'],flush=True)
            r['cases'].append(case);r['elapsed_seconds']=time.monotonic()-rs
            out['source_sha256']=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
            if len(r['cases'])%4==0 or any(s['threshold_one_sixteenth_reached'] for s in case['scopes']):OUT.write_text(json.dumps(out,indent=2)+'\n')
        r['maximum_certified_E1_over_X_exact']=str(roundbest)
        r['scope_count']=sum(len(c['scopes']) for c in r['cases']);OUT.write_text(json.dumps(out,indent=2)+'\n')
        print('ROUND_COMPLETE',rnd,'cases',len(r['cases']),'best',float(roundbest),'seconds',r['elapsed_seconds'],flush=True)
    scopes=[s for r in out['rounds'] for c in r['cases'] for s in c['scopes']]
    out.update(status='passed' if len(out['rounds'])==3 else 'partial_passed',new_input_count=sum(len(r['cases']) for r in out['rounds']),scope_count=len(scopes),
        threshold_reached_scope_count=sum(s['threshold_one_sixteenth_reached'] for s in scopes),
        uncertified_T1_screen_only_scope_count=sum(s['source_screen_only'] for s in scopes),elapsed_this_invocation=time.monotonic()-start)
    OUT.write_text(json.dumps(out,indent=2)+'\n');print('DONE',out['status'],'best',float(bestscore),'threshold_scopes',out['threshold_reached_scope_count'],flush=True)

if __name__=='__main__':main()
