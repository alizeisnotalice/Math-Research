#!/usr/bin/env python3
"""Exact, bounded postprocessing of saved real source posteriors only."""
from collections import defaultdict
from fractions import Fraction as F
from pathlib import Path
import datetime, gzip, hashlib, json, time

BASE=Path(__file__).resolve().parent
PREFIX='winner_radial_transport_guard_20261007'
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def write_new(p,obj):
    with Path(p).open('x') as h:json.dump(obj,h,ensure_ascii=False,indent=2);h.write('\n')

def main():
    started=time.perf_counter();regpath=BASE/f'{PREFIX}_registration.json'
    reg=json.loads(regpath.read_text())
    for e in reg['source_files']:assert sha(e['path'])==e['sha256']
    loaded=[];eligible=defaultdict(list)
    for fi,e in enumerate(reg['profiles_in_original_order']):
        assert sha(e['path'])==e['sha256']
        data=json.load(gzip.open(e['path'],'rt'));loaded.append(data)
        case=data['case'];n=case['n'];rid=case['round']
        for ri,record in enumerate(data['records']):
            if record['R'] is None or F(record['M'])==0:continue
            R=F(record['R'])
            winner=[(ci,c) for ci,c in enumerate(record['candidate_L']) if F(c['L'])==R]
            assert len(winner)==1
            wi,wc=winner[0];mr=F(wc['mass']);assert mr>0
            for ci,c in enumerate(record['candidate_L']):
                if F(c['L'])>R and F(c['mass'])>mr:
                    eligible[rid].append({'profile_index':fi,'receiver_index':ri,'candidate_index':ci,'winner_candidate_index':wi,'n':n,'case_id':case['case_id']})
    round_summaries=[];allselected=[]
    for rid in reg['round_ids']:
        selected=[];seen=set();limit=reg['max_pairs_per_round']
        for ticket in eligible[rid]:
            if ticket['n'] not in seen and len(selected)<limit:
                selected.append(ticket);seen.add(ticket['n'])
        for ticket in eligible[rid]:
            if len(selected)>=limit:break
            if ticket not in selected:selected.append(ticket)
        rs={'round':rid,'eligible_pairs':len(eligible[rid]),'selected_count':len(selected),'empty_selection':not selected,'selection_tickets':selected,'pairs':[]}
        for ordinal,ticket in enumerate(selected):
            data=loaded[ticket['profile_index']];case=data['case'];record=data['records'][ticket['receiver_index']]
            wc=record['candidate_L'][ticket['winner_candidate_index']];lc=record['candidate_L'][ticket['candidate_index']]
            n=case['n'];b=F(case['b']);R=F(record['R']);L=F(lc['L']);mr=F(wc['mass']);ml=F(lc['mass']);M=F(record['M'])
            c=M*L**n/ml;assert c==F(lc['c_exp_g']) and mr==M*R**n and ml>mr and L>R
            atoms=[tuple(map(F,y)) for y in case['source']['atoms']];weights=list(map(F,case['source']['weights']))
            pr=list(map(F,wc['posterior']));pl=list(map(F,lc['posterior']));dist=list(map(F,record['D']))
            assert len(pr)==len(pl)==len(atoms)==len(dist) and sum(pr)==sum(pl)==1
            assert all(p>=0 for p in pr+pl)
            physical=sorted(set(atoms));ix={y:i for i,y in enumerate(physical)}
            pw=[F(0)]*len(physical);pR=[F(0)]*len(physical);pL=[F(0)]*len(physical);pd=[None]*len(physical)
            for y,w,u,v,d in zip(atoms,weights,pr,pl,dist):
                j=ix[y];pw[j]+=w;pR[j]+=u;pL[j]+=v
                if pd[j] is not None:assert pd[j]==d
                pd[j]=d
            assert sum(pw)==1 and sum(pR)==sum(pL)==1
            beta=mr/ml;assert beta==c*(R/L)**n and beta<1
            assert all(v==beta*u for u,v in zip(pR,pL) if u>0)
            coordinate=[];sumD=F(0);sumphase=F(0);phase_cells=0;cdf_cells=0
            for i in range(n):
                sig=defaultdict(F)
                for y,u,v in zip(physical,pR,pL):sig[y[i]]+=u-v
                coords=sorted(sig);cumulative=F(0);Di=F(0);cdf=[]
                for lo,hi in zip(coords,coords[1:]):
                    cumulative+=sig[lo];contrib=(hi-lo)*abs(cumulative);Di+=contrib
                    cdf.append({'lo':str(lo),'hi':str(hi),'cdf_signed_difference':str(cumulative),'absolute_integral':str(contrib)})
                assert sum(sig.values())==0
                active=[y[i] for y,v in zip(physical,pL) if v>0]
                assert max(active)-min(active)<=L<=b
                cuts=sorted({F(0),2*b}|{z%(2*b) for z in coords}|{(z-b)%(2*b) for z in coords})
                phases=[];Ji=F(0)
                for lo,hi in zip(cuts,cuts[1:]):
                    mid=(lo+hi)/2
                    signs=[1 if (y[i]-mid)%(2*b)<b else -1 for y in physical]
                    u=sum((p*s for p,s in zip(pR,signs)),F(0));v=sum((p*s for p,s in zip(pL,signs)),F(0))
                    contribution=(hi-lo)*abs(u-v);Ji+=contribution
                    phases.append({'lo':str(lo),'hi':str(hi),'midpoint':str(mid),'physical_atom_signs':signs,'piR_v':str(u),'piL_v':str(v),'absolute_integral':str(contribution)})
                assert Ji==4*Di
                coordinate.append({'i':i,'D_i':str(Di),'phase_absolute_integral':str(Ji),'percoordinate_exact_identity_pass':True,'cdf_intervals':cdf,'phase_intervals':phases})
                sumD+=Di;sumphase+=Ji;phase_cells+=len(phases);cdf_cells+=len(cdf)
            periodic_mean=sumphase/(2*b*n);rhs_identity=2*sumD/(n*b);assert periodic_mean==rhs_identity
            edistance=sum((v*max(d-R,F(0))/2 for v,d in zip(pL,pd)),F(0));assert sumD>=edistance
            attempts=[];found=None
            for k in range(1,33):
                s=R+(L-R)/2**k;tail_lower=1-c*(s/L)**n
                attempts.append({'k':k,'s':str(s),'tail_lower':str(tail_lower),'positive_and_legal':s<=L and tail_lower>0})
                if s<=L and tail_lower>0:
                    local=(s-R)*tail_lower/2
                    actual_inside=sum((v for v,d in zip(pL,pd) if d<=s),F(0))
                    cap=c*(s/L)**n
                    assert actual_inside<=cap and sumD>=local
                    found={'k':k,'s':str(s),'tail_lower':str(tail_lower),'actual_inside_posterior_mass':str(actual_inside),'radial_cap':str(cap),'radial_cap_residual':str(cap-actual_inside),'transport_lower':str(local),'transport_residual':str(sumD-local),'exact_pass':True}
                    break
            pair={'round':rid,'ordinal':ordinal,'selection_ticket':ticket,'source_profile':reg['profiles_in_original_order'][ticket['profile_index']],'case':case,'saved_receiver':{k:v for k,v in record.items() if k!='candidate_L'},'saved_winner_candidate':wc,'saved_outer_candidate':lc,'physical_atoms':[list(map(str,y)) for y in physical],'physical_weights':list(map(str,pw)),'physical_piR':list(map(str,pR)),'physical_piL':list(map(str,pL)),'physical_saved_D':list(map(str,pd)),'R':str(R),'L':str(L),'mR':str(mr),'mL':str(ml),'c_exp_g':str(c),'beta':str(beta),'sum_coordinate_CDF_distance':str(sumD),'periodic_global_mean_abs_difference':str(periodic_mean),'identity_RHS':str(rhs_identity),'exact_global_identity_pass':True,'expected_distance_infty_to_inner_cube':str(edistance),'distance_transport_residual':str(sumD-edistance),'exact_distance_transport_pass':True,'local_s_attempts':attempts,'local_s_not_found':found is None,'local_s_found':found,'coordinate_results':coordinate}
            path=BASE/f'{PREFIX}_r{rid}_pair{ordinal+1}.json';write_new(path,pair)
            ps={'ordinal':ordinal,'n':n,'case_id':case['case_id'],'R':str(R),'L':str(L),'mR':str(mr),'mL':str(ml),'D_sum':str(sumD),'global_periodic_mean':str(periodic_mean),'distance_expectation':str(edistance),'local_s_found_k':None if found is None else found['k'],'local_s_not_found':found is None,'coordinate_checks':n,'phase_intervals':phase_cells,'cdf_intervals':cdf_cells,'pair_path':str(path),'pair_sha256':sha(path)}
            rs['pairs'].append(ps);allselected.append(ps)
        round_summaries.append(rs)
        print(json.dumps({'round':rid,'eligible':rs['eligible_pairs'],'selected':[x['n'] for x in rs['pairs']],'notfound':sum(x['local_s_not_found'] for x in rs['pairs'])}),flush=True)
    output={'status':'passed_fixed_saved_posterior_postprocessing','registration_sha256':sha(regpath),'script_sha256':sha(Path(__file__)),'rounds':round_summaries,'counts':{'selected_pairs':len(allselected),'coordinate_identity_checks':sum(x['coordinate_checks'] for x in allselected),'global_identity_checks':len(allselected),'distance_transport_checks':len(allselected),'positive_local_transport_checks':sum(not x['local_s_not_found'] for x in allselected),'local_s_not_found':sum(x['local_s_not_found'] for x in allselected),'phase_intervals':sum(x['phase_intervals'] for x in allselected),'cdf_intervals':sum(x['cdf_intervals'] for x in allselected)},'elapsed_seconds':time.perf_counter()-started,'source_or_winner_oracle_rerun':False,'old_probe_imported':False,'arithmetic':'Exact Fraction only; actual full source atoms and saved complete posterior weights.','scope':reg['scope']}
    resultpath=BASE/f'{PREFIX}_results.json';write_new(resultpath,output)
    receipt={'status':'completed','timestamp_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'attempts':1,'exit_status':0,'source_files':reg['source_files'],'source_profiles':reg['profiles_in_original_order'],'registration':{'path':str(regpath),'sha256':sha(regpath)},'script':{'path':str(Path(__file__)),'sha256':sha(Path(__file__))},'results':{'path':str(resultpath),'sha256':sha(resultpath)},'selected_pair_outputs':[{'path':x['pair_path'],'sha256':x['pair_sha256']} for x in allselected],'counts':output['counts'],'source_or_winner_oracle_rerun':False,'receiver_Lebesgue_integration':False,'nearmax_K_or_actual_history_claim':False,'scope':reg['scope']}
    write_new(BASE/f'{PREFIX}_receipt.json',receipt)
    print(json.dumps({'status':output['status'],'counts':output['counts'],'seconds':output['elapsed_seconds']}),flush=True)

if __name__=='__main__':main()
