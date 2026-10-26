#!/usr/bin/env python3
"""Self-contained archive audit: weighted A/B, certified order and all closed values."""
from pathlib import Path
from fractions import Fraction as Q
import sys,os,time,json,gzip
sys.dont_write_bytecode=True
P=Path(__file__).resolve().parent;ROUND=P.parent;ROOT=ROUND.parent
sys.path[:0]=[str(ROUND),str(ROOT/'lower')]
from certify_rectangular import field
from full_static import ENV,load,save,sha,rt,decode,encode,mul,addto,ZERO
def main():
    start=time.monotonic();reg=load(P/'registration.json');assert sha(P/'exact_cells.py')==reg['source_sha256']
    for x in reg['inputs']:
        original=Path(x['path']);fallbacks=[original,ROUND/original.name,ROOT/'lower_round025'/original.name,ROOT/'lower_round003'/original.name,ROOT/'lower'/original.name]
        assert any(f.exists() and sha(f)==x['sha256'] for f in fallbacks),x
    F,u,v,bounds,_,_=field();a=1+u**3/2;C=decode(load(ROOT/'lower_round003/a2_equals_a3_results.json')['coefficients'])
    forms=[(decode(x['constant']),x['slope']) for x in load(ROUND/'capture_forms.json')['forms']];fb=[bounds(tuple(c)) for c,m in forms]
    events=load(ROUND/'events.json')['events'];E=[decode(x['value']) for x in events]
    def cmp(x,y):
        z=[c-d for c,d in zip(x,y)]
        if not any(z):return 0
        lo,hi=bounds(tuple(z));assert lo>0 or hi<0,'unresolvedexactcomparison'
        return 1 if lo>0 else -1
    with gzip.open(P/'interface.json.gz','rt') as f:interface=json.load(f)
    assert sha(P/'interface.json.gz')==load(P/'state.json')['interface_sha256']
    with gzip.open(ROUND/'descriptors.json.gz','rt') as f:original_descriptors=json.load(f)['descriptors']
    descriptors=interface['descriptors'];assert len(descriptors)==len(original_descriptors)==91
    weights=[];fm={(tuple(c),m):i for i,(c,m) in enumerate(forms)};radii={r['N']:F(a=decode(r['radius'])) for r in load(ROOT/'lower_round025/tail_geometry.json')['families']};z=(1,0,-1)
    for row,old in zip(descriptors,original_descriptors):
        assert row['max_N_by_ell1']==old['max_N_by_ell1'] and row['physical_tail_weight']==old['physical_tail_weight'];weights.append(decode(row['physical_tail_weight']));pairs=[]
        for ell,N in enumerate(row['max_N_by_ell1'],1):
            if not N:continue
            r=radii[N]
            for j in range(-1,5):pairs.append([fm[(((j+ell-1)*a-r).a,z[(j+ell-1)%3])],fm[((j*a+r).a,z[j%3])]])
        assert row['raw_interval_form_pairs']==pairs
    inherited=load(P/'cell0000_inherited.json');assert inherited['status']=='INHERITED_exact_first_cell_no_union'
    assert decode(inherited['ratio_constant'])==C and decode(inherited['ratio_slope'])==decode(load(ROUND/'first_cell_certificate.json')['D_coefficients'])
    best=[ZERO]*25;best_event=0;max_noninitial=None;max_noninitial_event=None;closed=[];previous_right=None;next_cell=0;objects=0;actual_cells=0;pidrows=[];previous_receipt_sha=None;new_descs=0
    def audit_cell(cell):
        nonlocal best,best_event,max_noninitial,max_noninitial_event,previous_right,next_cell,objects
        index=cell['cell_index'];assert index==next_cell
        A=decode(cell['ratio_constant']);B=decode(cell['ratio_slope'])
        if index:
            rows=cell['descriptor_responses'];assert len(rows)==91;AA=[ZERO]*25;BB=[ZERO]*25
            for di,row in enumerate(rows):
                assert row['descriptor']==di;addto(AA,mul(weights[di],decode(row['length_constant'])));addto(BB,[v*row['length_slope'] for v in weights[di]]);objects+=1
            assert [v*Q(16,3) for v in AA]==A and [v*Q(16,3) for v in BB]==B
            q=Q(cell['ordering_only_interior_dyadic']);assert Q(events[index]['isolating_interval'][1])<q<Q(events[index+1]['isolating_interval'][0])
            lows=[l+m*q for (l,h),(c,m) in zip(fb,forms)];highs=[h+m*q for (l,h),(c,m) in zip(fb,forms)];order=sorted(range(len(forms)),key=lambda i:lows[i])
            assert all(highs[i]<lows[j] for i,j in zip(order,order[1:]));assert __import__('hashlib').sha256(json.dumps(order,separators=(',',':')).encode()).hexdigest()==cell['ordering_sha256']
        residuals=[]
        for event_id,key in ((index,'left_residual'),(index+1,'right_residual')):
            residual=[x-y for x,y in zip(A,C)];addto(residual,mul(B,E[event_id]));assert residual==decode(cell[key]);residuals.append(residual)
        if previous_right is not None:assert previous_right==residuals[0]
        else:closed.append(residuals[0])
        previous_right=residuals[1];closed.append(residuals[1]);next_cell+=1
        chosen=residuals[1] if cmp(residuals[1],residuals[0])>0 else residuals[0];assert chosen==decode(cell['cell_max_residual'])
    audit_cell(inherited)
    for rp in sorted((P/'batches').glob('batch_*.receipt.json')):
        receipt=load(rp);name=rp.name.replace('.receipt.json','.registration.json');br=load(rp.with_name(name))
        assert sha(rp.with_name(name))==receipt['batch_registration_sha256'] and br['source_sha256']==reg['source_sha256']
        assert receipt['first_cell']==next_cell and br['first_cell']==next_cell and br['last_cell']==receipt['last_requested_cell']
        assert br['last_cell']-br['first_cell']+1<=64
        if previous_receipt_sha is not None:assert br.get('previous_receipt_sha256')==previous_receipt_sha
        artifact=P/'batches'/receipt['cell_artifact'];assert sha(artifact)==receipt['cell_artifact_sha256']
        with gzip.open(artifact,'rt') as f:cells=json.load(f)['cells']
        assert [c['cell_index'] for c in cells]==receipt['completed_cell_indices']
        for cell in cells:audit_cell(cell)
        actual_cells+=receipt['new_cells_computed_this_process'];new_descs+=receipt['new_descriptor_dual_union_objects'];pidrows.append({'first_cell':br['first_cell'],'last_cell':br['last_cell'],'runtime':receipt['runtime']});previous_receipt_sha=sha(rp)
    assert next_cell==964 and len(closed)==965 and actual_cells==963 and objects==new_descs==87633
    negative=0;zero_events=[];positive=[]
    for event_id,residual in enumerate(closed):
        if not any(residual):zero_events.append(event_id);continue
        lo,hi=bounds(tuple(residual));assert lo>0 or hi<0
        if lo>0:positive.append(event_id)
        else:negative+=1
        if cmp(residual,best)>0:best=residual;best_event=event_id
        if event_id and (max_noninitial is None or cmp(residual,max_noninitial)>0):max_noninitial=residual;max_noninitial_event=event_id
    state=load(P/'state.json');assert best==decode(state['best_global_residual']) and state['next_cell']==964
    collision=[c+d for c,d in zip(C,closed[-1])];maxratio=[c+d for c,d in zip(C,best)];noninitialratio=[c+d for c,d in zip(C,max_noninitial)]
    result={'status':'COMPLETE_all964_exact_affine_cells_no_improvement' if not positive else 'POSITIVE_contact_requires_finitewitness_certification','closed_events':965,'open_cells':964,'inherited_C_D_cells':1,'new_cells_each_computed_once':actual_cells,'new_descriptor_dual_union_objects':objects,'all_963_adjacent_closed_equalities_exact':True,'all279_endpoint_outward_orders_audited_per_new_cell':True,'independent_saved_weighted_A_B_recomputed':objects,'zero_residual_event_indices':zero_events,'strict_negative_closed_events':negative,'strict_positive_closed_events':positive,'global_maximum_residual_coefficients':encode(best),'global_maximum_residual_interval':list(map(str,bounds(tuple(best)))),'global_maximum_ratio_coefficients':encode(maxratio),'global_maximum_ratio_interval':list(map(str,bounds(tuple(maxratio)))),'global_maximum_event':best_event,'maximum_equals_original_complete_C':not any(best),'strictly_below_C_for_every_0_less_epsilon_less_a':not positive and zero_events==[0],'best_noninitial_contact':{'event':max_noninitial_event,'epsilon_coefficients':events[max_noninitial_event]['value'],'residual_coefficients':encode(max_noninitial),'residual_interval':list(map(str,bounds(tuple(max_noninitial)))),'ratio_interval':list(map(str,bounds(tuple(noninitialratio)))),'ratio_display_only':[float(x) for x in bounds(tuple(noninitialratio))]},'collision_limit_ratio_interval':list(map(str,bounds(tuple(collision)))),'collision_limit_used_only_as_continuous_upper':True,'remaining_cells':0,'parameter_point_union_evaluations':0,'cell_calculation_runtime_seconds':sum(row['runtime']['wall_seconds'] for row in pidrows),'cell_calculation_peak_rss_bytes':max(row['runtime']['peak_rss_bytes'] for row in pidrows),'batch_runtimes':pidrows,'audit_runtime':rt(start),'source_sha256':sha(Path(__file__)),'producer_source_sha256':reg['source_sha256'],'scope':'fixedz=(1,0,-1),0<=epsilon<a1,registered1024lengthfinitefunctionalonly;notgeneralperiod3placements,Cn,orallscaleupper;mainproofL1collarreview'}
    save(P/'window_result.json',result);print({k:v for k,v in result.items() if k not in ('global_maximum_ratio_coefficients','global_maximum_ratio_interval','best_noninitial_contact','collision_limit_ratio_interval','batch_runtimes')},flush=True)
if __name__=='__main__':
    assert all(os.environ.get(k)=='1' for k in ENV) and os.getpriority(os.PRIO_PROCESS,0)>=10
    main()
