#!/usr/bin/env python3
"""One exact affine coefficient computation per registered topology cell."""
from pathlib import Path
from fractions import Fraction as Q
import sys,os,time,json,gzip,hashlib,traceback
sys.dont_write_bytecode=True
P=Path(__file__).resolve().parent;ROUND=P.parent;ROOT=ROUND.parent
sys.path[:0]=[str(ROUND),str(ROOT/'lower')]
from certify_rectangular import field
from full_static import ENV,load,save,sha,rt,decode,encode,mul,addto,ZERO

def main(stage,arg=None):
    start=time.monotonic();reg=load(P/'registration.json');assert sha(Path(__file__))==reg['source_sha256']
    for row in reg['inputs']:assert sha(Path(row['path']))==row['sha256']
    F,u,v,bounds,_,_=field();C=decode(load(ROOT/'lower_round003/a2_equals_a3_results.json')['coefficients']);a=1+u**3/2
    scratch=Path(reg['scratch_cells']);formsfile=load(ROUND/'capture_forms.json');forms=[(decode(x['constant']),x['slope']) for x in formsfile['forms']]
    events=load(ROUND/'events.json')['events'];E=[decode(x['value']) for x in events];N=len(events)-1;assert N==964
    def bnd(z):return bounds(tuple(z))
    def cmp(x,y):
        delta=[a-b for a,b in zip(x,y)]
        if not any(delta):return 0
        lo,hi=bnd(delta)
        if lo>0:return 1
        if hi<0:return -1
        raise ArithmeticError('rootisolationdidnotseparateexactcoefficientcomparison')
    def endpoint_residual(A,B,e):
        z=[x-y for x,y in zip(A,C)];addto(z,mul(B,e));return z
    def certificate(index,A,B,metadata):
        left=endpoint_residual(A,B,E[index]);right=endpoint_residual(A,B,E[index+1]);choose_right=cmp(right,left)>0
        m=right if choose_right else left;event_id=index+1 if choose_right else index;ratio=[c+d for c,d in zip(C,m)]
        return dict(metadata,cell_index=index,left_event=index,right_event=index+1,ratio_constant=encode(A),ratio_slope=encode(B),left_residual=encode(left),right_residual=encode(right),left_residual_interval=list(map(str,bnd(left))),right_residual_interval=list(map(str,bnd(right))),cell_max_residual=encode(m),cell_max_residual_interval=list(map(str,bnd(m))),cell_max_ratio_interval=list(map(str,bnd(ratio))),cell_max_event=event_id,source_sha256=reg['source_sha256'])
    interfacepath=P/'interface.json.gz';statepath=P/'state.json'
    if stage=='init':
        assert not statepath.exists(),'PWAalreadyinitialized'
        radii={x['N']:F(a=decode(x['radius'])) for x in load(ROOT/'lower_round025/tail_geometry.json')['families']}
        formmap={(tuple(c),m):i for i,(c,m) in enumerate(forms)}
        zero_id=formmap[(tuple([ZERO]*25),0)];domain_id=formmap[((3*a).a,0)];z=(1,0,-1)
        with gzip.open(ROUND/'descriptors.json.gz','rt') as f:descriptors=json.load(f)['descriptors']
        interface=[]
        for row in descriptors:
            pairs=[]
            for ell,n in enumerate(row['max_N_by_ell1'],1):
                if not n:continue
                r=radii[n]
                for j in range(-1,5):
                    pairs.append([formmap[(((j+ell-1)*a-r).a,z[(j+ell-1)%3])],formmap[((j*a+r).a,z[j%3])]])
            interface.append({'max_N_by_ell1':row['max_N_by_ell1'],'physical_tail_weight':row['physical_tail_weight'],'raw_interval_form_pairs':pairs})
        data={'descriptors':interface,'zero_form':zero_id,'domain_form':domain_id,'source_sha256':reg['source_sha256'],'capture_forms_sha256':sha(ROUND/'capture_forms.json'),'events_sha256':sha(ROUND/'events.json'),'descriptors_sha256':sha(ROUND/'descriptors.json.gz')}
        with interfacepath.open('wb') as raw:
            with gzip.GzipFile(fileobj=raw,mode='wb',mtime=0) as out:out.write((json.dumps(data,separators=(',',':'))+'\n').encode())
        first=load(ROUND/'first_cell_certificate.json');D=decode(first['D_coefficients'])
        cell0=certificate(0,C,D,{'status':'INHERITED_exact_first_cell_no_union','descriptor_dual_unions':0,'source_certificate_sha256':sha(ROUND/'first_cell_certificate.json')})
        assert cell0['right_residual']==first['right_contact_minus_C_coefficients']
        save(scratch/'cell0000.json',cell0)
        state={'status':'INITIALIZED_inherited_first_cell','next_cell':1,'total_cells':N,'new_cells_computed':0,'inherited_cells':1,'best_global_residual':{},'best_global_event':0,'interface_sha256':sha(interfacepath),'source_sha256':reg['source_sha256']}
        save(statepath,state);print({'stage':'initialized','pid':os.getpid(),'inherited_cell':0,'new_cells':0,'next_cell':1,'interface_sha256':state['interface_sha256'],'runtime':rt(start)},flush=True);return
    assert stage=='batch' and arg
    brpath=P/'batches'/arg;br=load(brpath);state=load(statepath)
    assert br['source_sha256']==reg['source_sha256'] and br['expected_state_sha256']==sha(statepath)
    assert br['interface_sha256']==sha(interfacepath)==state['interface_sha256']
    first,last=br['first_cell'],br['last_cell'];assert first==state['next_cell'] and 1<=first<=last<N and last-first+1<=64
    if (P/'positive_candidate_LOCK.json').exists() and not br.get('explicit_inventory_after_positive',False):raise ValueError('positivecandidateawaitsmaincertificationpriority')
    with gzip.open(interfacepath,'rt') as f:interface=json.load(f)
    descriptors=interface['descriptors'];assert len(descriptors)==91
    weights=[decode(x['physical_tail_weight']) for x in descriptors];zero_id=interface['zero_form'];domain_id=interface['domain_form']
    const_bounds=[bnd(c) for c,m in forms];actual_cells=0;new_objects=0;completed=[];max_new=None;max_new_event=None;fault_cell=first
    print({'stage':'registered_exact_PWA_batch','pid':os.getpid(),'first_cell':first,'last_cell':last,'descriptors_per_new_cell':91,'remaining_before':N-first,'runtime':rt(start)},flush=True)
    try:
        for index in range(first,last+1):
            fault_cell=index;outpath=scratch/f'cell{index:04d}.json';inflight=scratch/f'inflight{index:04d}.json'
            if outpath.exists():
                cell=load(outpath);assert cell['source_sha256']==reg['source_sha256'] and cell['interface_sha256']==state['interface_sha256']
            else:
                left_hi=Q(events[index]['isolating_interval'][1]);right_lo=Q(events[index+1]['isolating_interval'][0]);assert left_hi<right_lo
                denominator=1
                while True:
                    q=Q((left_hi*denominator).numerator//(left_hi*denominator).denominator+1,denominator)
                    if q<right_lo:break
                    denominator*=2
                # Disjoint outward rational value intervals certify a complete
                # strict ordering without repeatedly evaluating field differences.
                lows=[lo+m*q for (lo,hi),(c,m) in zip(const_bounds,forms)];highs=[hi+m*q for (lo,hi),(c,m) in zip(const_bounds,forms)]
                order=sorted(range(len(forms)),key=lambda i:lows[i]);assert all(highs[i]<lows[j] for i,j in zip(order,order[1:]))
                ranks=[0]*len(forms)
                for rank,fid in enumerate(order):ranks[fid]=rank
                orderhash=hashlib.sha256(json.dumps(order,separators=(',',':')).encode()).hexdigest()
                Avol=[ZERO]*25;Bvol=[ZERO]*25;dprefix=0;response=[]
                if inflight.exists():
                    ck=load(inflight);assert ck['source_sha256']==reg['source_sha256'] and ck['order_sha256']==orderhash
                    dprefix=ck['completed_descriptors'];Avol=decode(ck['partial_A']);Bvol=decode(ck['partial_B']);response=ck['descriptor_responses']
                def interval_responses(raw_pairs):
                    ints=[]
                    for lo,hi in raw_pairs:
                        lo=lo if ranks[lo]>ranks[zero_id] else zero_id;hi=hi if ranks[hi]<ranks[domain_id] else domain_id
                        if ranks[lo]<ranks[hi]:ints.append((lo,hi))
                    ints.sort(key=lambda p:(ranks[p[0]],ranks[p[1]]))
                    A=[ZERO]*25;B=0;l=h=None
                    def add_interval(left,right,target):
                        addto(target,[y-x for x,y in zip(forms[left][0],forms[right][0])])
                        return forms[right][1]-forms[left][1]
                    for lo,hi in ints:
                        if l is None:l,h=lo,hi
                        elif ranks[h]<ranks[lo]:B+=add_interval(l,h,A);l,h=lo,hi
                        elif ranks[h]<ranks[hi]:h=hi
                    if l is not None:B+=add_interval(l,h,A)
                    # Independent coordinate-cut coverage mask, entirely on the
                    # certified endpoint ordering; duplicates contribute once.
                    cuts=sorted({f for pair in ints for f in pair},key=lambda f:ranks[f]);AA=[ZERO]*25;BB=0
                    for lo,hi in zip(cuts,cuts[1:]):
                        if any(ranks[l]<=ranks[lo] and ranks[hi]<=ranks[h] for l,h in ints):BB+=add_interval(lo,hi,AA)
                    assert A==AA and B==BB
                    return A,B
                for di in range(dprefix,len(descriptors)):
                    dlA,dlB=interval_responses(descriptors[di]['raw_interval_form_pairs']);addto(Avol,mul(weights[di],dlA));addto(Bvol,[x*dlB for x in weights[di]])
                    response.append({'descriptor':di,'length_constant':encode(dlA),'length_slope':dlB});new_objects+=1;dprefix=di+1
                    if dprefix%16==0:save(inflight,{'source_sha256':reg['source_sha256'],'order_sha256':orderhash,'completed_descriptors':dprefix,'partial_A':encode(Avol),'partial_B':encode(Bvol),'descriptor_responses':response})
                A=[x*Q(16,3) for x in Avol];B=[x*Q(16,3) for x in Bvol]
                cell=certificate(index,A,B,{'status':'COMPLETE_exact_PWA_cell_dual_oneD_engines','ordering_only_interior_dyadic':str(q),'ordering_sha256':orderhash,'outward_value_intervals_strictly_ordered':True,'descriptor_dual_unions':91,'descriptor_responses':response,'interface_sha256':state['interface_sha256'],'batch_registration_sha256':sha(brpath),'no_parameter_point_union_or_interpolation':True})
                save(outpath,cell);actual_cells+=1
                if inflight.exists():inflight.unlink()
            residual=decode(cell['cell_max_residual']);event_id=cell['cell_max_event']
            if max_new is None or cmp(residual,max_new)>0:max_new=residual;max_new_event=event_id
            best=decode(state['best_global_residual'])
            if cmp(residual,best)>0:state['best_global_residual']=encode(residual);state['best_global_event']=event_id
            state['next_cell']=index+1;state['new_cells_computed']+=int(cell['status']!='INHERITED_exact_first_cell_no_union');state['status']='EXACT_PWA_prefix_saved'
            save(statepath,state);completed.append(index)
            lo,hi=bnd(residual)
            if lo>0:
                lock={'status':'STRICT_POSITIVE_nominal_finitefunctional_candidate_LOCK','cell':index,'event':event_id,'epsilon_coefficients':events[event_id]['value'],'collision_limit_only':event_id==len(events)-1,'cell_max_residual':cell['cell_max_residual'],'residual_interval':list(map(str,(lo,hi))),'cell_path':str(outpath),'strict_level_and_distinct_finite_epsilon_not_yet_certified':True,'source_sha256':reg['source_sha256']}
                save(P/'positive_candidate_LOCK.json',lock);print({'NEW_STRICT_POSITIVE_CANDIDATE_LOCK':lock},flush=True);break
        artifact=P/'batches'/br['cell_artifact']
        with artifact.open('wb') as raw:
            with gzip.GzipFile(fileobj=raw,mode='wb',mtime=0) as out:out.write((json.dumps({'cells':[load(scratch/f'cell{i:04d}.json') for i in completed]},separators=(',',':'))+'\n').encode())
        ratio_new=[c+d for c,d in zip(C,max_new)] if max_new is not None else C
        receipt={'status':'POSITIVE_LOCK_stop_for_certification' if (P/'positive_candidate_LOCK.json').exists() else 'COMPLETE_registered_batch','first_cell':first,'last_requested_cell':last,'completed_cell_indices':completed,'new_cells_computed_this_process':actual_cells,'new_descriptor_dual_union_objects':new_objects,'max_new_residual_coefficients':encode(max_new) if max_new is not None else {},'max_new_residual_interval':list(map(str,bnd(max_new))) if max_new is not None else ['0','0'],'max_new_ratio_interval':list(map(str,bnd(ratio_new))),'max_new_ratio_display_only':[float(x) for x in bnd(ratio_new)],'max_new_event':max_new_event,'covered_prefix_cells':state['next_cell'],'remaining_cells':N-state['next_cell'],'global_best_residual':state['best_global_residual'],'cell_artifact':artifact.name,'cell_artifact_sha256':sha(artifact),'batch_registration_sha256':sha(brpath),'source_sha256':reg['source_sha256'],'runtime':rt(start)}
        save(P/'batches'/br['receipt'],receipt);print({k:v for k,v in receipt.items() if k not in ('max_new_residual_coefficients','max_new_residual_interval','max_new_ratio_interval')},flush=True)
    except Exception as ex:
        save(P/f'FAULT_cell{fault_cell:04d}.json',{'status':'ACTUAL_FAULT_no_automatic_restart','cell':fault_cell,'type':type(ex).__name__,'message':str(ex),'traceback':traceback.format_exc(),'exact_state_sha256':sha(statepath),'source_sha256':reg['source_sha256'],'runtime':rt(start)});raise

if __name__=='__main__':
    assert all(os.environ.get(k)=='1' for k in ENV) and os.getpriority(os.PRIO_PROCESS,0)>=10
    if len(sys.argv)<2:raise SystemExit('usage: exact_cells.py init|batch registration_filename')
    main(sys.argv[1],sys.argv[2] if len(sys.argv)>2 else None)
