#!/usr/bin/env python3
"""Exact a2 event inventory, first actual raw event, one static certificate."""
from fractions import Fraction as Q
from functools import lru_cache
from itertools import product
from math import prod
from pathlib import Path
import hashlib
import importlib.util
import json
import os
import sys
import time

sys.dont_write_bytecode=True
HERE=Path(__file__).resolve().parent
LOWER=HERE.parents[1]
ARITHMETIC=LOWER/'round001_rectangular_certificate.py'
CURRENT=Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube6-lower-round005-endpointfix/cube-dim-6-bounds/audit/lower007/candidate_number_field/candidate_results.json')
JET=LOWER/'round007/shape_jets/jet_results.json'
spec=importlib.util.spec_from_file_location('chamber_field',ARITHMETIC)
A=importlib.util.module_from_spec(spec)
spec.loader.exec_module(A)
A.D=10**80
A.RU,A.RV=A.floor_root(2*A.D**6,6),A.floor_root(3*A.D**6,6)
A.UL,A.UH=Q(A.RU,A.D),Q(A.RU+1,A.D)
A.VL,A.VH=Q(A.RV,A.D),Q(A.RV+1,A.D)
assert A.UL**6<=2<A.UH**6 and A.VL**6<=3<A.VH**6
A.bounds.cache_clear()
F,u,v=A.F,A.u,A.v
ZERO,ONE=F(0),F(1)
STARTED=time.monotonic()
NODES=0
NODE_CAP,WALL_CAP,EVENT_CAP=1000000,1200,12

def budget():
    if time.monotonic()-STARTED>WALL_CAP:
        raise TimeoutError('global 1200-second cap')
    if NODES>NODE_CAP:
        raise RuntimeError('global million union-state cap')

def coeff(x,full=False):
    return {f'{i},{j}':str(x.a[6*i+j]) for i in range(6) for j in range(6)
            if full or x.a[6*i+j]}

def field(values):
    result=[Q(0)]*36
    for ij,x in values.items():
        i,j=map(int,ij.split(',')); result[6*i+j]=Q(x)
    return F(a=result)

def interval(x):
    return [str(z) for z in A.bounds(x.a)]

@lru_cache(None)
def root_radius(m):
    r=A.root(m)/2
    assert (2*r)**6==F(m)
    return r

def build_family(spacing):
    raw=[]; no_other=0; no_moving=0
    for ell in product(range(1,5),repeat=6):
        budget()
        r=root_radius(prod(ell))
        starts=tuple(-(l//2) if l%2 else 1-l//2 for l in ell)
        box=[None]*12
        for j in (0,2,3,4,5):
            aj,q,l=spacing[j],starts[j],ell[j]
            lo,hi=max(ZERO,(q+l-1)*aj-r),min(aj/2,q*aj+r)
            if hi<=lo:
                no_other+=1; break
            box[2*j:2*j+2]=(lo,hi)
        else:
            s,l=spacing[1],ell[1]
            lo,hi=((ZERO,min(s/2,r-(l-1)*s/2)) if l%2 else
                   (max(ZERO,l*s/2-r),s/2))
            if hi<=lo:
                assert l>=2 and 2*r<=(l-1)*s
                no_moving+=1; continue
            assert l==1 or (l-1)*s<2*r
            box[2:4]=(lo,hi)
            raw.append({'ell':ell,'starts':starts,'radius':r,'base_box':tuple(box)})
    assert len(raw)+no_other+no_moving==4096
    return raw,{'length_vectors':4096,'raw_nonempty':len(raw),
                'empty_in_fixed_coordinates':no_other,
                'empty_in_moving_coordinate_never_born_for_larger_s':no_moving}

def make_functions(raw):
    functions=[]; index={}
    def add(c,m,label):
        key=(c,Q(m))
        if key not in index:
            index[key]=len(functions)
            functions.append({'c':c,'m':Q(m),'labels':[]})
        j=index[key]; functions[j]['labels'].append(label)
        return j
    zero=add(ZERO,0,{'face':'zero'})
    cap=add(ZERO,Q(1,2),{'face':'half_cell'})
    for k,record in enumerate(raw):
        l,r=record['ell'][1],record['radius']
        label={'raw_index':k,'lengths':record['ell'],'mass':prod(record['ell'])}
        if l%2:
            moving=add(r,-Q(l-1,2),{**label,'face':'odd_upper'})
            record['low_options']=(zero,); record['high_options']=(cap,moving)
        else:
            moving=add(-r,Q(l,2),{**label,'face':'even_lower'})
            record['low_options']=(zero,moving); record['high_options']=(cap,)
    return functions

def state_at_event(raw,functions,event,side):
    keys=[(f['c']+f['m']*event,side*f['m']) for f in functions]
    selections={}; active=set()
    for k,r in enumerate(raw):
        low=max(r['low_options'],key=lambda i:keys[i])
        high=min(r['high_options'],key=lambda i:keys[i])
        if keys[low]<keys[high]:
            selections[k]=(low,high); active.update((low,high))
    return selections,tuple(sorted(active,key=lambda i:keys[i]))

def event_inventory(functions,current):
    events={}
    for i,f in enumerate(functions):
        budget()
        for j,g in enumerate(functions[:i]):
            if f['m']==g['m']: continue
            crossing=(f['c']-g['c'])/(g['m']-f['m'])
            if current<crossing:
                events.setdefault(crossing,[]).append((j,i))
    return sorted(events),events

def make_interface(raw,functions,current,event,old_upper):
    delta=event-current
    chosen,_=state_at_event(raw,functions,current+delta/2,1)
    assert len(chosen)==len(raw)
    endpoints={}; endpoint_ids={}; by_axis=[set() for _ in range(6)]; boxes=[]
    def add_endpoint(c0,c1):
        key=(c0,c1)
        if key not in endpoint_ids:
            name=f'e{len(endpoint_ids):04d}'
            endpoint_ids[key]=name
            endpoints[name]={'c0':coeff(c0,True),'c1':coeff(c1,True)}
        return endpoint_ids[key]
    for k,r in enumerate(raw):
        coordinates=[]
        for axis in range(6):
            pair=[]
            for side in (0,1):
                if axis==1:
                    f=functions[chosen[k][side]]
                    c0,c1=f['c']+f['m']*current,f['m']*delta
                else:
                    c0,c1=r['base_box'][2*axis+side],ZERO
                name=add_endpoint(c0,c1)
                pair.append(name); by_axis[axis].add(name)
            coordinates.append(pair)
        boxes.append(coordinates)
    reverse={name:key for key,name in endpoint_ids.items()}
    orders=[]; checks=0; boundary_ties=0
    for names in by_axis:
        ordered=sorted(names,key=lambda name:reverse[name][0]+reverse[name][1]/2)
        for left,right in zip(ordered,ordered[1:]):
            lc,ld=reverse[left]; rc,rd=reverse[right]
            g0,g1=rc-lc,rc+rd-lc-ld
            assert ZERO<=g0 and ZERO<=g1
            assert ZERO<g0 or ZERO<g1
            checks+=1; boundary_ties+=g0==ZERO or g1==ZERO
        orders.append([[name] for name in ordered])
    return {'relations':['u^6=2, u>0','v^6=3, v>0'],'interval':['0','1'],
            'dimension':6,'normalization_factor':'64','candidate_parameter':'1',
            'static_check_parameter':'1/100',
            'source':'unit rectangular lattice; all 127 current raw boxes',
            'spacing2':{'c0':coeff(current,True),'c1':coeff(delta,True)},
            'current_Rideal_complete_upper':str(old_upper),
            'endpoints':endpoints,'boxes':boxes,'box_lengths':[r['ell'] for r in raw],
            'order_certificate':orders,'order_checks':checks,
            'adjacent_boundary_ties':boundary_ties,
            'tie_convention':'globally deduplicated identical affine functions; singleton groups suffice'}

def candidate_union(raw,spacing):
    boxes=[]; endpoint_checks=0
    for r in raw:
        box=[]
        for aj,q,l in zip(spacing,r['starts'],r['ell']):
            lo,hi=max(ZERO,(q+l-1)*aj-r['radius']),min(aj/2,q*aj+r['radius'])
            if hi<=lo: break
            for x in (lo,hi):
                for index in (q,q+l-1):
                    assert -r['radius']<=x-index*aj<=r['radius']
                    endpoint_checks+=1
            box.extend((lo,hi))
        else: boxes.append(tuple(box))
    @lru_cache(None)
    def union(B):
        global NODES
        NODES+=1; budget()
        if not B: return ZERO
        if len(B[0])==2:
            ordered=sorted(B); lo,hi=ordered[0]; value=ZERO
            for c,d in ordered[1:]:
                if hi<c:
                    value=value+hi-lo; lo,hi=c,d
                else: hi=max(hi,d)
            return value+hi-lo
        cuts=sorted({x for b in B for x in b[:2]}); value=ZERO
        for lo,hi in zip(cuts,cuts[1:]):
            child=tuple(sorted({b[2:] for b in B if b[0]<=lo and hi<=b[1]}))
            value=value+(hi-lo)*union(child)
        return value
    value=64*union(tuple(sorted(set(boxes))))
    return value,{'positive_raw_at_event':len(boxes),'unique_at_event':len(set(boxes)),
                  'source_endpoint_checks':endpoint_checks}

def main():
    inputs=[CURRENT,ARITHMETIC,JET,HERE/'registration.md',Path(__file__)]
    hashes={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs}
    old=json.loads(CURRENT.read_text()); old_value=field(old['coefficients'])
    old_upper=Q(old['exact_interval'][1])
    derivative=field(json.loads(JET.read_text())['positive_direction']['derivative_coefficients'])
    a0=(u**3+u**4*v**3)/4; current=a0*Q(1000001,1000000); w=u*v**3
    inverse_a0=-2*u**3*sum(((-w)**k for k in range(6)),ZERO)/53
    assert a0*inverse_a0==ONE
    slope=derivative*inverse_a0
    assert ZERO<slope
    spacing=(ONE+u**4/2,current,*((ONE+u)/2 for _ in range(4)))
    out={'status':'running','pid':os.getpid(),'input_sha256':hashes,
         'current_Rideal_exact_interval':old['exact_interval'],
         'current_Rideal_coefficients':old['coefficients'],
         'current_spacing2_coefficients':coeff(current),
         'spacing_derivative_coefficients':coeff(slope),
         'limits':{'event_candidates':EVENT_CAP,'full_candidate_unions':1,
                   'union_states':NODE_CAP,'wall_seconds':WALL_CAP},
         'events_examined':[],'full_candidate_unions':0}
    def save():
        out['elapsed_seconds']=time.monotonic()-STARTED; out['union_states']=NODES
        (HERE/'chamber_results.json').write_text(json.dumps(out,indent=2)+'\n')
    save()
    try:
        raw,family=build_family(spacing)
        assert len(raw)==old['nonempty_raw_boxes']==127
        out['family']=family; functions=make_functions(raw)
        events,partners=event_inventory(functions,current)
        out['distinct_unclipped_functions']=len(functions)
        out['potential_future_event_locations']=len(events); save()
        selected=None
        for event in events[:EVENT_CAP]:
            left,lorder=state_at_event(raw,functions,event,-1)
            right,rorder=state_at_event(raw,functions,event,1)
            kinds=[]
            if set(left)!=set(right): kinds.append('raw_box_death')
            if any(left[k]!=right[k] for k in left.keys()&right.keys()):
                kinds.append('clipping_branch_switch')
            if lorder!=rorder: kinds.append('effective_endpoint_order_change')
            out['events_examined'].append({'spacing2_coefficients':coeff(event),
                'spacing2_interval':interval(event),
                'spacing2_float':float(sum(A.bounds(event.a),Q(0))/2),
                'actual_raw_event':bool(kinds),'event_types':kinds,
                'intersecting_function_pairs':partners[event],
                'raw_left':len(left),'raw_right':len(right)})
            save()
            if kinds: selected=event; break
        if selected is None:
            raise RuntimeError('no actual raw event within 12 examined locations')
        event=selected
        interface=make_interface(raw,functions,current,event,old_upper)
        (HERE/'chamber_interface.json').write_text(json.dumps(interface,indent=2)+'\n')
        out['interface_sha256']=hashlib.sha256((HERE/'chamber_interface.json').read_bytes()).hexdigest()
        out['chamber']={'parameter':'s=s_current+(s_event-s_current)*t','interval':['0','1'],
            'current_spacing2_interval':interval(current),'event_spacing2_interval':interval(event),
            'event_spacing2_coefficients':coeff(event),'order_checks':interface['order_checks'],
            'adjacent_boundary_ties':interface['adjacent_boundary_ties'],
            'all_current_raw_boxes_remain_positive_interior':True,
            'births_for_larger_s':'none by exact interval-width monotonicity',
            'event_is_first_union_volume_kink':'not asserted; first actual raw-arrangement event'}
        out['selected_event_functions']=[{
            'left':{'c':coeff(functions[i]['c']),'slope':str(functions[i]['m']),
                    'labels':functions[i]['labels']},
            'right':{'c':coeff(functions[j]['c']),'slope':str(functions[j]['m']),
                     'labels':functions[j]['labels']}} for i,j in partners[event]]
        predicted=old_value+slope*(event-current)
        out['full_candidate_unions']=1
        actual,counts=candidate_union(raw,(spacing[0],event,*spacing[2:]))
        assert actual==predicted,'static candidate differs from chamber affine identity'
        lo,hi=A.bounds(actual.a)
        out['candidate']={'spacing2_coefficients':coeff(event),
            'other_spacings':['1+u^4/2','(1+u)/2 x4'],
            'parameter':'t=1 in chamber_interface.json','exact_interval':[str(lo),str(hi)],
            'coefficients':coeff(actual),'ratio_float':float((lo+hi)/2),
            'strictly_above_current_complete_Rideal_upper':lo>old_upper,
            'margin_exact_lower':str(lo-old_upper),'margin_float':float(lo-old_upper),
            'full_static_equals_affine_prediction':True,**counts}
        out['status']=('promoted_sole_candidate_stop_pending_E_L1_and_independent_audit'
                       if lo>old_upper else 'sole_candidate_did_not_promote_stop')
    except (RuntimeError,TimeoutError,ArithmeticError,AssertionError) as exc:
        out['status']='stopped'; out['reason']=repr(exc)
    out['actual_process_status']='finished; no MATLAB, child process, or continuing event search'
    save()
    print(json.dumps({'status':out['status'],'events_examined':len(out['events_examined']),
        'full_candidate_unions':out['full_candidate_unions'],'union_states':NODES,
        'elapsed_seconds':out['elapsed_seconds'],
        'candidate_ratio_float':out.get('candidate',{}).get('ratio_float'),
        'candidate_margin_float':out.get('candidate',{}).get('margin_float'),
        'reason':out.get('reason')}),flush=True)

if __name__=='__main__': main()
