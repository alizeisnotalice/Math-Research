#!/usr/bin/env python3
"""B030 pure saved certificate checker; no engine call or new Green value."""
import sys
sys.dont_write_bytecode=True
sys.set_int_max_str_digits(0)
import hashlib,importlib.util,itertools,json,os,time
from fractions import Fraction as Q
from pathlib import Path

HERE=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('frozen_preflight',HERE.parents[2]/'upper/round015/spatial_minimum_search/prepare.py')
prep=importlib.util.module_from_spec(spec);spec.loader.exec_module(prep)

def read(path):return json.loads(Path(path).read_text())
def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def interval(row):
    assert set(row)=={'lower_mantissa','upper_mantissa','common_binary_exponent'}
    a,b,e=map(int,(row['lower_mantissa'],row['upper_mantissa'],row['common_binary_exponent']))
    assert a<=b and abs(e)<10000
    return Q(a)*Q(2)**e,Q(b)*Q(2)**e
def round_even(q):
    n,r=divmod(q.numerator,q.denominator)
    return n+(2*r>q.denominator or (2*r==q.denominator and n%2==1))
def canonical(nums,den):return tuple(str(Q(n,den))for n in nums)

def check():
    tic=time.perf_counter();reg=read(HERE/'registration.json');run=read(HERE/'run_receipt.json');res=read(HERE/'results.json')
    for path,digest in reg['source_sha256'].items():assert sha(path)==digest
    assert sha(HERE/'diagnose.py')==reg['driver_sha256']
    if(HERE/'manifest.json').exists():
        for path,digest in read(HERE/'manifest.json')['whitelist_sha256'].items():assert sha(HERE/path)==digest
    assert run['registration_sha256']==sha(HERE/'registration.json') and run['results_sha256']==sha(HERE/'results.json')
    assert run['status']=='COMPLETED' and run['exit_code']==0 and run['all_child_handles_closed']
    assert run['elapsed_seconds']<1200 and len(run['native_calls'])==len(res['points'])==2
    assert res['status']=='TWO_POINT_DIAGNOSTIC_COMPLETE_GLOBAL_G_UNRESOLVED'
    assert reg['bounds_snapshot']['upper_exact']==reg['bounds_snapshot']['public_upper_exact']==reg['bounds_snapshot']['internal_upper_exact']=='64'
    assert reg['g0']=='1/60000' and reg['max_native']==2
    _,rays=prep.fan_data();forms=prep.all_zero_forms(rays)
    assert len(forms)==reg['zero_form_count']==9308
    assert hashlib.sha256(json.dumps(forms,separators=(',',':')).encode()).hexdigest()==reg['zero_forms_sha256']
    oldkeys=set(tuple(row)for row in reg['old_G_exact_point_keys']);values=[];gradients=[]
    summaries=[]
    for i,(record,row)in enumerate(zip(run['native_calls'],res['points']),1):
        nums,den=record['numerators'],record['denominator']
        assert nums==row['numerators'] and den==row['denominator']==10**14
        key=canonical(nums,den);assert key not in oldkeys;oldkeys.add(key)
        assert record['exact_key']==list(key) and record['exit_code']==0
        assert record['elapsed_seconds']<120 and record['status']=='COMPLETE_VALIDATED'
        assert record['command']==[reg['binary_path'],'full',reg['fan_path'],str(den),*map(str,nums)]
        assert sha(HERE/record['stdout'])==record['stdout_sha256'] and sha(HERE/record['stderr'])==record['stderr_sha256']
        assert (HERE/record['stderr']).read_bytes()==b''
        pf=prep.preflight(nums,den,forms,rays)
        assert pf==record['preflight'] and pf['feasible'] and pf['int64_safe'] and pf['zero_form_count']==0
        out=read(HERE/record['stdout'])
        assert out['status']=='generic_nonzero_supported' and out['precision_bits']==512
        assert out['full_group_terms']==967680 and out['computed_sign_pairs']==483840 and out['zero_node_terms']==0
        assert out['confluent_terms']==0 and len(out['gradient'])==6
        g=interval(out['G']);grad=[interval(r)for r in out['gradient']];ev=interval(out['Euler_x_grad_plus5G'])
        assert g[0]>Q(1,60000) and ev[0]<=0<=ev[1]
        el,eu=5*g[0],5*g[1]
        for x,(a,b)in zip(nums,grad):
            el+=Q(x,den)*a;eu+=Q(x,den)*b
        assert el<=0<=eu
        for a,b in [g,*grad]:assert b-a<=max(Q(1),abs((a+b)/2))*Q(2)**-256
        saved=row['check']
        assert [str(v)for v in g]==saved['G_interval']
        assert [[str(v)for v in pair]for pair in grad]==saved['gradient_intervals']
        assert [str(el),str(eu)]==saved['independent_Euler_interval']
        assert not saved['rejects_g0_strictly'] and saved['strictly_above_g0']
        assert Q(saved['G_lower_minus_g0'])==g[0]-Q(1,60000)>0
        assert Q(saved['G_upper_minus_g0'])==g[1]-Q(1,60000)>0
        x=[Q(n,den)for n in nums]
        distance=max(min(abs(x[a]-x[b])for a,b in itertools.combinations(triple,2))/2 for triple in itertools.combinations(range(6),3))
        assert distance==Q(saved['exact_dist_infinity_to_S2_lower_bound'])>0
        assert Q(saved['necessary_B_ceiling_for_cost_below64'])==Q(16,3)*g[1]
        scale=10**20
        display_bounds=[Q((g[0]*scale).__floor__(),scale),Q((g[1]*scale).__ceil__(),scale)]
        assert display_bounds[0]<g[0]<=g[1]<display_bounds[1]
        summary={'index':i,'strict_rational_G_display_bounds':list(map(str,display_bounds)),
                 'strict_margin_above_g0_lower':str(Q(((g[0]-Q(1,60000))*10**18).__floor__(),10**18)),
                 'distance_to_S2_lower':str(distance),'confluent_terms':out['confluent_terms']}
        summaries.append(summary);values.append(g);gradients.append(grad)
    # Reconstruct the sole candidate from the first native GREEN gradient.
    candidate=read(HERE/'candidate.json');seed=reg['seed_numerators'];den=reg['seed_denominator']
    assert run['native_calls'][0]['numerators']==seed
    mid=[(a+b)/2 for a,b in gradients[0]];maximum=max(map(abs,mid[1:]))
    assert maximum>0 and list(map(Q,candidate['gradient_exact_dyadic_midpoints']))==mid
    direction=[Q(0)]+[-v/maximum for v in mid[1:]]
    raw=[Q(n,den)+v/20 for n,v in zip(seed,direction)]
    assert list(map(Q,candidate['negative_normalized_direction']))==direction
    assert list(map(Q,candidate['raw_rational_proposal']))==raw and candidate['fixed_step']=='1/20'
    gap=Q(1,10**8);clipped=[Q(1)]
    for i in range(1,6):
        value=raw[i]
        if value<(6-i)*gap:value=(6-i)*gap
        if value>1-i*gap:value=1-i*gap
        if value>clipped[-1]-gap:value=clipped[-1]-gap
        clipped.append(value)
    assert list(map(Q,candidate['clipped_rational_proposal']))==clipped
    quantized=[round_even(q*10**14)for q in clipped]
    assert quantized==candidate['quantized_numerators']
    assert candidate['initial_preflight']==prep.preflight(quantized,10**14,forms,rays)
    assert candidate['initial_preflight']['zero_forms']==[[1,-2,1,0,0,0]]
    # Only k=1 was needed; neither this rejected input nor k>=2 reached native.
    assert len(candidate['repairs'])==1
    repair=candidate['repairs'][0];gaps=[quantized[i]-quantized[i+1]for i in range(5)]+[quantized[-1]]
    j=max(range(6),key=lambda i:gaps[i]);delta=[11**i if i!=j else 0 for i in range(6)];delta[j]=-sum(delta)
    assert repair['k']==1 and repair['largest_gap_index']==j==4
    assert delta==repair['gap_delta_numerators'] and sum(delta)==0
    newg=[g+d for g,d in zip(gaps,delta)]
    point=[10**14]+[10**14-sum(newg[:i])for i in range(1,6)]
    assert point==repair['numerators']==candidate['final_numerators']==run['native_calls'][1]['numerators']
    assert candidate['final_preflight']==prep.preflight(point,10**14,forms,rays)
    assert candidate['final_preflight']['zero_form_count']==0
    displacement=[Q(a,10**14)-Q(b,den)for a,b in zip(point,seed)]
    linear=[Q(0),Q(0)]
    for dx,(a,b)in zip(displacement,gradients[0]):
        products=sorted([dx*a,dx*b]);linear[0]+=products[0];linear[1]+=products[1]
    assert linear[1]<0
    assert Q(candidate['clipped_step_gradient_dot_displacement'])==sum(a*b for a,b in zip(mid,displacement))
    oldc=interval(read(reg['old_c_native_path'])['G'])
    const=read(reg['old_x0_constants_path'])['G_from_old_Gamma']
    oldx=(Q(int(const['lower_mantissa']))*Q(2)**int(const['lower_binary_exponent']),Q(int(const['upper_mantissa']))*Q(2)**int(const['upper_binary_exponent']))
    assert res['old_G_intervals']=={'old_c':list(map(str,oldc)),'x0':list(map(str,oldx))}
    assert oldx[1]<oldc[0]
    for g,row in zip(values,res['points']):
        assert g[1]<oldx[0]<oldc[0]
        assert row['strictly_below_old_c'] and row['strictly_below_x0'] and not row['strictly_above_x0']
    assert values[1][1]<values[0][0] and res['second_strictly_below_first']
    delta_g=[values[1][0]-values[0][1],values[1][1]-values[0][0]]
    assert list(map(Q,res['second_minus_first_interval']))==delta_g
    assert not res['new_C6_bound'] and not run['new_C6_bound']
    return dict(status='PASS_BOUNDED_GREEN_DIAGNOSTIC_G0_UNRESOLVED',pid=os.getpid(),elapsed_seconds=time.perf_counter()-tic,
                native_calls_in_original_batch=2,new_native_calls_in_checker=0,new_G_evaluations_in_checker=0,
                complete_terms_each=967680,zero_forms_checked_each=9308,points=summaries,
                actual_displacement_first_gradient_dot_interval=list(map(str,linear)),
                G2_strictly_less_G1_strictly_less_Gx0=True,g0_strictly_less_G2=True,
                public_upper_unchanged='64',internal_upper_unchanged='64',new_C6_bound=False)

if __name__=='__main__':
    assert sys.argv[1:]in([],['--write-receipt'])
    result=check()
    if sys.argv[1:]==['--write-receipt']:
        assert not(HERE/'checked_results.json').exists()
        (HERE/'checked_results.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items()if k!='actual_displacement_first_gradient_dot_interval'}),flush=True)
