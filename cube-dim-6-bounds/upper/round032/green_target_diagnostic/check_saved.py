#!/usr/bin/env python3
"""Exact saved B032 first-stop certificate. Never runs a kernel."""
import sys
sys.dont_write_bytecode=True
sys.set_int_max_str_digits(0)
import hashlib,importlib.util,itertools,json,os,time
from fractions import Fraction as Q
from pathlib import Path
HERE=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('B015_preflight',HERE.parents[2]/'upper/round015/spatial_minimum_search/prepare.py')
prep=importlib.util.module_from_spec(spec);spec.loader.exec_module(prep)

def read(path):return json.loads(Path(path).read_text())
def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def interval(row):
    assert set(row)=={'lower_mantissa','upper_mantissa','common_binary_exponent'}
    a,b,e=map(int,(row['lower_mantissa'],row['upper_mantissa'],row['common_binary_exponent']))
    assert a<=b and abs(e)<10000
    return Q(a)*Q(2)**e,Q(b)*Q(2)**e
def rounded(q):
    n,r=divmod(q.numerator,q.denominator)
    if 2*r>q.denominator or (2*r==q.denominator and n%2):n+=1
    return n
def key(nums,den):return tuple(str(Q(n,den))for n in nums)

def check():
    tic=time.perf_counter();reg=read(HERE/'registration.json');run=read(HERE/'run_receipt.json');res=read(HERE/'results.json')
    for path,digest in reg['source_sha256'].items():assert sha(path)==digest
    assert sha(HERE/'run_batch.py')==reg['driver_sha256']
    if(HERE/'manifest.json').exists():
        for name,digest in read(HERE/'manifest.json')['whitelist_sha256'].items():assert sha(HERE/name)==digest
    assert run['registration_sha256']==sha(HERE/'registration.json') and run['results_sha256']==sha(HERE/'results.json')
    assert run['status']=='COMPLETED' and run['exit_code']==0 and run['all_child_handles_closed']
    assert run['elapsed_seconds']<1200 and not run['starting_point_recomputed']
    assert len(run['native_calls'])==len(res['points'])==1
    assert res['status']=='STRICT_G0_COUNTEREXAMPLE_IMMEDIATE_STOP' and res['counterexample_index']==1
    assert not(HERE/'candidate_02.json').exists() and not(HERE/'native_outputs/02.stdout.json').exists()
    assert not res['starting_point_recomputed'] and not res['new_C6_bound']
    assert reg['bounds_snapshot']['public_upper']==reg['bounds_snapshot']['internal_upper']=='64'
    record=run['native_calls'][0];row=res['points'][0]
    assert record['exit_code']==0 and record['elapsed_seconds']<120 and record['status']=='COMPLETE_VALIDATED'
    assert sha(HERE/record['stdout'])==record['stdout_sha256'] and sha(HERE/record['stderr'])==record['stderr_sha256']
    assert (HERE/record['stderr']).read_bytes()==b''
    assert sha(HERE/record['candidate_path'])==record['candidate_sha256']
    nums,den=record['numerators'],record['denominator']
    assert nums==reg['first_candidate_numerators']==row['numerators'] and den==row['denominator']==10**14
    assert record['command']==[reg['binary_path'],'full',reg['fan_path'],str(den),*map(str,nums)]
    oldkeys=set(tuple(p)for p in reg['old_G_exact_point_keys']);assert len(oldkeys)==158
    assert key(nums,den)not in oldkeys and key(reg['starting_numerators'],reg['starting_denominator'])in oldkeys
    _,rays=prep.fan_data();forms=prep.all_zero_forms(rays)
    assert len(forms)==9308 and hashlib.sha256(json.dumps(forms,separators=(',',':')).encode()).hexdigest()==reg['zero_forms_sha256']
    pf=prep.preflight(nums,den,forms,rays)
    assert pf==reg['first_final_preflight'] and pf['feasible'] and pf['int64_safe'] and pf['zero_form_count']==0
    old=read(reg['starting_native_path']);new=read(HERE/record['stdout'])
    assert new['status']=='generic_nonzero_supported' and new['precision_bits']==512
    assert new['full_group_terms']==967680 and new['computed_sign_pairs']==483840
    assert new['zero_node_terms']==new['confluent_terms']==0 and len(new['gradient'])==6
    g=interval(new['G']);grad=[interval(v)for v in new['gradient']];euler=interval(new['Euler_x_grad_plus5G'])
    assert 0<g[0]<=g[1]<Q(1,60000) and euler[0]<=0<=euler[1]
    for a,b in [g,*grad]:assert b-a<=max(Q(1),abs((a+b)/2))*Q(2)**-256
    lo,hi=5*g[0],5*g[1]
    for n,(a,b)in zip(nums,grad):lo+=Q(n,den)*a;hi+=Q(n,den)*b
    assert lo<=0<=hi
    saved=row['check'];assert saved['rejects_g0_strictly'] and not saved['strictly_above_g0']
    assert list(map(Q,saved['G_interval']))==list(g)
    assert [list(map(Q,p))for p in saved['gradient_intervals']]==[list(p)for p in grad]
    assert list(map(Q,saved['independent_Euler_interval']))==[lo,hi]
    assert Q(saved['G_upper_minus_g0'])==g[1]-Q(1,60000)<0
    # Independently reconstruct the one fixed descent, clipping and exact repair.
    c=read(HERE/'candidate_01.json');start=reg['starting_numerators'];startden=reg['starting_denominator']
    oldgrad=[interval(v)for v in old['gradient']];mid=[(a+b)/2 for a,b in oldgrad];maximum=max(map(abs,mid[1:]))
    assert maximum>0 and list(map(Q,c['gradient_exact_dyadic_midpoints']))==mid
    direction=[Q(0)]+[-a/maximum for a in mid[1:]]
    raw=[Q(n,startden)+d/20 for n,d in zip(start,direction)]
    assert list(map(Q,c['negative_normalized_direction']))==direction and list(map(Q,c['raw_rational_proposal']))==raw
    assert c['fixed_step']=='1/20' and c['starting_native_sha256']==sha(reg['starting_native_path'])
    gap=Q(1,10**8);clip=[Q(1)]
    for i in range(1,6):
        value=raw[i]
        if value<(6-i)*gap:value=(6-i)*gap
        if value>1-i*gap:value=1-i*gap
        if value>clip[-1]-gap:value=clip[-1]-gap
        clip.append(value)
    assert list(map(Q,c['clipped_rational_proposal']))==clip
    quantized=[rounded(v*den)for v in clip];assert quantized==c['quantized_numerators']
    assert c['initial_preflight']==prep.preflight(quantized,den,forms,rays)==reg['first_initial_preflight']
    assert c['initial_preflight']['zero_forms']==[[1,-2,1,0,0,0]]
    assert len(c['repairs'])==1 and c['repairs'][0]['k']==1
    gaps=[quantized[i]-quantized[i+1]for i in range(5)]+[quantized[-1]]
    j=max(range(6),key=lambda i:gaps[i]);deltas=[11**i if i!=j else 0 for i in range(6)];deltas[j]=-sum(deltas)
    assert j==4 and deltas==c['repairs'][0]['gap_delta_numerators'] and sum(deltas)==0
    newgaps=[a+b for a,b in zip(gaps,deltas)]
    rebuilt=[den]+[den-sum(newgaps[:i])for i in range(1,6)]
    assert rebuilt==nums==c['final_numerators'] and c['final_preflight']==pf
    # True S2 distance and comparison with the entire frozen old Green chain.
    x=[Q(n,den)for n in nums]
    distance=max(min(abs(x[a]-x[b])for a,b in itertools.combinations(triple,2))/2 for triple in itertools.combinations(range(6),3))
    assert distance==Q(saved['exact_dist_infinity_to_S2_lower_bound'])>0
    oldg=interval(old['G']);assert g[1]<oldg[0] and row['strictly_below_previous']
    assert list(map(Q,row['G_minus_previous_interval']))==[g[0]-oldg[1],g[1]-oldg[0]]
    old_res=read(Path(reg['starting_native_path']).parents[1]/'results.json')
    oldz1=list(map(Q,old_res['points'][0]['check']['G_interval']))
    oldx=list(map(Q,old_res['old_G_intervals']['x0']));oldc=list(map(Q,old_res['old_G_intervals']['old_c']))
    assert oldg[1]<oldz1[0] and oldz1[1]<oldx[0] and oldx[1]<oldc[0]
    defect=Q(1,60000)-g[1];simple=Q((defect*10**20).__floor__(),10**20)
    assert 0<simple<defect
    rational_lo=Q((g[0]*10**20).__floor__(),10**20);rational_hi=Q((g[1]*10**20).__ceil__(),10**20)
    assert 0<rational_lo<g[0]<=g[1]<rational_hi<Q(1,60000)
    return dict(status='PASS_STRICT_G0_TARGET_COUNTEREXAMPLE_FIRST_STOP',pid=os.getpid(),elapsed_seconds=time.perf_counter()-tic,
                original_native_calls=1,new_checker_native_calls=0,starting_point_recomputed=False,
                G_rational_display_bounds=[str(rational_lo),str(rational_hi)],
                strict_g0_minus_G_lower=str(simple),exact_g0_minus_Gupper=str(defect),
                exact_distance_to_true_S2_lower=str(distance),complete_logical_terms=967680,
                full_zero_form_count=9308,confluent_terms=0,first_stop_observed=True,
                full_old_G_chain_verified=True,refutes_global_g_at_least_1_over_60000=True,
                excludes_true_I_over_g_model=False,public_C6_upper_unchanged='64',internal_C6_upper_unchanged='64',new_C6_bound=False)

if __name__=='__main__':
    assert sys.argv[1:]in([],['--write-receipt'])
    result=check()
    if sys.argv[1:]==['--write-receipt']:
        assert not(HERE/'checked_results.json').exists()
        (HERE/'checked_results.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items()if k!='exact_g0_minus_Gupper'}),flush=True)
