#!/usr/bin/env python3
"""Pure saved-data checker for A038; no native or other kernel invocation."""
import sys
sys.dont_write_bytecode=True
sys.set_int_max_str_digits(0)
import hashlib,importlib.util,itertools,json,time,os
from fractions import Fraction as Q
from pathlib import Path
HERE=Path(__file__).resolve().parent
PREP=Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube6-upper-astra/cube-dim-6-bounds/upper/round015/spatial_minimum_search/prepare.py')
spec=importlib.util.spec_from_file_location('frozen_zero_preflight',PREP)
prep=importlib.util.module_from_spec(spec);spec.loader.exec_module(prep)
def read(p):return json.loads(Path(p).read_text())
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def interval(a):
    assert set(a)=={'lower_mantissa','upper_mantissa','common_binary_exponent'}
    lo,hi,e=map(int,(a['lower_mantissa'],a['upper_mantissa'],a['common_binary_exponent']))
    assert lo<=hi and abs(e)<10000
    return Q(lo)*Q(2)**e,Q(hi)*Q(2)**e
def nearest(a):
    q,r=divmod(a.numerator,a.denominator)
    return q+(2*r>a.denominator or (2*r==a.denominator and q%2))
def key(nums,den):return tuple(str(Q(n,den))for n in nums)

def check():
    tic=time.monotonic();reg=read(HERE/'registration.json');run=read(HERE/'run_receipt.json');res=read(HERE/'results.json')
    assert reg['g0']=='1/65000' and reg['max_native']==3 and reg['total_timeout_seconds']==600
    for p,h in reg['source_sha256'].items():assert sha(p)==h,p
    assert sha(HERE/'run_batch.py')==reg['driver_sha256']
    assert sha(HERE/'registration.md')==reg['registration_markdown_sha256']
    assert run['registration_sha256']==sha(HERE/'registration.json')
    assert run['results_sha256']==sha(HERE/'results.json')
    assert run['status']=='COMPLETED' and run['exit_code']==0 and run['all_child_handles_closed']
    assert run['start_epoch']+run['elapsed_seconds']<reg['deadline_epoch']
    assert not run['starting_point_recomputed'] and run['new_H_K_p']==0 and not run['new_C6_bound']
    assert res['status']=='THREE_NEW_POINTS_COMPLETED_G0_STILL_UNRESOLVED'
    assert len(run['native_calls'])==len(res['points'])==3
    assert not(HERE/'candidate_04.json').exists() and not(HERE/'native_outputs/04.stdout.json').exists()
    _,rays=prep.fan_data();forms=prep.all_zero_forms(rays)
    assert len(forms)==9308 and hashlib.sha256(json.dumps(forms,separators=(',',':')).encode()).hexdigest()==reg['zero_forms_sha256']
    seen=set(map(tuple,reg['old_G_exact_point_keys']));assert len(seen)==159
    nums,den=reg['starting_numerators'],reg['starting_denominator']
    assert key(nums,den)in seen and den==10**14
    oldpath=Path(reg['starting_native_path']);old=read(oldpath)
    assert sha(oldpath)==reg['starting_native_sha256']
    summary=[]
    for i,(rec,row)in enumerate(zip(run['native_calls'],res['points']),1):
        c=read(HERE/rec['candidate_path']);out=read(HERE/rec['stdout'])
        assert rec['index']==row['index']==c['step_index']==i
        assert rec['status']=='COMPLETE_VALIDATED' and rec['exit_code']==0 and rec['elapsed_seconds']<120
        assert sha(HERE/rec['candidate_path'])==rec['candidate_sha256']
        assert sha(HERE/rec['stdout'])==rec['stdout_sha256'] and sha(HERE/rec['stderr'])==rec['stderr_sha256']
        assert (HERE/rec['stderr']).read_bytes()==b''
        assert c['starting_numerators']==nums and c['starting_denominator']==den
        assert c['starting_native_path']==str(oldpath) and c['starting_native_sha256']==sha(oldpath)
        mid=[sum(interval(v))/2 for v in old['gradient']]
        assert list(map(Q,c['gradient_exact_dyadic_midpoints']))==mid
        maximum=max(map(abs,mid[1:]));assert maximum>0
        direction=[Q(0)]+[-v/maximum for v in mid[1:]]
        raw=[Q(n,den)+v/20 for n,v in zip(nums,direction)]
        assert Q(c['normalization_max_free_abs'])==maximum and c['fixed_step']=='1/20'
        assert list(map(Q,c['negative_normalized_direction']))==direction
        assert list(map(Q,c['raw_rational_proposal']))==raw
        gap=Q(1,10**8);clip=[Q(1)]
        for j in range(1,6):
            value=max((6-j)*gap,min(1-j*gap,raw[j]))
            clip.append(min(value,clip[-1]-gap))
        assert list(map(Q,c['clipped_rational_proposal']))==clip
        rounded=[nearest(x*den)for x in clip]
        assert rounded==c['quantized_numerators']
        assert list(map(Q,c['quantization_errors']))==[Q(n,den)-v for n,v in zip(rounded,clip)]
        pf0=prep.preflight(rounded,den,forms,rays);assert pf0==c['initial_preflight']
        assert pf0['feasible'] and pf0['int64_safe'] and pf0['zero_form_count']>0
        # Every actual candidate uses the first registered gap repair, reconstructed here.
        assert len(c['repairs'])==1 and c['repairs'][0]['k']==1
        gaps=[rounded[j]-rounded[j+1]for j in range(5)]+[rounded[5]]
        largest=max(range(6),key=lambda j:gaps[j]);delta=[11**j if j!=largest else 0 for j in range(6)];delta[largest]=-sum(delta)
        moved=[a+b for a,b in zip(gaps,delta)]
        newnums=[den]+[den-sum(moved[:j])for j in range(1,6)]
        repair=c['repairs'][0]
        assert repair['largest_gap_index']==largest and repair['gap_delta_numerators']==delta
        assert repair['base_gap_numerators']==gaps and repair['candidate_gap_numerators']==moved
        assert newnums==repair['numerators']==c['final_numerators']==rec['numerators']==row['numerators']
        assert rec['denominator']==row['denominator']==den and key(newnums,den)not in seen
        if i==1:assert newnums==reg['first_candidate_numerators']
        pf=prep.preflight(newnums,den,forms,rays)
        assert pf==c['final_preflight']==repair['preflight']
        assert pf['feasible'] and pf['int64_safe'] and pf['zero_form_count']==0
        assert Q(c['clipped_step_gradient_dot_displacement'])==sum(mid[j]*Q(newnums[j]-nums[j],den)for j in range(1,6))
        assert rec['command']==[reg['binary_path'],'full',reg['fan_path'],str(den),*map(str,newnums)]
        assert out['status']=='generic_nonzero_supported' and out['precision_bits']==512
        assert out['full_group_terms']==967680 and out['computed_sign_pairs']==483840
        assert out['zero_node_terms']==out['confluent_terms']==0 and len(out['gradient'])==6
        g=interval(out['G']);gradient=[interval(v)for v in out['gradient']]
        assert g[0]>Q(1,65000)
        for a,b in [g,*gradient]:assert b-a<=max(Q(1),abs((a+b)/2))*Q(2)**-256
        native_euler=interval(out['Euler_x_grad_plus5G']);assert native_euler[0]<=0<=native_euler[1]
        lo,hi=5*g[0],5*g[1]
        for n,(a,b)in zip(newnums,gradient):lo+=Q(n,den)*a;hi+=Q(n,den)*b
        assert lo<=0<=hi
        saved=row['check']
        assert list(map(Q,saved['G_interval']))==list(g)
        assert [list(map(Q,p))for p in saved['gradient_intervals']]==[list(p)for p in gradient]
        assert list(map(Q,saved['independent_Euler_interval']))==[lo,hi]
        assert saved['strictly_above_g0'] and not saved['rejects_g0_strictly']
        assert Q(saved['G_upper_minus_g0'])==g[1]-Q(1,65000)
        oldg=interval(old['G']);assert g[1]<oldg[0] and row['strictly_below_previous']
        assert list(map(Q,row['G_minus_previous_interval']))==[g[0]-oldg[1],g[1]-oldg[0]]
        x=[Q(n,den)for n in newnums]
        distance=max(min(abs(x[a]-x[b])for a,b in itertools.combinations(t,2))/2 for t in itertools.combinations(range(6),3))
        assert distance==Q(saved['exact_dist_infinity_to_S2_lower_bound'])>0
        rational_lo=Q((g[0]*10**20).__floor__(),10**20);rational_hi=Q((g[1]*10**20).__ceil__(),10**20)
        assert Q(1,65000)<rational_lo<g[0]<=g[1]<rational_hi
        summary.append(dict(index=i,numerators=newnums,denominator=den,G_simple_rational_bounds=[str(rational_lo),str(rational_hi)],strictly_above_target=True,strictly_below_previous=True,clipped_free_coordinates=[j for j in range(1,6)if raw[j]!=clip[j]],zero_forms_before_repair=pf0['zero_form_count'],repair_k=1,min_abs_primitive_form=pf['min_abs_primitive_form'],distance_to_true_S2_lower_bound=str(distance)))
        seen.add(key(newnums,den));nums=newnums;oldpath=HERE/rec['stdout'];old=out
    return dict(status='PASS_THREE_NEW_POINTS_TARGET_UNRESOLVED',checker_pid=os.getpid(),elapsed_seconds=time.monotonic()-tic,original_native_calls=3,new_native_calls=0,seed_recomputed=False,full_terms_per_point=967680,zero_forms_checked_per_point=9308,source_hashes_unchanged=True,independent_proposal_and_Euler_reconstruction=True,points=summary,global_g_lower_bound_proved=False,new_C6_bound=False,checker_source_sha256=sha(__file__))

if __name__=='__main__':
    assert sys.argv[1:] in ([],['--write-receipt'])
    result=check()
    if sys.argv[1:]==['--write-receipt']:
        assert not(HERE/'checked_results.json').exists()
        (HERE/'checked_results.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))
