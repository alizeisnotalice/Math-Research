#!/usr/bin/env python3
"""B036 pure exact saved-data check. No Arb configuration or kernel is run."""
import sys
sys.dont_write_bytecode=True
sys.set_int_max_str_digits(0)
import hashlib,json,os,time
from fractions import Fraction as Q
from pathlib import Path
HERE=Path(__file__).resolve().parent
def read(p):return json.loads(Path(p).read_text())
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def bounds(row):
    a,b=Q(int(row['lower_mantissa']))*Q(2)**int(row['lower_binary_exponent']),Q(int(row['upper_mantissa']))*Q(2)**int(row['upper_binary_exponent'])
    assert a<=b and abs(int(row['lower_binary_exponent']))<10000 and abs(int(row['upper_binary_exponent']))<10000
    return a,b
def check():
    tic=time.perf_counter();reg=read(HERE/'registration.json');rel=read(HERE/'release.json');rec=read(HERE/'run_receipt.json');out=read(HERE/'results.json')
    for p,h in reg['source_sha256'].items():assert sha(p)==h
    for p,h in rel['gate_source_sha256'].items():assert sha(p)==h
    assert rel['explicit_release'] and rel['F035_PASS_or_main_release']
    assert rel['registration_sha256']==rec['registration_sha256']==sha(HERE/'registration.json')
    assert rec['release_sha256']==sha(HERE/'release.json')
    for n,h in reg['code_sha256'].items():assert sha(HERE/n)==h
    if(HERE/'manifest.json').exists():
        for n,h in read(HERE/'manifest.json')['whitelist_sha256'].items():assert sha(HERE/n)==h
    assert rec['status']=='COMPLETED' and rec['exit_code']==0 and rec['elapsed_seconds']<1200
    assert rec['attempted_configurations']==rec['completed_configurations']==len(out['configurations'])==12
    assert rec['results_sha256']==sha(HERE/'results.json')
    assert all(rec[k]==0 for k in ['spatial_kernel_calls','grid_points','native_G_H_K_calls','matlab_calls'])
    rows=out['configurations'];assert [r['index']for r in rows]==list(range(1,13))
    assert [r['name']for r in rows]==reg['configuration_names']
    peaks={r['id']:r for r in read(reg['peaks_path'])['cases']}
    r1=Q(peaks['r1']['R_exact']);first=rows[0]
    assert Q(first['R1'])==r1==Q(67320683,3251298631680)
    assert Q(first['new_over_old_Cf_exact'])==7776*r1<1
    pl,pu=bounds(first['pi']);assert Q(3141592653,10**9)<pl<=pu<Q(3141592654,10**9)
    cflo,cfhi=bounds(first['Cf']);assert 0<cflo<=cfhi
    assert cflo*pl**6<=32*r1<=cfhi*pu**6
    qlo,qhi=bounds(first['Q5']);assert qlo*qlo*(2*pl)**10<=6*r1*r1<=qhi*qhi*(2*pu)**10
    weighted=Q(0);counts=[480,1920,480,640,1440];dets=[64,96,128,160,192]
    for row,n,d in zip(rows[1:6],counts,dets):
        original=peaks[row['name']];r=Q(original['R_exact'])
        assert Q(row['R3'])==r and row['unordered_count']==original['triple_orbit_count']==n
        assert row['Gram_determinant']==original['Gram_determinant']==d and row['ordered_count']==6*n
        lo,hi=bounds(row['Q3_over_sqrt_Gram']);assert lo*(2*pl)**3<=r<=hi*(2*pu)**3
        weighted+=n*r
    aggregate=rows[6];assert aggregate['unordered_total']==sum(counts)==4960 and aggregate['ordered_total']==29760
    assert Q(aggregate['weighted_R3_sum'])==weighted and Q(aggregate['B3_pi6_exact'])==24576*weighted
    blo,bhi=bounds(aggregate['B3']);assert blo*pl**6<=24576*weighted<=bhi*pu**6
    geometry=rows[7];assert geometry['delta']=='1/16'
    vlo,vhi=bounds(geometry['volume_upper']);assert 0<vlo<=vhi<1
    old=read(reg['geometry_results_path'])['rows'][0]['new_volume']
    assert max(vlo,Q(old['lower_exact']))<=min(vhi,Q(old['upper_exact']))
    tail=rows[8];assert tail['epsilon']=='1/1000';tl,tu=bounds(tail['tail_upper_budget']);assert tl>0
    off=rows[9];assert off['delta']=='1/16' and off['epsilon']=='1/1000';assert bounds(off['off_S2_tail_budget'])[0]>0
    cost=rows[10];assert cost['g_target']=='1/65000' and cost['H_norm_goal']=='1/14000' and cost['B_critical_exact']=='4/48750'
    gl,gu=bounds(cost['B_goal']);cl,cu=bounds(cost['conditional_cost']);assert gl<=Q(1,14000)+tu and gu>=Q(1,14000)+tl
    assert cl**5<=(6*gu*65000)**6 and cu**5>=(6*gl*65000)**6
    flag=(gu<Q(4,48750) and cu<64)
    assert flag==cost['strictly_below64_if_both_global_premises_proved']==out['conditional_cost_below64']
    assert not cost['both_global_premises_proved']
    tube=rows[11];assert tube['delta']=='1/16' and tube['derived_epsilon_selection_rule']=='(Cf*pi*V)^(1/5)'
    el,eu=bounds(tube['epsilon_star']);assert el**5<=cfhi*pu*vhi and eu**5>=cflo*pl*vlo
    assert bounds(tube['R'])[0]>0
    ul,uu=bounds(tube['allocated_upper_budget']);ll,lu=bounds(tube['best_two_envelope_lower_budget'])
    assert 0<ll<=lu<ul<=uu
    fails=ll>Q(4,48750);assert fails==tube['two_envelope_method_fails_new_target']
    assert Q(tube['safe_lower_method_ratio'])<=ll/Q(4,48750)
    assert tube['old_cap_formula_stopped'] and not tube['actual_tube_G_norm_lower_bound']
    assert not out['new_I_certificate'] and not out['new_global_g_lower_bound'] and not out['new_global_H_norm_upper_bound'] and not out['new_C6_bound']
    return dict(status='PASS_EXACT_SAVED_PEAK_BUDGET_CHECK',pid=os.getpid(),elapsed_seconds=time.perf_counter()-tic,
                source_and_release_hashes_valid=True,ordered_triples=29760,configurations_in_original_batch=12,
                new_Arb_configurations_in_checker=0,physical_kernel_calls=0,
                tube_method_fails_new_target=fails,safe_lower_method_ratio=tube['safe_lower_method_ratio'],
                conditional_cost_below64_if_both_global_premises_proved=flag,
                exact_B3_pi6_coefficient=str(24576*weighted),new_I_certificate=False,new_C6_bound=False)

if __name__=='__main__':
    assert sys.argv[1:]in([],['--write-receipt'])
    result=check()
    if sys.argv[1:]==['--write-receipt']:
        assert not(HERE/'checked_results.json').exists()
        (HERE/'checked_results.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result),flush=True)
