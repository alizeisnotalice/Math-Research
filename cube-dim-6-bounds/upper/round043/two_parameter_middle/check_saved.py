#!/usr/bin/env python3
"""Exact saved-data verification only; no new Arb or spatial configurations."""
import sys
sys.dont_write_bytecode=True
sys.set_int_max_str_digits(0)
from pathlib import Path
from fractions import Fraction as Q
from itertools import product
import json,hashlib,time,os
HERE=Path(__file__).resolve().parent
def read(p):return json.loads(Path(p).read_text())
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def bounds(r):
    lo=Q(int(r['lower_mantissa']))*Q(2)**int(r['lower_binary_exponent'])
    hi=Q(int(r['upper_mantissa']))*Q(2)**int(r['upper_binary_exponent'])
    assert lo<=hi and hi-lo<=max(Q(1),abs((lo+hi)/2))*Q(2)**-450
    return lo,hi
def dot(x,y):return sum(a*b for a,b in zip(x,y))
def check():
    tic=time.monotonic();reg=read(HERE/'registration.json');run=read(HERE/'run_receipt.json');res=read(HERE/'results.json')
    for p,h in reg['source_sha256'].items():assert sha(p)==h,p
    for f,h in reg['code_sha256'].items():assert sha(HERE/f)==h,f
    assert run['status']=='COMPLETED' and run['exit_code']==0
    assert run['attempted_configurations']==run['completed_configurations']==6
    assert run['registration_sha256']==sha(HERE/'registration.json') and run['results_sha256']==sha(HERE/'results.json')
    assert run['native_calls']==run['spatial_grid_points']==run['matlab_calls']==0
    assert res['status']=='COMPLETE_SIX_CONFIGURATIONS_TRUE_SHORT_MIDDLE_UPPER'
    rows=res['configurations'];assert [r['id']for r in rows]==list(range(1,7))
    assert sum(abs(1+sum(s))for s in product((-1,1),repeat=5))==rows[0]['psi_ones']==60
    assert rows[0]['missing_symbol_fraction']=='9/10' and rows[0]['dimension3_peak_loss']=='1000/729'
    b,c=rows[1]['other_pair'];x=[Q(1)]*3+[Q(1,4)]*3;y=[z-Q(1,4)for z in x]
    assert y==[Q(3,8)*(a+d)for a,d in zip(b,c)]
    assert b not in [(1,)*6,[-1]*6] and len(b)==len(c)==6
    assert dot(b,c)==0 and dot(b,b)==dot(c,c)==6
    assert dot(y,y)-(dot(y,b)**2+dot(y,c)**2)/6==0
    assert dot(x,x)-(dot(x,b)**2+dot(x,c)**2)/6==Q(1,6)
    v=[1]*6;w=[1]*3+[-1]*3
    assert dot(v,v)*dot(w,w)-dot(v,w)**2==rows[1]['selected_pair_Gram_determinant']==36
    assert rows[1]['selected_pair_coarea_J2']==6 and rows[1]['uniform_fiber_distance_claim_refuted']
    geo=rows[2];rho=Q(geo['rho']);eps=Q(geo['epsilon']);assert rho==Q(1,64) and eps==Q(1,1000)
    assert geo['middle_interval']==['1/4','3/4'] and geo['thresholds']==['1/2','1/16']
    piL,piU=bounds(geo['pi']);assert 3<piL<piU<Q(22,7)
    etaL,etaU=bounds(geo['eta']);assert 0<etaL and etaL**2<=3*rho*rho<=etaU**2
    lowers={}
    for name,square in [('hv0',Q(3,4)),('hv1',Q(3,32)),('hw0',Q(675,256)),('hw1',Q(3,128)),('g',Q(1,32))]:
        lo,hi=bounds(geo[name]);lo*=1-Q(2)**-400
        assert 0<lo and (lo+etaU)**2<=square
        lowers[name]=lo
    factorL,factorU=bounds(geo['cauchy_average_factor']);factorU*=1+Q(2)**-400
    assert factorL>1 and (factorU-1)**2>=Q(6,928)
    assert geo['deleted_frame_Gram_bounds']==[26,32] and geo['true_time_coefficients']==[992,960,928]
    R2=Q(rows[3]['R2_exact']);assert R2==Q(2051,47038464)
    assert Q(rows[3]['coefficient_C2_times_pi6'])==64*R2 and rows[3]['quotient_dimension']==4
    coef=Q(rows[4]['Crem_exact_numerator_over_pi6']);assert coef==Q(44743079,2970)
    old=read(HERE.parents[1]/'round041/one_large_jump_cap/results.json')
    assert coef==Q(old['configurations'][1]['Crem_times_pi6_exact'])
    assert rows[4]['ordered_triples_per_retained_direction']==26970
    assert rows[4]['missing_peak_loss_exact']=='1000/729' and rows[4]['omitted_parameter_events']==0 and not rows[4]['Q2_used']
    assert not rows[4]['third_distance_uniformly_positive_assumed'] and rows[4]['actual_Cauchy_density_integrated']
    commonU=Q(1000,729)*coef/piL**6*factorU/lowers['g']**2
    MU=commonU*(1/(lowers['hv0']**2*lowers['hv1']**2)+1/(lowers['hw0']**2*lowers['hw1']**2))
    enlarge=1+Q(2)**-384
    good=bounds(rows[3]['big_branch_short_norm_upper'])[1]*enlarge
    bad=bounds(rows[5]['combined_actual_complement_short_norm_upper'])[1]*enlarge
    assert good**6>=Q(2,3)**6*64*R2/piL**11*eps**9
    assert bad**6>=Q(6,11)**6*MU/piL**5*eps**11
    simple=Q(rows[5]['simple_rational_upper'])
    assert good+bad<simple<Q(4,48750)
    assert bounds(rows[5]['total_short_norm_upper'])[1]<simple
    assert rows[5]['single_actual_complement_L1_budget'] and not rows[5]['whole_face_or_full_time_claim']
    return dict(status='PASS_SAVED_ZERO_GEOMETRY_AND_RATIONAL_SHORT_NORM',checker_pid=os.getpid(),elapsed_seconds=time.monotonic()-tic,new_Arb_configurations=0,new_kernel_calls=0,explicit_fiber_zero_verified=True,positive_first_two_gaps_verified=True,missing_peak_loss_accounted=True,complete_remainders_accounted=True,rational_short_middle_norm_upper=str(simple),new_I_certificate=False,new_C6_bound=False,checker_source_sha256=sha(__file__))
if __name__=='__main__':
    assert sys.argv[1:]in([],['--write-receipt'])
    result=check()
    if sys.argv[1:]==['--write-receipt']:
        assert not(HERE/'checked_results.json').exists()
        (HERE/'checked_results.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))
