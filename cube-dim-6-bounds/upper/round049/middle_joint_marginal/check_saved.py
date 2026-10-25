#!/usr/bin/env python3
"""Saved-data rational certificate: no Arb, kernel, or new parameter choice."""
import sys
sys.dont_write_bytecode=True
sys.set_int_max_str_digits(0)
from pathlib import Path
from fractions import Fraction as Q
import hashlib,json,time,os
HERE=Path(__file__).resolve().parent
OLD=HERE.parents[1]/'round043/two_parameter_middle'
def read(p):return json.loads(Path(p).read_text())
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def bounds(r):
    lo=Q(int(r['lower_mantissa']))*Q(2)**int(r['lower_binary_exponent'])
    hi=Q(int(r['upper_mantissa']))*Q(2)**int(r['upper_binary_exponent'])
    assert lo<=hi and hi-lo<=max(Q(1),abs((lo+hi)/2))*Q(2)**-450
    return lo,hi
def ceil_at(x,n):return Q((x*10**n).__ceil__(),10**n)
def floor_at(x,n):return Q((x*10**n).__floor__(),10**n)
def F_bounds(z,N=96):
    assert 0<z<1
    acc=Q(0);term=z;z2=z*z
    for k in range(N):
        acc+=term/(2*k+1);term*=z2
    tail=term/((2*N+1)*(1-z2))
    rational=z/(2*(1-z2))
    return rational+acc/2,rational+(acc+tail)/2
def determinant(m):
    a=[[Q(x)for x in row]for row in m];n=len(a);value=Q(1)
    for j in range(n):
        k=next(k for k in range(j,n)if a[k][j])
        if k!=j:a[j],a[k]=a[k],a[j];value=-value
        pivot=a[j][j];value*=pivot
        for i in range(j+1,n):
            ratio=a[i][j]/pivot
            for k in range(j+1,n):a[i][k]-=ratio*a[j][k]
            a[i][j]=0
    return value
def check():
    tic=time.monotonic();reg=read(HERE/'registration.json');run=read(HERE/'run_receipt.json');res=read(HERE/'results.json')
    for p,h in reg['source_sha256'].items():assert sha(p)==h,p
    for f,h in reg['code_sha256'].items():assert sha(HERE/f)==h,f
    assert run['status']=='COMPLETED' and run['exit_code']==0 and run['completed_configurations']==run['attempted_configurations']==3
    assert run['new_kernel_calls']==run['full_enumerations']==run['matlab_calls']==0
    assert run['registration_sha256']==sha(HERE/'registration.json') and run['results_sha256']==sha(HERE/'results.json')
    rows=res['configurations'];assert [r['id']for r in rows]==[1,2,3]
    a=rows[0];columns=a['basis_columns'];matrix=list(map(list,zip(*columns)))
    assert abs(determinant(matrix))==32
    assert abs(determinant([[columns[j][i]for j in range(2)]for i in[0,3]]))==2
    assert a['free_schur_abs_det']==16 and a['block_counts']==[2**4,2**4]
    assert Q(a['inverse_joint_Jacobian'])==Q(1,2) and Q(a['remaining_positive_orthant_factor'])==Q(1,16)
    assert a['two_density_numerator']==Q(1,2)*16**2*4**2==2048
    assert a['joint_L1_numerator']==2048//16==128
    assert Q(a['B_deviation_functional_norm_squared'])==Q(4,9)+2*Q(1,9)==Q(2,3)
    b=rows[1];assert b['branch_symmetry_assumed'] is False and b['actual_disjoint_complement_domination']
    assert -2*Q(1,6)+2*Q(5,6)==Q(4,3) and 2*Q(5,6)==Q(5,3)
    assert b['integrated_prefactors']==['3/7','3/8'] and b['integrated_time_powers']==['7/3','8/3']
    r=rows[2];old=read(OLD/'results.json')['configurations']
    assert sha(OLD/'results.json')==r['source_results_sha256']
    assert r['source_D_big_interval']==old[3]['C2'] and r['source_M_bad_interval']==old[4]['M_combined'] and r['source_old_budget_interval']==old[5]['total_short_norm_upper']
    assert not r['old_budget_recomputed'] and r['source_values_from_complete_dyadic_endpoints']
    eps=Q(r['epsilon']);rho=Q(r['rho']);assert eps==Q(1,1000) and rho==Q(1,64)
    deltaL,deltaU=bounds(r['coordinate_deviation']);assert 0<deltaL and deltaL**2<=Q(2,3)*rho*rho<=deltaU**2
    deltaU=ceil_at(deltaU,30);L=Q(1,4)-deltaU;U=Q(3,4)+deltaU;assert 0<L<U<1
    piL,piU=bounds(r['pi']);piL=floor_at(piL,30);piU=ceil_at(piU,30)
    assert 3<piL<piU<Q(22,7)
    # Independent positive atanh series,96 terms and a proved geometric tail.
    # This verifies the same fixed integral; it is not a new range or Arb call.
    Flo=F_bounds(L)[0];Fhi=F_bounds(U)[1]
    C_upper=ceil_at(128*(Fhi-Flo)/piL**2,20)
    assert bounds(r['C_joint'])[1]<C_upper
    assert C_upper*piU*eps<1
    D_upper=ceil_at(bounds(r['source_D_big_interval'])[1],30)
    M_upper=ceil_at(bounds(r['source_M_bad_interval'])[1],30)
    enlargement=1+Q(1,10**18)
    good=bounds(r['new_big_short_upper'])[1]*enlargement
    bad=bounds(r['new_bad_short_upper'])[1]*enlargement
    assert good**6>=Q(3,7)**6*D_upper*C_upper**5*eps**14
    assert bad**6>=Q(3,8)**6*M_upper*C_upper**5*eps**16
    simple=Q(r['simple_rational_upper']);assert good+bad<simple<Q(4,48750)
    assert bounds(r['new_total_short_upper'])[1]<simple
    oldlo,oldhi=bounds(r['source_old_budget_interval']);assert oldlo>18*simple
    return dict(status='PASS_TRUE_JOINT_MARGINAL_SAVED_RATIONAL_CERTIFICATE',checker_pid=os.getpid(),elapsed_seconds=time.monotonic()-tic,new_Arb_configurations=0,new_kernel_calls=0,new_parameter_choices=0,atanh_series_terms_per_endpoint=96,positive_tail_accounted=True,C_joint_rational_upper=str(C_upper),short_middle_norm_rational_upper=str(simple),saved_old_upper_budget_greater_than_18_times_new_bound=True,old_budget_formula_recomputed=False,branch_symmetry_assumed=False,new_I_certificate=False,new_C6_bound=False,checker_source_sha256=sha(__file__))
if __name__=='__main__':
    assert sys.argv[1:]in([],['--write-receipt'])
    result=check()
    if sys.argv[1:]==['--write-receipt']:
        assert not(HERE/'checked_results.json').exists()
        (HERE/'checked_results.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))
