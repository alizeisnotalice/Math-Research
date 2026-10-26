#!/usr/bin/env python3
"""One finite rational-only saved-data audit; no Arb or kernel calls."""
import sys
sys.dont_write_bytecode=True
sys.set_int_max_str_digits(0)
from pathlib import Path
from fractions import Fraction as Q
import hashlib,json,time,math
H=Path(__file__).resolve().parent;U=H.parents[1]
def rd(p):return json.loads(Path(p).read_text())
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def b(r):return Q(int(r['lower_mantissa']))*Q(2)**int(r['lower_binary_exponent']),Q(int(r['upper_mantissa']))*Q(2)**int(r['upper_binary_exponent'])
def ceilq(x,d):return Q(math.ceil(x*d),d)
def floorq(x,d):return Q(math.floor(x*d),d)
def prim(z,N=96):
    assert 0<z<1
    s=z/(2*(1-z*z))+sum((z**(2*j+1)/Q(2*(2*j+1)) for j in range(N)),Q(0))
    e=z**(2*N+1)/Q(2*(2*N+1))/(1-z*z)
    return s,s+e

def main():
    tic=time.monotonic();reg=rd(H/'registration.json');r=rd(H/'results.json');rc=rd(H/'run_receipt.json')
    assert rc['status']=='COMPLETED' and rc['completed_configurations']==rc['attempted_configurations']==7 and rc['exit_code']==0
    assert rc['registration_sha256']==sha(H/'registration.json') and rc['results_sha256']==sha(H/'results.json')
    for p,h in reg['source_sha256'].items():assert sha(p)==h,p
    for p,h in reg['code_sha256'].items():assert sha(H/p)==h,p
    assert len(r['configurations'])==7
    weights=[math.comb(5,k-1) for k in range(1,6)]
    assert weights==r['union']['weights']==[1,5,10,10,5] and sum(weights)==31
    assert r['union']['overlap_permitted'] and not r['union']['domains_claimed_disjoint']
    rho=Q(1,64);eps=Q(1,1000)
    for row in r['configurations'][0]['cases']:
        k=row['k'];m=6-k
        assert row['m']==m and row['pair_Gram']==4*k*m
        expected={'c0v':Q(1,2)-rho,'c0w':Q(15,16)-rho,'c1v_squared_reduced':Q(m,96)-rho*rho,'c1w_squared_reduced':Q(m,384)-rho*rho,'g_squared':Q(1,32)-k*rho*rho}
        assert {a:Q(c)for a,c in row['positive_gate_margins'].items()}==expected and all(x>0 for x in expected.values())
    old43=rd(U/'round043/two_parameter_middle/results.json')['configurations']
    old49=rd(U/'round049/middle_joint_marginal/results.json')['configurations'][2]
    k3=r['cases']['3']
    for a,c in [('short_big_upper','new_big_short_upper'),('short_bad_upper','new_bad_short_upper'),('short_total_upper','new_total_short_upper'),('C_joint','C_joint')]:assert k3[a]==old49[c]
    assert k3['D_big']==old43[3]['C2'] and k3['M_combined']==old43[4]['M_combined']
    assert k3['source043_sha256']==sha(U/'round043/two_parameter_middle/results.json') and k3['source049_sha256']==sha(U/'round049/middle_joint_marginal/results.json')
    peaks={x['id']:x for x in rd(U/'round033/exact_projection_peaks/results.json')['cases']}
    piL=floorq(b(r['shared_constants']['pi'])[0],10**30)
    assert Q(3141592653589793238462643383279,10**30)<=piL<Q(22,7)
    Fup=ceilq(b(r['shared_constants']['cauchy_average_factor'])[1],10**30)
    assert (Fup-1)**2>=Q(6,928)
    cert={3:b(k3['short_total_upper'])[1]};checks=[]
    for k in [1,2,4,5]:
        row=r['cases'][str(k)];m=6-k
        assert row['n']==weights[k-1] and row['m']==m and row['all_complementary_events_paid']
        R=Q(peaks[f'r2_k{min(k,m)}']['R_exact']);assert R==Q(row['R2_exact'])
        dl,du=b(row['joint_coordinate_deviation']);dev2=Q(m-1,m)*rho*rho
        assert 0<=dl<=du and dl*dl<=dev2<=du*du
        if m==1:assert dl==du==0
        devU=ceilq(du,10**30);L=Q(1,4)-devU;V=Q(3,4)+devU
        assert 0<L<V<1
        Clow,Cup=prim(L);Vlow,Vup=prim(V)
        C=ceilq(128*(Vup-Clow)/piL**2,10**20)
        assert C>=b(row['C_joint'])[1]
        etaU=ceilq(b(row['eta'])[1],10**30)
        assert etaU*etaU>=k*rho*rho
        gaps={}
        target={'hv0':Q(k,4),'hv1':Q(k*m,96),'hw0':Q(225*k,256),'hw1':Q(k*m,384),'g':Q(1,32)}
        for a,t in target.items():
            low=floorq(b(row[a])[0],10**30)-Q(4,10**30)
            assert low>0 and (low+etaU)**2<=t,(k,a)
            gaps[a]=low
        D=64*R/piL**6
        common=Q(1000,729)*Q(44743079,2970)/piL**6*Fup/gaps['g']**2
        M=common*(1/(gaps['hv0']**2*gaps['hv1']**2)+1/(gaps['hw0']**2*gaps['hw1']**2))
        inflation=1+Q(1,10**18)
        good=b(row['short_big_upper'])[1]*inflation;bad=b(row['short_bad_upper'])[1]*inflation
        assert good**6>=Q(3,7)**6*D*C**5*eps**14
        assert bad**6>=Q(3,8)**6*M*C**5*eps**16
        cert[k]=good+bad
        assert Q(row['simple_rational_case_upper'])>cert[k]
        checks.append(dict(k=k,m=m,C_joint_rational_upper=str(C),simple_case_upper=row['simple_rational_case_upper'],sqrt_gap_checks=True,log_series_terms=96,log_positive_remainder=True,sixth_power_budget_checks=True))
    S=Q(0)
    for term in r['union']['weighted_terms']:
        k=term['k'];n=weights[k-1];assert term['n']==n
        upper=b(term['weighted_power'])[1]*(1+Q(1,10**16))
        assert (upper/n)**5>=cert[k]**6
        S+=upper
    final=Q(r['union']['simple_rational_union_upper'])
    assert final**6>S**5 and final<Q(r['union']['conditional_Bcrit'])
    result=dict(status='PASS_RATIONAL_SAVED_DATA_AUDIT',input_results_sha256=sha(H/'results.json'),registration_sha256=sha(H/'registration.json'),producer_configurations=7,saved_data_audit_configurations=1,total_configurations=8,new_Arb_or_kernel_configurations=0,new_kernel_calls=0,old_k3_scalar_recomputations=0,all_source_hashes_unchanged=True,k3_saved_intervals_identical=True,new_types=checks,union_weight_checks=weights,union_power_check='simple_upper^6 > (sum certified weighted powers)^5',simple_rational_union_upper=str(final),elapsed_seconds=time.monotonic()-tic)
    (H/'checked_results.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))
if __name__=='__main__':main()
