#!/usr/bin/env python3
"""One pure-rational audit; fresh Machin bounds and saved algebraic candidates."""
import sys
sys.dont_write_bytecode=True
sys.set_int_max_str_digits(0)
from pathlib import Path
from fractions import Fraction as Q
import json,hashlib,time
H=Path(__file__).resolve().parent
def rd(p):return json.loads(Path(p).read_text())
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def bounds(r):
    if 'common_binary_exponent'in r:
        e=int(r['common_binary_exponent']);return Q(int(r['lower_mantissa']))*Q(2)**e,Q(int(r['upper_mantissa']))*Q(2)**e
    return Q(int(r['lower_mantissa']))*Q(2)**int(r['lower_binary_exponent']),Q(int(r['upper_mantissa']))*Q(2)**int(r['upper_binary_exponent'])
def wide(row):
    l,u=bounds(row);eps=Q(1,10**40);return l*(1-eps),u*(1+eps)
def atan_bounds(d,n):
    assert n%2==0
    l=sum((Q((-1)**k,(2*k+1)*d**(2*k+1))for k in range(n)),Q(0))
    return l,l+Q(1,(2*n+1)*d**(2*n+1))
def main():
    tic=time.monotonic();assert not(H/'checked_results.json').exists(),'one saved-data check only'
    reg=rd(H/'registration.json');binding=rd(H/'execution_binding.json');receipt=rd(H/'run_receipt.json');r=rd(H/'results.json')
    assert receipt['status']=='COMPLETED'and receipt['completed_configurations']==receipt['attempted_configurations']==9 and receipt['exit_code']==0
    assert time.time()<receipt['deadline_epoch']
    assert receipt['registration_sha256']==binding['registration_sha256']==sha(H/'registration.json')and receipt['binding_sha256']==sha(H/'execution_binding.json')and receipt['results_sha256']==sha(H/'results.json')
    for f,x in reg['source_sha256'].items():assert sha(f)==x,f
    for f,x in binding['code_sha256'].items():assert sha(H/f)==x,f
    old_path=next(p for p in reg['source_sha256']if p.endswith('/round037/pair_resolved_green_envelope/results.json'))
    old=rd(old_path);raw={str(x['k']):x['M']for x in old['configurations'][:3]};assert raw==r['source_M_intervals']
    assert r['configurations'][0]['background_counts']==[180,450,300]and r['configurations'][0]['ones_first_counts']==[6,15,10]
    aa,bb=atan_bounds(5,120),atan_bounds(239,40);piL,piU=16*aa[0]-4*bb[1],16*aa[1]-4*bb[0];assert 3<piL<piU<4
    c=r['constants'];sL,sU=wide(c['s0']);RL,RU=wide(c['R']);zL,zU=wide(c['sqrt6']);tL,tU=wide(c['sqrt2'])
    rho=Q(1,64);a=Q(47,64);R2=Q(385,4096);r02=1+5*a*a;D2=Q(10,3)*a*a
    assert sL**4<=31*rho**4<=sU**4 and RL*RL<=R2<=RU*RU and zL*zL<=6<=zU*zU and tL*tL<=2<=tU*tU
    rbL,rbU=wide(c['rho_power_7_over5']);sbL,sbU=wide(c['s0_power_7_over5']);RbL,RbU=wide(c['R_power_7_over5'])
    assert rbL**5<=rho**7<=rbU**5
    assert sbL**5<=sL**7 and sbU**5>=sU**7
    assert RbL**5<=RL**7 and RbU**5>=RU**7
    sphereL,sphereU=Q(8,3)*piL*piL,Q(8,3)*piU*piU
    KL,KU=sphereL/zU*Q(10,7),sphereU/zL*Q(10,7)
    J3L=KL*Q(5,7)*(2*sbL-rbU-sbU*sbU/RbL)
    J3U=KU*Q(5,7)*(2*sbU-rbL-sbL*sbL/RbU)
    J1U=KU*(sbU*sbU*RU-Q(14,19)*sbL*sbL*sL-Q(5,19)*rbL*rbL*rho)
    VU=sphereU/zL*(31*rho**4*(RU-Q(4,5)*sL)-rho**5/5)
    assert 0<J3L<J3U and J1U>0 and VU>0
    n=[6,15,10];WL=sum((n[k-1]*bounds(raw[str(k)])[0]for k in range(1,4)),Q(0));WU=sum((n[k-1]*bounds(raw[str(k)])[1]for k in range(1,4)),Q(0))
    coefs=[WU/r02,WU/(r02*D2),30*WU/(r02*D2*a*tL)];integrals=[J3U,J1U,VU];names=['ones_first_joint','reverse_order','all930_background'];certs=[]
    for name,coef,j in zip(names,coefs,integrals):
        v=wide(r['branches'][name]['norm_upper'])[1]
        assert v**6>coef**6*j**5
        certs.append(v)
    final=Q(r['union']['simple_rational_upper']);assert sum(certs,Q(0))<final
    target=Q(4,48750);first_exceeds=(WL/r02)**6*J3L**5>target**6
    assert first_exceeds and r['union']['comparison']=='FIXED_UPPER_FORMULA_EXCEEDS_TARGET'
    assert r['union']['all992_terms_paid']and r['union']['all31_plane_overlap_retained']and not r['union']['actual_G_lower_bound']and not r['union']['original_Rc_excluded']
    finished=time.time();assert finished<receipt['deadline_epoch']
    out=dict(status='PASS_PURE_RATIONAL_HIGH_BUDGET_AUDIT',results_sha256=sha(H/'results.json'),registration_sha256=sha(H/'registration.json'),binding_sha256=sha(H/'execution_binding.json'),producer_configurations=9,pure_saved_data_configurations=1,total_configurations=10,exact_Machin_terms=[120,40],Machin_pi_enclosure_proved=True,algebraic_root_power_checks=7,complete_M_source_intervals_identical=True,all_three_branch_sixth_power_checks=True,strict_rational_high_upper=str(final),first_branch_formula_alone_strictly_exceeds_target=first_exceeds,comparison_direction='lower bound on this fixed upper expression only; no lower bound on actual G',conditional_target=str(target),old_M_recomputations=0,new_Arb_calls=0,new_kernel_calls=0,full_enumerations=0,matlab_calls=0,new_global_I=False,new_C6_bound=False,elapsed_seconds=time.monotonic()-tic,finished_epoch=finished,total_wall_from_producer_seconds=finished-receipt['started_epoch'])
    (H/'checked_results.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
if __name__=='__main__':main()
