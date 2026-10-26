#!/usr/bin/env python3
"""Rational-only saved endpoint audit: no M reevaluation or new Arb calls."""
import sys
sys.dont_write_bytecode=True
sys.set_int_max_str_digits(0)
from pathlib import Path
from fractions import Fraction as Q
from math import comb,ceil,floor
import hashlib,json,time
H=Path(__file__).resolve().parent
OLD=Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube6-audit/cube-dim-6-bounds/upper/round037/pair_resolved_green_envelope')
def rd(p):return json.loads(Path(p).read_text())
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def bounds(r):
    if 'common_binary_exponent'in r:
        e=int(r['common_binary_exponent']);return Q(int(r['lower_mantissa']))*Q(2)**e,Q(int(r['upper_mantissa']))*Q(2)**e
    return Q(int(r['lower_mantissa']))*Q(2)**int(r['lower_binary_exponent']),Q(int(r['upper_mantissa']))*Q(2)**int(r['upper_binary_exponent'])
def ceilq(x,d=10**30):return Q(ceil(x*d),d)
def floorq(x,d=10**30):return Q(floor(x*d),d)
def main():
    tic=time.monotonic();reg=rd(H/'registration.json');r=rd(H/'results.json');receipt=rd(H/'run_receipt.json')
    assert receipt['status']=='COMPLETED' and receipt['completed_configurations']==receipt['attempted_configurations']==7 and receipt['exit_code']==0
    assert receipt['registration_sha256']==sha(H/'registration.json') and receipt['results_sha256']==sha(H/'results.json')
    for p,v in reg['source_sha256'].items():assert sha(p)==v,p
    for p,v in reg['code_sha256'].items():assert sha(H/p)==v,p
    old=rd(OLD/'results.json');oldM={str(x['k']):x['M']for x in old['configurations'][:3]}
    assert oldM==r['source_M_intervals']
    M={int(j):bounds(x)[1]for j,x in oldM.items()};assert all(x>0 for x in M.values())
    N=[192,480,320];weights=[comb(5,k-1)for k in range(1,6)]
    assert weights==[1,5,10,10,5] and sum(weights)==31 and r['union']['weights']==weights
    rho=Q(1,64);piU=ceilq(bounds(r['shared_constants']['pi'])[1]);assert Q(3)<piU<Q(22,7)
    checks=[];cert={}
    for k in range(1,6):
        row=r['cases'][str(k)];struct=r['configurations'][0]['cases'][k-1];m=6-k;j=min(k,m);c=2**m
        bg=[a-(c if ell==j else 0)for ell,a in enumerate(N,1)]
        assert row['n']==weights[k-1] and row['pair_class']==j and row['incident_planes']==2**(m-1) and row['incident_ordered_pairs']==c
        assert row['background_by_class']==struct['background_by_class']==bg and all(a>=0 for a in bg) and row['background_ordered_pairs']==sum(bg)==992-c
        assert Q(struct['A_normal_Gram_determinant'])==Q(1,k)==Q(struct['cone_fraction']) and struct['normal_dimension']==4
        assert Q(struct['line_gap_squared_margin'])==Q(3*m,32)-rho*rho>0
        assert Q(struct['pair_gap_squared_margin'])==Q(9,32)-k*rho*rho>0
        el,eu=bounds(row['eta']);assert el*el<=k*rho*rho<=eu*eu
        etaU=ceilq(eu)
        dL=floorq(bounds(row['d0'])[0])-Q(4,10**30);gL=floorq(bounds(row['g'])[0])-Q(4,10**30)
        assert dL>0 and gL>0 and (dL+etaU)**2<=Q(3*k*m,32) and (gL+etaU)**2<=Q(9,32)
        r2=k*(1-rho)**2;rl,ru=bounds(row['r0']);assert rl*rl<=r2<=ru*ru
        IU=c*M[j]/(r2*dL*dL);BU=sum((a*M[ell]for ell,a in enumerate(bg,1)),Q(0))/(r2*dL*dL*gL)
        fl,fu=bounds(row['face_measure_factor']);assert fl*fl<=Q(m,16*k)<=fu*fu
        fU=ceilq(fu);VU=fU*piU**2/2*rho**4
        # No fractional-power reevaluation: J^5 has integer rho exponent14.
        J5U=fU**5*(Q(5,7)*piU**2)**5*rho**14
        inflation=1+Q(1,10**18)
        sing=bounds(row['singular_norm_upper'])[1]*inflation;regular=bounds(row['background_norm_upper'])[1]*inflation
        assert sing**6>=IU**6*J5U
        assert regular**6>=BU**6*VU**5
        cert[k]=sing+regular
        simple=Q(row['simple_rational_upper']);assert simple>cert[k]
        assert row['all_pairs_paid'] and row['time_range']=='all t>0'
        checks.append(dict(k=k,m=m,incident_ordered_pairs=c,background_by_class=bg,simple_rational_case_upper=str(simple),true_radial_power='J^5 = factor^5*(5*pi^2/7)^5*rho^14',sixth_power_checks=True,source_M_unchanged=True))
    S=Q(0)
    for t in r['union']['weighted_terms']:
        k=t['k'];n=weights[k-1];assert t['n']==n
        u=bounds(t['weighted_power'])[1]*(1+Q(1,10**16))
        assert (u/n)**5>=cert[k]**6
        S+=u
    final=Q(r['union']['simple_rational_union_upper']);target=Q(r['union']['conditional_Bcrit'])
    assert final**6>S**5 and final<Q(251,1000)*target
    assert r['union']['overlap_permitted']and not r['union']['domains_claimed_disjoint']and not r['union']['actual_G_lower_bound']
    out=dict(status='PASS_RATIONAL_SAVED_DATA_AUDIT',input_results_sha256=sha(H/'results.json'),registration_sha256=sha(H/'registration.json'),total_scientific_configurations=8,producer_configurations=7,saved_data_audit_configurations=1,new_Arb_or_kernel_calls=0,old_M_recomputations=0,source_hashes_unchanged=True,cases=checks,weights=weights,domain_count=31,simple_rational_union_upper=str(final),union_power_check='final^6 > (sum certified weighted powers)^5',strict_relative_upper='251/1000',conditional_target=str(target),new_C6_bound=False,elapsed_seconds=time.monotonic()-tic)
    (H/'checked_results.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
if __name__=='__main__':main()
