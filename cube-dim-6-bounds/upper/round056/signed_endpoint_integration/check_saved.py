#!/usr/bin/env python3
"""Pure rational saved-data audit of signed endpoint ratios and norm powers."""
import sys
sys.dont_write_bytecode=True
sys.set_int_max_str_digits(0)
from pathlib import Path
from fractions import Fraction as Q
from math import comb
import json,hashlib,time
H=Path(__file__).resolve().parent;OLD=H.parents[1]/'round054/lower_endpoint_green_budget'
def rd(p):return json.loads(Path(p).read_text())
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def b(r):return Q(int(r['lower_mantissa']))*Q(2)**int(r['lower_binary_exponent']),Q(int(r['upper_mantissa']))*Q(2)**int(r['upper_binary_exponent'])
def main():
    tic=time.monotonic();reg=rd(H/'continuation_registration.json');legacy=rd(H/'registration.json');release=rd(H/'gate_release.json');r=rd(H/'results.json');receipt=rd(H/'run_receipt.json');old=rd(OLD/'results.json')
    assert release['final_F054_PASS']is True
    assert release['parent_final_release']is True
    assert not(H/'checked_results.json').exists(),'one pure saved-data audit only'
    assert reg['predecessor']['registration_sha256']==sha(H/'registration.json')
    assert reg['fixed_plan']==legacy['fixed_plan'] and reg['max_configurations']==10
    assert release['continuation_registration_sha256']==sha(H/'continuation_registration.json')
    assert receipt['scientific_budget_seconds']==reg['scientific_execution']['max_wall_seconds']==1200
    assert receipt['scientific_deadline_epoch']==receipt['started_epoch']+1200
    assert time.time()<receipt['scientific_deadline_epoch']
    assert receipt['status']=='COMPLETED' and receipt['attempted_configurations']==receipt['completed_configurations']==9 and receipt['exit_code']==0
    assert receipt['registration_sha256']==sha(H/'continuation_registration.json') and receipt['predecessor_registration_sha256']==sha(H/'registration.json') and receipt['gate_release_sha256']==sha(H/'gate_release.json') and receipt['results_sha256']==sha(H/'results.json')
    for p,x in reg['source_sha256'].items():assert sha(p)==x,p
    for p,x in release['audit_source_sha256'].items():assert sha(p)==x,p
    for p,x in release['production_code_sha256'].items():assert sha(H/p)==x,p
    assert set(release['production_code_sha256'])=={'theorem.md','evaluate.py','check_saved.py'}
    assert reg['theorem_sha256']==sha(H/'theorem.md')
    assert release['bound_A054_manifest_sha256']==sha(OLD/'manifest.json')
    q=Q(6,5);aa=Q(84,19);bb=Q(7,2)
    assert (4-q)*q/((q-1)*(5-q))==aa and (4-q)/(4*(q-1))==bb
    assert aa-bb==Q(35,38) and aa-q*bb==Q(21,95)
    for k,row in old['cases'].items():
        assert r['source_054_intervals'][k]=={f:row[f]for f in ['singular_norm_upper','background_norm_upper','total_norm_upper']}
    powers={};ratio_checks=[]
    for h in [1,2]:
        row=r['R'][str(h)];al,au=b(row['alpha']);assert 0<al<=au<1 and al*al<=Q(1,256*h)<=au*au
        RU=b(row['R'])[1]*(1+Q(1,10**18));assert 0<RU<1
        # R <= aa*alphaU - bb*alphaL^(6/5); prove by fifth powers.
        z=(aa*au-RU)/bb
        if z>0:assert z**5<=al**6
        PU=b(row['R_power_5_over6'])[1]*(1+Q(1,10**16));assert PU**6>=RU**5 and 0<PU<1
        powers[h]=PU;ratio_checks.append(dict(h=h,alpha_square_check=True,R_upper_by_fifth_power=True,R_power_upper_by_sixth_power=True))
    weights=[1,5,10,10,5];cert={};rows=[]
    for k in range(1,6):
        m=6-k;row=r['cases'][str(k)];prior=old['cases'][str(k)];N={}
        assert row['n']==comb(5,k-1)==weights[k-1]
        for h in range(1,m//2+1):
            n=comb(m,h)
            if h*2==m:assert n%2==0;n//=2
            N[str(h)]=n
        assert 1+sum(N.values())==2**(m-1) and prior['incident_ordered_pairs']==2**m
        assert row['background_norm_upper']==prior['background_norm_upper']
        if k==5:
            for f in ['singular_norm_upper','background_norm_upper','total_norm_upper']:assert row[f]==prior[f]
            cert[k]=b(prior['total_norm_upper'])[1]
        else:
            assert row['N_mh']==N and not row['background_recomputed'] and row['all_incident_pairs_and_background_paid']
            FU=(1+sum((n*powers[int(h)]for h,n in N.items()),Q(0)))/2**(m-1);assert 0<FU<1
            cert[k]=FU*b(prior['singular_norm_upper'])[1]+b(prior['background_norm_upper'])[1]
            assert cert[k]<Q(row['simple_rational_upper'])<b(prior['total_norm_upper'])[0]
        rows.append(dict(k=k,m=m,N_mh=N,planes=2**(m-1),background_identical=True,k5_no_reevaluation=(k==5)))
    total=Q(0)
    for t in r['union']['weighted_terms']:
        k=t['k'];n=weights[k-1];assert t['n']==n
        up=b(t['weighted_power'])[1]*(1+Q(1,10**14));assert (up/n)**5>=cert[k]**6;total+=up
    final=Q(r['union']['simple_rational_union_upper']);assert final**6>total**5
    assert r['union']['source_old_union_interval']==old['union']['full_time_union_norm_upper']
    assert final<b(old['union']['full_time_union_norm_upper'])[0]
    assert r['union']['overlap_permitted']and not r['union']['domains_claimed_disjoint'] and r['union']['weights']==weights
    finished=time.time();assert finished<receipt['scientific_deadline_epoch']
    out=dict(status='PASS_RATIONAL_SAVED_DATA_AUDIT',input_results_sha256=sha(H/'results.json'),registration_sha256=sha(H/'continuation_registration.json'),predecessor_registration_sha256=sha(H/'registration.json'),gate_release_sha256=sha(H/'gate_release.json'),ratio_checks=ratio_checks,type_checks=rows,domain_weights=weights,domain_count=31,source_M_recomputations=0,old_complete_budget_recomputations=0,k5_new_evaluations=0,new_Arb_or_kernel_calls=0,producer_configurations=9,saved_data_audit_configurations=1,total_configurations=10,simple_rational_union_upper=str(final),union_power_check='final^6 > sum(certified weighted powers)^5',strictly_below_old_formula_lower_endpoint=True,source_hashes_unchanged=True,new_C6_bound=False,elapsed_seconds=time.monotonic()-tic,finished_epoch=finished,scientific_total_wall_seconds=finished-receipt['started_epoch'],scientific_deadline_epoch=receipt['scientific_deadline_epoch'])
    (H/'checked_results.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
if __name__=='__main__':main()
