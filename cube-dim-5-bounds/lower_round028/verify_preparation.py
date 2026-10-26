#!/usr/bin/env python3
"""Audit preparation only; no actual epsilon functional evaluations."""
from pathlib import Path
from fractions import Fraction as Q
from itertools import combinations
from collections import defaultdict
import sys,os,time,json,gzip
sys.dont_write_bytecode=True
P=Path(__file__).resolve().parent;ROOT=P.parent;OLD=ROOT/'lower_round025'
sys.path[:0]=[str(ROOT/'lower'),str(P)]
from certify_rectangular import field
from full_static import ENV,load,save,sha,rt,decode,encode,mul,addto,ZERO
from verify_pairs import stream_histogram
def main():
    start=time.monotonic();reg=load(P/'preparation_verification_registration.json');assert sha(Path(__file__))==reg['source_sha256']
    for x in reg['inputs']:assert sha(Path(x['path']))==x['sha256']
    F,u,v,bounds,_,_=field();zero=F(0);a=1+u**3/2
    def dec(c):return F(a=decode(c))
    for oldname in ('lower_round025','lower_round026'):
        old=ROOT/oldname
        for x in load(old/'archive_whitelist.json')['files']:assert sha(old/x['path'])==x['sha256']
    geom=load(OLD/'tail_geometry.json');fams={x['family_id']:x for x in geom['families']}
    for fid,row in fams.items():
        r=dec(row['radius']);ell=row['lengths'][0]
        assert 32*r**5==F(row['N']) and r<=F(2) and r<ell*a
        if ell==4:assert 2*r<=F(4)<3*a
    formsfile=load(P/'capture_forms.json');assert len(formsfile['all6144_raw_inventory'])==6144
    forms=[(dec(x['constant']),x['slope']) for x in formsfile['forms']];events=load(P/'events.json')['events'];eventvalues=[dec(x['value']) for x in events];assert len(forms)==279 and len(events)==965
    assert eventvalues[0]==zero and eventvalues[-1]==a
    assert all(x<y for x,y in zip(eventvalues,eventvalues[1:]))
    eventset={x.a for x in eventvalues};seen={zero.a,a.a};pairs=0
    for (c,m),(d,n) in combinations(forms,2):
        pairs+=1
        if m==n:continue
        eps=(d-c)/(m-n)
        if zero<=eps<=a:assert eps.a in eventset;seen.add(eps.a)
    assert seen==eventset and pairs==38781
    # Independently stream the archived static histogram to reproduce all91
    # descriptor weights; no SQL primary-state dependency and no geometrical replay.
    desc=defaultdict(lambda:[ZERO]*25);maskrows=0
    for maskstr,coeff in stream_histogram(P/'full_tail_histogram.json.gz'):
        bits=int(maskstr);best=[0,0,0]
        while bits:
            bit=bits&-bits;fid=bit.bit_length()-1;bits-=bit
            row=fams[fid];ell=row['lengths'][0]
            if ell<4:best[ell-1]=max(best[ell-1],row['N'])
        addto(desc[tuple(best)],decode(coeff));maskrows+=1
    with gzip.open(P/'descriptors.json.gz','rt') as f:saved=json.load(f)['descriptors']
    assert len(saved)==len(desc)==91 and maskrows==139843
    for row in saved:assert decode(row['physical_tail_weight'])==desc[tuple(row['max_N_by_ell1'])]
    total=[ZERO]*25
    for value in desc.values():addto(total,value)
    assert F(a=total)==((1+u)/4)**4
    upperrows=load(P/'relaxation_checkpoint.json')['rows'];assert len(upperrows)==91
    numerator=[ZERO]*25
    for row in upperrows:
        key=tuple(row['max_N_by_ell1']);length=dec(row['swept_interval_union_length_upper']);assert zero<=length<=3*a
        addto(numerator,mul(desc[key],length.a))
    ratio=F(a=numerator)*Q(16,3);rr=load(P/'relaxation_result.json');assert ratio==dec(rr['upper_ratio_coefficients'])
    C=load(ROOT/'lower_round003/a2_equals_a3_results.json');Cexact=dec(C['coefficients']);residual=ratio-Cexact
    assert residual==dec(rr['upper_minus_complete_C_coefficients']) and bounds(residual.a)[0]>0
    # No new epsilon union: inherit oldC and its accurate right derivative, then
    # extend that affine identity to the first complete global topology cell.
    D=dec(load(OLD/'plus_derivative.json')['directional_derivative_coefficients']);first=eventvalues[1]
    assert zero<first and bounds(D.a)[1]<0
    contact_difference=D*first;assert bounds(contact_difference.a)[1]<0
    firstcert={'status':'INHERITED_exact_first_open_cell_exclusion','left_epsilon':'0','right_event_index':1,'right_epsilon_coefficients':encode(first.a),'right_epsilon_interval':list(map(str,bounds(first.a))),'ratio_expression':'R(epsilon)=C+epsilon*D025 on0<=epsilon<=first_positive_global_event','D_coefficients':encode(D.a),'right_contact_minus_C_coefficients':encode(contact_difference.a),'right_contact_minus_C_interval':list(map(str,bounds(contact_difference.a))),'strictly_below_C_for_0_less_epsilon_le_first_contact':True,'no_new_baseline_or_epsilon_union':True,'justification':'completeglobalendpointarrangementgivesPWA;025suppliesCidentityandD;continuityextendsfirstcelltorightcontact','source_sha256':reg['source_sha256']}
    save(P/'first_cell_certificate.json',firstcert)
    pair=load(P/'pair_verification_result.json');assert pair['status']=='PASS_entire_full213_histogram_independent_pair_factorization'
    out={'status':'PASS_complete_nonlocal_preparation_not_window_exclusion','event_pairs_rechecked':pairs,'complete_closed_events':965,'artifact_masks_streamed':maskrows,'descriptor_weights_reproduced':91,'descriptor_total_physical_mass_exact':True,'safe_upper_exactly_above_C':True,'first_cell_inherited_exclusion_only':True,'full213_pairfactorization_all_coefficients_equal':139843,'old025026_archive_files_unchanged':True,'actual_epsilon_functional_evaluations':0,'no_unique_contact_selected':True,'runtime':rt(start),'source_sha256':reg['source_sha256']}
    save(P/'verification.json',out);print(out,flush=True)
if __name__=='__main__':
    assert all(os.environ.get(k)=='1' for k in ENV) and os.getpriority(os.PRIO_PROCESS,0)>=10
    main()
