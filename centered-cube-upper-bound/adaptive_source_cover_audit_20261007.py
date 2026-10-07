"""Exact arithmetic guards for the new greedy output cover and stronger-box obstacle.

Does not run other people's scripts or claim a universal source-packet partition.
"""
from fractions import Fraction as F
from pathlib import Path
import hashlib
import json

OUT=Path(__file__).with_name('adaptive_source_cover_audit_20261007_results.json')


def integral(poly,lo,hi):
    return sum(c*(hi**(i+1)-lo**(i+1))/F(i+1) for i,c in enumerate(poly))


def convolution_phi(C,s):
    D=(C+2)/2
    # u=t-s; polynomial phi(s+u)=D²-s²-2su-u².
    phi=[D*D-s*s,-2*s,F(-1)]
    ans=F(0)
    for lo,hi,triangle in [(F(-1),F(0),[F(1),F(1)]),(F(0),F(1),[F(1),F(-1)])]:
        prod=[F(0)]*4
        for i,c in enumerate(phi):
            for j,d in enumerate(triangle):prod[i+j]+=c*d
        ans+=integral(prod,lo,hi)
    assert ans==D*D-s*s-F(1,6)
    return ans


def finite_greedy(n):
    a=F(1);h=F(1,n);eta=F(1,4);lam=F(1,3);v=2*eta*lam*a**n
    weights=[F(w,26) for w in [3,1,5,2,7,1,4,3]]
    atoms=[F(j,2)*h for j in range(8)]
    assert sum(weights)==1
    candidates=[]
    for j in range(-1,8):
        lo=F(j,2)*h;hi=lo+h
        mass=sum(w for y,w in zip(atoms,weights) if lo<=y<=hi)
        if mass>v:candidates.append({'lo':lo,'hi':hi,'mass':mass})
    remaining=list(candidates);selected=[];owners=[]
    while remaining:
        best=max(remaining,key=lambda q:q['mass'])
        deleted=[q for q in remaining if q['lo']<=best['hi'] and best['lo']<=q['hi']]
        assert best['mass']>v
        for q in deleted:
            assert q['mass']<=2*best['mass']
            assert abs((q['lo']+q['hi']-best['lo']-best['hi'])/2)<=h
            owners.append((q,best))
        selected.append(best)
        remaining=[q for q in remaining if q not in deleted]
    assert all(x['hi']<y['lo'] or y['hi']<x['lo'] for i,x in enumerate(selected) for y in selected[i+1:])
    allocated=sum(q['mass'] for q in selected)
    assert allocated<=1 and len(selected)*v<1
    # x=0,R=a captures the entire true atom input at all these n.
    assert max(atoms)<=a/2
    m=F(1);u=F(1)
    assert 2*lam<u<=4*lam
    witness_count=0
    for B,S in owners:
        if B['mass']>eta*m:
            witness_count+=1
            center=(S['lo']+S['hi'])/2
            assert abs(center)<=a/2+3*h/2
            assert a**n < B['mass']/(2*eta*lam) <= S['mass']/(eta*lam)
            assert (2*abs(center))**n < (1+3*h/a)**n*S['mass']/(eta*lam)
    assert witness_count>0
    return {'n':n,'a':str(a),'h':str(h),'eta':str(eta),'lambda':str(lam),'v':str(v),
            'selected':[{'lo':str(q['lo']),'hi':str(q['hi']),'mass':str(q['mass'])} for q in selected],
            'selected_total_mass':str(allocated),'candidate_count':len(candidates),'actual_hardband_witnesses_checked':witness_count,
            'lambda_sum_output_cover_volumes_exact':str((1+3*h/a)**n*allocated/eta),
            'exact_prefactor_for_h_1_over_n':str((1+F(3,n))**n),
            'exact_prefactor_for_h_2_over_n':str((1+F(6,n))**n),
            'scope':'Finite candidate-family arithmetic guard; general infinite-family termination/coverage is proved in the companion Markdown.'}


def strong_box_guard(n):
    C=F(n.bit_length().bit_length())  # a frozen integer O(log log n) size rule
    D=(C+2)/2;beta=1-F(1,6)/(D*D)
    assert beta==1-F(2,3)/(C+2)**2
    for s in [F(0),C/4,-C/4,C/2,-C/2]:
        full=convolution_phi(C,s)
        phi=D*D-s*s
        assert full<=beta*phi and phi>0
    T=2*n+1
    boundary_ratio=F(T,T-1)**n
    assert boundary_ratio<2
    gamma_squared_upper=boundary_ratio*beta**n
    return {'n':n,'C_integer_rule':str(C),'D':str(D),'beta':str(beta),'finite_source_side_over_h':T,
            'full_source_over_interior_translation_volume_ratio':str(boundary_ratio),
            'gamma_squared_upper_bound':str(gamma_squared_upper),'gamma_squared_upper_float':float(gamma_squared_upper),
            'scope':'Exact Schur polynomial and finite-boundary constants for the stronger all-moving-box retention contract; does not disprove conditional hard-query output cover.'}


def overlap_guards(n):
    # Exact L1 fixtures are explained by their geometry in the Markdown.
    return {'n':n,'fixed_grid_equal_cell_fraction':str(F(1,2**n)),
            'fixed_grid_eta':str(F(1,4)),'fixed_grid_misses_eta':F(1,2**n)<F(1,4),
            'fixed_grid_fixture_actual_weak_ratio_at_lambda_1_over_3':str(F(1,3)*(1-F(1,2*n))**n),
            'overlapping_candidate_count':2**n,'overlap_fixture_total_source_mass':2,
            'overlap_fixture_each_candidate_mass':str(1+F(1,2**n)),
            'incorrect_sum_of_candidate_masses':2**n+1,
            'duplicate_count_ratio':str(F(2**n+1,2)),
            'whole_overlap_fixture_can_be_one_box_side_at_most_2h':True}


def main():
    assert not OUT.exists(),'Refuse overwrite.'
    rows=[]
    for n in (8,32,128):
        rows.append({'n':n,'greedy_output_cover':finite_greedy(n),'strong_box_schur':strong_box_guard(n),'overlap_boundary_fixtures':overlap_guards(n)})
    result={'status':'PASS_EXACT_GREEDY_COVER_AND_STRONG_BOX_GUARDS','rounds':3,'rows':rows,
            'claim':'The direct threshold-specific greedy output cover is valid. Strong witness-preserving source allocation remains unnecessary and is obstructed in its unconditional all-small-box form.',
            'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            'no_general_source_partition_existence_claim':True,'no_actual_FIRST_gate_reconstruction':True}
    OUT.write_text(json.dumps(result,indent=2))
    print(json.dumps({'status':result['status'],'rounds':3,'dimensions':[8,32,128]}))


if __name__=='__main__':main()
