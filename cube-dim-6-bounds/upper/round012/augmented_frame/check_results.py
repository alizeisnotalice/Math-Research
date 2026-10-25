#!/usr/bin/env python3
"""Independent saved-data reconstruction; no kernel calls or new points."""
import sys
sys.dont_write_bytecode=True
from pathlib import Path
from fractions import Fraction as Q
from time import perf_counter
import hashlib,json,os
from flint import arb,ctx,fmpq,fmpz
ctx.prec=768
HERE=Path(__file__).resolve().parent
B010=Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube6-upper-astra/cube-dim-6-bounds/upper/round010/short_minimum_diagnostic')
T=Q(2,45)
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def A(q):return arb(q.numerator)/q.denominator
def Gvalue(row):return arb(fmpq(fmpz(row['G_times_pi5_numerator']),fmpz(row['G_times_pi5_denominator'])))/arb.pi()**5
def display(a):return a.str(110,more=True)

start=perf_counter()
result=json.loads((HERE/'results.json').read_text())
assert result['status']=='COMPLETE_TWO_POINTS_NECESSARY_CONDITIONS_ONLY'and result['exit_code']==0
assert result['new_full_passes']==4 and result['new_G_passes']==result['new_H_passes']==2
hash_count=0
for name in ['input_manifest.json','kernel_manifest.json']:
    for item in json.loads((HERE/name).read_text())['inputs']:
        assert sha(item['path'])==item['sha256'],item['path'];hash_count+=1
assert result['kernel_manifest_sha256']==sha(HERE/'kernel_manifest.json')
assert result['source_sha256']==sha(HERE/'evaluate.py')
old1=json.loads((B010/'point1_results.json').read_text())
old32=json.loads((B010/'other_time_results.json').read_text())
old2=json.loads((B010/'batch2_results.json').read_text())
oldGs=[Gvalue(old1['G']),Gvalue(old2['G'])]
oldKs=[oldGs[0]-arb(old32['H_arb384']),oldGs[1]-arb(old2['H_arb384'])]
rows=[]
for idx,row in enumerate(result['points']):
    x=[Q(n,row['denominator'])for n in row['point_numerators']]
    assert x[0]==1 and all(x[j]>x[j+1]>0 for j in range(5))
    oldrow=old1 if idx==0 else old2
    assert row['denominator']==oldrow['denominator'] and row['point_numerators']==oldrow['point_numerators']
    assert row['new_full_passes']==2
    g=row['G'];assert g['complete_group_terms']==967680
    assert g['summed_divergent_Laurent_coefficients']==['0','0','0']
    G=Gvalue(g);assert (G-arb(g['G_arb512']['ball'])).contains(0)
    h=row['H'];assert h['complete_group_terms']==967680
    if idx==0:
        raw=json.loads((HERE/h['native_file']).read_text())
        assert sha(HERE/h['native_file'])==h['native_sha256']
        assert raw['full_group_terms']==967680 and raw['precision_bits']==512
        assert raw['zero_node_terms']==g['zero_node_terms']
        assert raw['confluent_terms']==g['repeated_node_terms']
        dy=raw['dyadic_interval'];factor=Q(2)**int(dy['common_binary_exponent'])
        lo=Q(int(dy['lower_mantissa']))*factor;hi=Q(int(dy['upper_mantissa']))*factor
        assert lo<=hi
        Hlo=A(lo);Hhi=A(hi)
    else:
        assert h['native_engine_called']is False and h['actual_precision_bits']==512
        original=arb(h['original_serialized_H']);Hlo=original.lower();Hhi=original.upper()
    Klo=G-Hhi;Khi=G-Hlo
    assert Klo>A(T/4000)>A(T/4096)
    stored=arb(row['K_beta_T']['ball'])
    assert (stored-Klo).contains(0)and(stored-Khi).contains(0)
    costlo=A(T)/(64*Khi);costhi=A(T)/(64*Klo)
    assert costhi<arb(125)/2
    assert row['necessary_condition_for_below_64_passes']
    assert not row['rules_out_below_64_candidate']and not row['rules_out_62_5_sufficient_target']
    # Coarse exact-rational bounds provide short, auditable displayed intervals.
    coarse=[(Q(1119361581454720,10**20),Q(1119361581454721,10**20)),
            (Q(1115606229724259,10**20),Q(1115606229724260,10**20))][idx]
    assert Klo>A(coarse[0])and Khi<A(coarse[1])
    rows.append({'point_index':idx+1,'K_lower_reconstructed':display(Klo),'K_upper_reconstructed':display(Khi),
                 'coarse_exact_rational_K_enclosure':list(map(str,coarse)),
                 'K_minus_1_over_90000_lower':display(Klo-arb(1)/90000),
                 'K_minus_1_over_92160_lower':display(Klo-arb(1)/92160),
                 'formal_point_cost_lower':display(costlo),'formal_point_cost_upper':display(costhi),
                 'G_zero_node_terms':g['zero_node_terms'],'G_repeated_node_terms':g['repeated_node_terms'],
                 'principal_parts_exactly_zero':True,
                 'new_G':display(G),'old_body32_G_same_exact_point':display(oldGs[idx]),
                 'new_G_minus_old_G':display(G-oldGs[idx]),
                 'old_body32_T_4_over_125_K':display(oldKs[idx]),
                 'old_body32_T_4_over_125_formal_cost':display(arb(1)/(2000*oldKs[idx])),
                 'comparison_changes_frame_and_time':True})
out={'status':'PASS_SAVED_CERTIFICATES_TWO_POINTS_ONLY','pid':os.getpid(),'precision_bits':768,
     'rows':rows,'hash_checks':hash_count,'new_kernel_calls':0,'new_spatial_points':0,
     'cumulative_augmented_full_passes':4,'global_m_certified':False,'C6_established_upper_remains':64,
     'actual_run_status':result['status'],'actual_kernel_batch_seconds':result['elapsed_seconds'],
     'source_sha256':sha(__file__),'results_sha256':sha(HERE/'results.json'),
     'elapsed_seconds':perf_counter()-start}
(HERE/'checked_results.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
