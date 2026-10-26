"""Bind primary fixed-functional review to the unchanged full C certificate."""
from pathlib import Path
from fractions import Fraction as F
import hashlib
import json
import os
import time
H=Path(__file__).resolve().parent
P=H.parent


def main():
    t=time.monotonic()
    names=['lower_round003/a2_equals_a3_results.json',
           'lower_round028/archive_whitelist.json',
           'lower_round028/pwa/archive_whitelist.json',
           'lower_round028/pwa/window_result.json',
           'lower_round028/pwa/exact_cells.py',
           'audit_lower028_106/manifest.json',
           'audit_lower028_106/results.json','audit_lower028_106/report.md']
    locks={x:hashlib.sha256((P/x).read_bytes()).hexdigest() for x in names}
    c=json.loads((P/names[0]).read_text())
    w=json.loads((P/'lower_round028/pwa/window_result.json').read_text())
    a=json.loads((P/'audit_lower028_106/results.json').read_text())
    assert locks[names[0]]=='f385811e5460af63abc117804965461eed88d8ecb5b9d2cd30268387646f45d2'
    assert w['global_maximum_ratio_coefficients']==c['coefficients']
    assert w['zero_residual_event_indices']==a['closed_zero']==[0]
    assert w['closed_events']==a['closed_events']==965
    assert w['open_cells']==a['completed_cells']==964
    assert w['new_descriptor_dual_union_objects']==a['weighted_response_rows']==87633
    assert a['closed_negative']==964 and a['adjacent_exact_continuities']==963
    assert a['static_masks_read']==139843 and a['descriptor_weights']==91
    assert F(a['old_D_upper'])<0
    assert w['strictly_below_C_for_every_0_less_epsilon_less_a']
    assert not a['uniform_positive_epsilon_gap_claim'] and not a['C5_upper_bound_claim']
    out={'status':'primary_accepts_fixed_z_period3_position_finite_functional_exclusion',
         'scope':'z=(1,0,-1),0<epsilon<a1,registered1024length functional only',
         'maximum_on_closed_extension':'C_cert exactly,only epsilon=0',
         'coefficients_compared_to_full_C':25,'positive_epsilon_uniform_gap':False,
         'collision_endpoint':'continuous limit only; not distinct-atom input',
         'primary_read':'028 preparation+PWA report,exact_cells source and106 full mathematical report',
         'independence_limits':'106 independent sparse arithmetic/root isolation and response aggregation; frozen fine/pair geometry and producer dual-union responses retained as shared dependencies',
         'general_L1_all_scales':'originalC proof retained; this exclusion does not upper-bound this input at every scale or arbitrary inputs',
         'new_C5_lower_or_upper':False,'difference_from_continuation_lower':'0',
         'dual_breakthrough_time':None,'post_dual_rounds':0,
         'new_union_or_old_baseline_calls':0,'pid':os.getpid(),'seconds':time.monotonic()-t,
         'input_hashes':locks,'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    p=H/'primary_lower028_review.json';assert not p.exists()
    p.write_text(json.dumps(out,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps({k:out[k] for k in ['status','coefficients_compared_to_full_C','pid','seconds']}))


if __name__=='__main__':main()
