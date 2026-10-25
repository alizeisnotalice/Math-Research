"""Saved-results-only final exact crosscheck. No full evaluator execution."""
from pathlib import Path
from fractions import Fraction as Q
import json,hashlib,sys
from flint import arb,ctx,fmpq
sys.set_int_max_str_digits(0);ctx.prec=1024
HERE=Path(__file__).resolve().parent
native=json.loads((HERE/'native_result.json').read_text());ref=json.loads((HERE/'gradient_exact_reference.json').read_text());nv=json.loads((HERE/'native_validation.json').read_text())

def fraction_bounds(row):
    e=int(row['common_binary_exponent']);return[Q(int(row[k]))*Q(2)**e for k in ['lower_mantissa','upper_mantissa']]
def ball(v):return arb(fmpq(v.numerator,v.denominator))
def enclosure(row):
    l,u=fraction_bounds(row);return arb(ball((l+u)/2),ball((u-l)/2))

exact=arb(fmpq(int(ref['gradient_Q_over_pi5_numerator']),int(ref['gradient_Q_over_pi5_denominator'])))/arb.pi()**5
lo,hi=fraction_bounds(native['gradient'][1]);assert ball(lo)<exact.lower()and ball(hi)>exact.upper(),'STOP: independent exact x2 gradient disagreement'
x=[Q(v,ref['denominator'])for v in ref['point_numerators']]
euler=5*enclosure(native['G'])+sum((ball(xi)*enclosure(g)for xi,g in zip(x,native['gradient'])),arb(0))
assert euler.contains(0)
assert ref['full_group_terms']==native['full_group_terms']==967680
out={'status':'PASS_GENERIC_G_AND_ALL_GRADIENT_CHAIN; EXACT_X2_AND_EULER_VALIDATED','full_evaluator_executions':2,
    'unique_full_group_points':1,'new_points':0,'native_G_contains_old_exact_Q_over_pi5':nv['native_G_contains_exact_old_Q_over_pi5'],
    'native_gradient_x2_contains_independent_fmpq_DD6_exact_Q_over_pi5':True,'Euler_saved_full_intervals_contains_zero':True,
    'native_external_seconds':nv['native_external_seconds'],'native_kernel_seconds':nv['native_kernel_seconds'],
    'independent_single_component_fmpq_reference_seconds':ref['elapsed_seconds'],
    'reference_one_component_over_native_G_plus6gradient_time':ref['elapsed_seconds']/nv['native_external_seconds'],
    'old_G_value_only_over_native_G_plus6gradient_time':nv['old_value_only_over_native_value_plus6gradient_time'],
    'timing_comparison_warning':'Recorded workloads differ; not an equal-work or uniform throughput guarantee.',
    'G_display_only':enclosure(native['G']).str(30),'gradient_display_only':[enclosure(g).str(30)for g in native['gradient']],
    'native_source_sha256':hashlib.sha256((HERE/'generic_green.cpp').read_bytes()).hexdigest(),
    'exact_reference_source_sha256':hashlib.sha256((HERE/'exact_reference.py').read_bytes()).hexdigest(),
    'gradient_exact_result_sha256':hashlib.sha256((HERE/'gradient_exact_reference.json').read_bytes()).hexdigest(),
    'stopped':True,'live_handles':[]}
(HERE/'validation_results.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
