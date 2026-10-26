"""Primary exact complete-I and qualified-g budget binding, without heat reruns."""
from pathlib import Path
from fractions import Fraction as F
from math import isqrt
import hashlib
import json
import os
import time

HERE = Path(__file__).resolve().parent
P = HERE.parent


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def fourth_enclosure(q):
    D = 1 << 256
    a = isqrt(isqrt(q.numerator*D**4//q.denominator))
    assert a**4*q.denominator <= q.numerator*D**4 < (a+1)**4*q.denominator
    return [str(F(a,D)), str(F(a+1,D))]


def main():
    started = time.monotonic()
    rp = HERE/'primary_I095_registration.json'
    reg = json.loads(rp.read_text())
    for path, expected in reg['locks'].items():
        assert sha(P/path) == expected, path
    a = json.loads((P/'audit_adaptive102/results.json').read_text())
    assert a['status'] == 'pass' and a['new_children'] == a['all280_guards'] == 13406
    assert a['new_endpoint_inputs'] == 978638 and a['exact_H64_reused'] == 32
    assert a['all_new_range_mean_bins_exact'] and a['rational_sum_fallbacks'] == 0
    assert a['new_H_evaluated_by_auditor'] == a['new_H_log_exp_Arb_calls'] == 0
    assert a['full_ledger']['final_nodes'] == 69545
    assert a['full_ledger']['physical_frontier'] == 68201
    assert a['full_ledger']['volume'] == '1' and a['full_ledger']['all_parent_prefixes'] == 1024
    s = json.loads((P/'upper_adaptive095/summary.json').read_text())
    b = s['after']
    U, B, E, I = map(F, (b['UH_truncated_integral'],b['B_upper'],
                          b['source_E076_upper'],b['I_upper']))
    cap, old = F(54702517,10**12), F(6010267,10**11)
    assert F(a['full_ledger']['I_simple']) == cap
    assert F(a['full_ledger']['budget']['UH']) == U
    assert B > E > 0 and (B-E)**5 >= U**4 and I**4 >= B**5
    assert I <= cap < old and b['no_tube_added']
    gproof = json.loads((HERE/'primary_globalg105_acceptance.json').read_text())
    assert gproof['status'] == 'primary_accepts_complete_same_model_global_g105'
    g = F(gproof['qualified_g'])
    costs = {}
    for label, gg in [('g100',F(1,10000)),('g105',g)]:
        q = 3125*cap**4/gg**5
        enclosure = fourth_enclosure(q)
        assert q > 32**4 and F(enclosure[0]) > 32
        costs[label] = {'g':str(gg),'paired_cap_cost_fourth_power':str(q),
                        'cost_enclosure':enclosure,
                        'display_only':float(F(enclosure[1])),
                        'this_budget_proves_cost_below32':False}
    out = {'status':'primary_accepts_complete_same_alpha1_I095_and_g105_pair',
           'I_upper':str(cap),'old092_I_cap':str(old),
           'old092_cap_minus_new_exact':str(old-cap),
           'old088_cap_minus_new_exact':str(F(70301737,10**12)-cap),
           'old010_cap_minus_new_exact':str(F(71555899,10**12)-cap),
           'UH':str(U),'B_upper':str(B),'shorttime_E_norm_upper':str(E),
           'complete_I_formula':'I <= (UH^(4/5)+E_norm)^(5/4)',
           'full_face_and_all_positive_time':True,'extra_tube_added':False,
           'qualified_g':str(g),'paired_cap_costs':costs,
           'new_C5_bound':False,'formal_lower_exact_difference_from_ref':'0',
           'formal_upper_exact_difference_from_ref':'0','formal_C5_upper':'32',
           'public_opening_continuation_comparison':'all three formal upper comparisons remain equality32; original full algebraic lower unchanged',
           'dual_breakthrough_time':None,'post_dual_rounds':0,
           'primary_review':'095 report/mechanism/kernel/run/portable and102 protocol/check/report all read; convex chord, full parent orbit, short-time norm and029 retained',
           'independence_limits':'017 fan/theory/H Arb source shared;102 new full range/mean/chord/geometry/ledger independently rederived;11 PF stencil samples',
           'new_H_or_log_calls':0,'pid':os.getpid(),
           'seconds':time.monotonic()-started,'registration_sha256':sha(rp),
           'source_sha256':sha(__file__),'input_locks':reg['locks']}
    path = HERE/'primary_I095_acceptance.json'
    assert not path.exists()
    path.write_text(json.dumps(out,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps({k:out[k] for k in ('status','I_upper','qualified_g',
        'old092_cap_minus_new_exact','pid','seconds','new_C5_bound')},ensure_ascii=False))
    print(json.dumps({k:v['display_only'] for k,v in costs.items()}))


if __name__ == '__main__':
    main()
