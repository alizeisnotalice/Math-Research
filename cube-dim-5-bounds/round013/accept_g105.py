"""Primary finite evidence binding after the full analytic and code review.

No H, log, or historical baseline is evaluated. The global extension is proved
in primary_globalg105_review.md, rather than inferred from audit status labels.
"""
from pathlib import Path
from fractions import Fraction as F
import hashlib
import json
import os
import time

HERE = Path(__file__).resolve().parent
P = HERE.parent


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    started = time.monotonic()
    regpath = HERE/'primary_g105_registration.json'
    reg = json.loads(regpath.read_text())
    for name, expected in reg['locks'].items():
        assert sha(P/name) == expected, name
    audit = json.loads((P/'audit_mixed100/closed_gate_result.json').read_text())
    assert audit['status'] == 'mixed_H64_H128_complete_closed_cover_numeric_and_structural_gates_passed_for_main'
    assert audit['total_passed'] == 100998 and audit['retained_H64_passed'] == 100804
    assert audit['new_H128_passed'] == audit['fresh_H128_Delta_checks'] == 194
    assert audit['fresh_H128_closed_corners'] == 3104 and audit['penalty'] == 1024
    assert F(audit['minimum_new_corner_margin']) > 0
    assert audit['journal_events'] == 583 and audit['overlay_exact'] and audit['final_every_field_exact']
    assert all(audit[k] == 0 for k in ('queue','inflight','faults','H128_obstructions'))
    assert audit['cache_H64_records_retained'] == 184796
    assert audit['cache_H128_records'] == 194 and audit['typed_cache_total'] == 184990
    assert audit['final_state_sha256'] == reg['full098_state_sha256']
    summary = json.loads((P/'upper_portable096/run098_v3/checker_summary.json').read_text())
    assert summary['status'] == 'PASS_portable_compact_closed_root_packet'
    assert (summary['roots'], summary['leaves'], summary['splits']) == (52360,100998,48638)
    assert summary['actual_unique_epsilon_point_pairs'] == summary['actual_unique_leaf_nodes'] == 100998
    assert summary['epsilon_counts'] == {'1/64':100804,'1/128':194}
    assert summary['quantization_bits'] == {'160':100998,'192':0,'256':0}
    assert summary['g'] == '21/200000'
    assert all(F(x) > 0 for x in summary['minimum_margin_by_mode'].values())
    packet = P/'upper_portable096/run098_v3'
    manifest = json.loads((packet/'manifest.json').read_text())
    assert summary['packet_manifest_sha256'] == sha(packet/'manifest.json')
    assert manifest['source_locks']['state'] == reg['full098_state_sha256']
    assert manifest['source_locks']['source_manifest'] == reg['production098_manifest_sha256']
    assert len(manifest['chunks']) == 52
    total = 0
    for row in manifest['chunks']:
        data = (packet/row['file']).read_bytes()
        assert len(data) == row['bytes'] and hashlib.sha256(data).hexdigest() == row['sha256']
        total += len(data)
    assert total == 13764568
    g, old = F(21,200000), F(1,10000)
    assert g-old == F(1,200000)
    out = {'status':'primary_accepts_complete_same_model_global_g105',
           'qualified_g':str(g), 'previous_global_g':str(old),
           'exact_g_improvement':str(g-old),
           'closed_ordered_face_roots':52360,'closed_leaves':100998,
           'leaf_epsilons':{'1/64':100804,'1/128':194},
           'new_H128_independent_components':970,'new_closed_corners':3104,
           'portable_quantized_intervals_outward_bits':160,
           'primary_analytic_review':'primary_globalg105_review.md',
           'same_epsilon_H_homogeneity_used':False,
           'all_singular_regions_and_closed_boundaries_included':True,
           'same029_general_L1_and_all_positive_common_scales':True,
           'model':'n5 alpha1 body16 unnormalized weight1 common time Fourier (2pi)^-5',
           'full098_state_sha256':reg['full098_state_sha256'],
           'new_C5_bound':False,'formal_C5_upper':'32',
           'dual_breakthrough_time':None,'post_dual_rounds':0,
           'new_H_or_log_calls':0,'pid':os.getpid(),
           'seconds':time.monotonic()-started,
           'registration_sha256':sha(regpath),'source_sha256':sha(__file__),
           'input_locks':reg['locks']}
    target = HERE/'primary_globalg105_acceptance.json'
    assert not target.exists()
    target.write_text(json.dumps(out,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps({k:v for k,v in out.items() if k!='input_locks'},ensure_ascii=False))


if __name__ == '__main__':
    main()
