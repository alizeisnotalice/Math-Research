#!/usr/bin/env python3
"""Pure saved-data scope and exact rational-sign audit; no new geometric point."""
import sys
sys.dont_write_bytecode=True
from pathlib import Path
from fractions import Fraction as Q
import json,hashlib,time,os
H=Path(__file__).resolve().parent
def rd(p):return json.loads(Path(p).read_text())
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def main():
    tic=time.monotonic();assert not(H/'checked_results.json').exists()
    reg=rd(H/'registration.json');bind=rd(H/'execution_binding.json');r=rd(H/'results.json');p=rd(H/'run_receipt.json');pg=rd(H/'producer_process_receipt.json');cg=rd(H/'checker_process_receipt.json')
    assert cg['status']=='RUNNING'and cg['owned_pid']==os.getpid()and cg['work_timeout_seconds']==170
    assert pg['status']=='COMPLETED'and pg['exit_code']==0 and pg['elapsed_seconds']<180 and pg['within180_including_cleanup']and pg['child_reaped']
    assert p['status']=='COMPLETED'and p['completed_configurations']==p['attempted_configurations']==8
    assert p['registration_sha256']==bind['registration_sha256']==sha(H/'registration.json')and p['binding_sha256']==sha(H/'execution_binding.json')and p['results_sha256']==sha(H/'results.json')
    for f,x in reg['source_sha256'].items():assert sha(f)==x,f
    for f,x in bind['code_sha256'].items():assert sha(H/f)==x,f
    assert [c['id']for c in r['configurations']]==list(range(1,9))
    cs={c['id']:c for c in r['configurations']}
    assert cs[2]['point']==[1,2,4,0,0,0]and cs[6]['point']==[1,1,1,0,0,0]
    for ident,count in [(2,15),(6,28)]:
        rows=cs[ident]['pairs'];assert len(rows)==count and len({tuple(z['pair'])for z in rows})==count
        vals=[Q(z['residual_squared'])for z in rows];assert all(x>0 for x in vals)and min(vals)==Q(cs[ident]['min_pair_residual_squared'])
    assert Q(cs[4]['value_at_zero_times_pi5'])==Q(7,1440)
    assert Q(cs[4]['right_derivative_times_pi5'])==-Q(1,128)and Q(cs[4]['left_derivative_times_pi5'])==Q(1,128)
    assert Q(cs[4]['derivative_jump_times_pi5'])==-Q(1,64)and Q(cs[4]['right_second_derivative_limit_times_pi5'])>0
    assert cs[6]['W1_rank']+cs[6]['W2_rank']-cs[6]['union_rank']==cs[6]['intersection_dimension']==2
    assert cs[6]['Green_evaluations']==cs[6]['Hessian_evaluations']==0 and not cs[6]['semiconcavity_claim_for_this_frame']
    assert cs[7]['second_derivative_limit']=='+infinity'and Q(cs[7]['derivative_jump'])==-2 and not cs[7]['Cauchy_frame_counterexample']and not cs[7]['body32_counterexample']
    assert r['actual_work']['exact_tiny_pair_residuals']==43 and r['actual_work']['true_tiny_Green_symbolic_paths']==1
    assert not r['body32_semiconcavity_proved']and not r['body32_counterexample_found']and not r['new_Cn_bound']
    result=dict(status='PASS_SAVED_SCOPED_RESULTS_AND_BOUNDARIES',results_sha256=sha(H/'results.json'),registration_sha256=sha(H/'registration.json'),binding_sha256=sha(H/'execution_binding.json'),producer_configurations=8,saved_data_audit_configurations=1,total_configurations=9,saved_pair_residual_signs_checked=43,new_pair_residual_evaluations=0,new_Green_evaluations=0,axis_negative_cusp_exact=True,axis_one_sided_second_limit_finite=True,overlap_fixture_only_geometry=True,abstract_counterexample_not_frame_counterexample=True,axis_independent6_author_theorem=True,body32_semiconcavity_proved=False,body32_counterexample_found=False,new_Cn_bound=False,source_hashes_unchanged=True,elapsed_seconds=time.monotonic()-tic)
    (H/'checked_results.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
if __name__=='__main__':main()
