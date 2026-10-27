#!/usr/bin/env python3
"""Freeze completed metadata receipts only. No endpoint/math audit is called."""
from pathlib import Path
import hashlib,json,os,shutil,time
H=Path(__file__).resolve().parent
T=Path('/Users/zhengzhihao/Desktop/T')
P=T/'Math-Research-cube5-upper/cube-dim-5-bounds/upper_bandcurvature173'
S=T/'cube-finite-dim-scratch/upper_bandcurvature173/reviewed_v3_synthetic'
A=T/'cube-finite-dim-scratch/audit_bandcurvature174/binding020_run001'
def sha(p):
    h=hashlib.sha256()
    with Path(p).open('rb') as f:
        for b in iter(lambda:f.read(1<<20),b''):h.update(b)
    return h.hexdigest()
def read(p):return json.loads(Path(p).read_text())
def put(p,x):
    with Path(p).open('x') as f:
        json.dump(x,f,sort_keys=True,indent=2);f.write('\n');f.flush();os.fsync(f.fileno())
def ref(p):return {'path':str(p),'sha256':sha(p)}
def copy_receipt(p,q):
    q.parent.mkdir(parents=True,exist_ok=True)
    assert not q.exists()
    shutil.copyfile(p,q);assert sha(p)==sha(q)
    return {'original':ref(p),'copy':q.relative_to(H).as_posix(),'sha256':sha(q),'bytes':q.stat().st_size}
def main():
    began=time.monotonic()
    assert not (H/'manifest.json').exists()
    root=H.parent;old=read(root/'manifest.json')
    assert sha(root/'manifest.json')=='05c00900b867bad231394285f74167687977f09a8e3cb777824fbd07266a1179'
    for p,s in old['files'].items():assert sha(root/p)==s,p
    assert sha(H/'qualified_donor_endpoint_provenance_v2.json')=='87fa22239a24b1c2c818348648c80899b6246e1a3830eaa7512884647495834c'
    versions={};manifests={}
    for name,expected in [('', 'b218200d7ae0d185509ecb39bc0f11abc52a12842219bd885bbc57dd4bf1b3dd'),('reviewed_v2','db85430d700ddefb04c002c899ada1d1aed05547f3777497649ffef793340448'),('reviewed_v3','cceaa29e847719bc90cbbc1ef7005f8cf2c62c1924772d942003740233b773eb')]:
        base=P/name;m=read(base/'prepared_manifest.json');assert sha(base/'prepared_manifest.json')==expected
        counts={}
        for group in ('repository_files','source_locks','prepared_external_files'):
            for p,s in m[group].items():assert sha(base/p if group=='repository_files' else p)==s,(name,group,p)
            counts[group]=len(m[group])
        versions[name or 'original']={'manifest':ref(base/'prepared_manifest.json'),'all_locks_match':True,'lock_counts':counts};manifests[name]=m
    m=manifests['reviewed_v3'];v1=manifests[''];v2=manifests['reviewed_v2']
    assert v1['source_locks'].items()<=m['source_locks'].items() and v1['prepared_external_files']==m['prepared_external_files']
    for p in ('core173.py','worker173.py','util173.py','fault_cleanup173.py'):assert sha(P/'reviewed_v2'/p)==sha(P/'reviewed_v3'/p)
    r=read(S/'test_receipt.json');assert sha(S/'test_receipt.json')=='cbe2996aa5962f13178e30e2cac2a2ac31d09a872302eeea2d072d547d9d6ee7'
    O=S/'fault_report_and_stuck_peer/work/run_v3_001'
    for field,name in [('stage_receipt_sha256','actual_attempt0001_mass_exit_receipt.json'),('cleanup_receipt_sha256','actual_attempt0001_mass_fault_cleanup.json'),('master_summary_sha256','actual_attempt0001_summary.json')]:assert sha(O/name)==r[field]
    e=read(O/'actual_attempt0001_mass_exit_receipt.json');c=read(O/'actual_attempt0001_mass_fault_cleanup.json');s=read(O/'actual_attempt0001_summary.json')
    assert len(e['summaries'])==1 and e['summaries']['0']['status']=='synthetic_reported_fault' and 'enter_cleanup_now' in e['fault']['exception']
    assert c['all_joined'] and len(c['workers'])==15 and all(not p['alive'] and p['actual_join_called'] for p in c['workers'])
    assert next(p['exitcode'] for p in c['workers'] if p['pid']==r['stuck_peer_pid'])==-9
    assert s['status']=='fault_prefix_preserved' and s['candidate'] is None and r['master_actual_exitcode']==2
    assert r['new_real_H_jet_donor_band_q_corner_calls']==0 and not r['old_v2_three_toys_repeated']
    copies=[]
    for p in sorted(A.rglob('*')):
        if p.is_file():copies.append(copy_receipt(p,H/'evidence/source_audit'/p.relative_to(A)))
    for p in (S/'test_receipt.json',O/'actual_attempt0001_mass_pool_registration.json',O/'actual_attempt0001_mass_exit_receipt.json',O/'actual_attempt0001_mass_fault_cleanup.json',O/'actual_attempt0001_summary.json'):
        copies.append(copy_receipt(p,H/'evidence/runtime_v3'/p.name))
    put(H/'receipt_copies.json',{'method':'byte_copy_only_no_source_or_math_replay','files':copies})
    runtime={'status':'PASS_theory_interface_and_v3_runtime_review_only','versions':versions,'v2_runtime_accepted':False,'v2_gap':'fault stop could remain in summaries loop before finally;stuck peer blocks cleanup','v3_change':'two stop checks raise immediately into fault-only cleanup;new output namespace','same_v2_math_worker_util_and_cleanup_bytes':True,'original_fixed_source_and_external_inputs_unchanged':True,'v3_code_hashes':m['repository_files'],'fault_contract':{'normal_wall_cap':None,'after_observed_fault_seconds':[5,2,2],'preserve_prefix':True,'actual_join':True,'fault_exitcode':2,'OS_unresponsive_exitcode':3,'automatic_restart':False},'integrated_synthetic_read_only':{'receipt':ref(S/'test_receipt.json'),'test_pid':r['test_pid'],'master_pid':r['synthetic_master_pid'],'master_actual_exitcode':r['master_actual_exitcode'],'summaries_at_cleanup_entry':len(e['summaries']),'stuck_peer_pid':r['stuck_peer_pid'],'stuck_peer_actual_exitcode':-9,'all15_joined':True,'fake_prefixes_preserved':len(c['prefix_before_cleanup']),'wall_seconds':r['wall_seconds'],'old_or_new_toy_reruns_by_this_auditor':0,'synthetic_is_not_mathematical_qualification':True},'machine_gate_contract':{'primary_authorization_sha_required':True,'exact_reviewed_manifest_sha_required':True,'fixed_counts':[455,35872,15],'all_three_gate_hashes_and_accepted_declarations_required':True,'known168':m['known_qualified168'],'primary_must_choose_exact_gate_hashes_not_just_names':True},'source_gate':ref(H/'qualified_donor_endpoint_provenance_v2.json'),'new_H_band_q_corner_evaluations':0,'protected_original174_files_unchanged':True}
    put(H/'runtime_review.json',runtime)
    theory=T/'Math-Research-cube-dim-5/cube-dim-5-bounds/primary/band_curvature172.md'
    assert sha(theory)=='0b5f37e251d5e019120ab7f60ccb7d47c237f7f71c785962b6a247b0a6cb9131'
    gate={'gate':'qualified172_174','accepted':True,'scope':'fixed172 theory/interface and reviewed173_v3 prepared runtime control only;not new band/q/corner output or H enclosure truth','reviewed_manifest':ref(P/'reviewed_v3/prepared_manifest.json'),'reviewed_manifest_sha256':sha(P/'reviewed_v3/prepared_manifest.json'),'rejected_runtime_manifest_sha256':[sha(P/'prepared_manifest.json'),sha(P/'reviewed_v2/prepared_manifest.json')],'theory_source':ref(theory),'independent174_theory_manifest':ref(root/'manifest.json'),'runtime_review':ref(H/'runtime_review.json'),'qualified_donor_endpoint_provenance':ref(H/'qualified_donor_endpoint_provenance_v2.json'),'fixed_input_registry':ref(P/'packet/source_registry.json.gz'),'fixed_input_registry_uncompressed_sha256':'0e8af6a4f1ca7dcb4f252dceb681547203629477b4821b068471529d6963ed8d','fixed_finite_plan':ref(P/'packet/finite_plan.json'),'code_hashes':m['repository_files'],'known_qualified168':m['known_qualified168'],'source_count':455,'weighted_band_terms':22295,'maximum_pairs':35872,'maximum_workers':15,'positive_K_denominator':'16*max(original_g_H64_lower,donor_H64_lower)','153_negative_term_denominator':'original_g_H64_upper_unchanged','require_nonempty_same_cutoff_intersection':True,'strict_positive_band_upper_diagnostic':True,'unchanged_original_Delta_and_physical_gamma':True,'old71_pair_guard':'evaluate only on strict decrease of same original qualified q','same_anchor_all16_and_complete_domain_min_required':True,'shared_H_endpoint_and_jet_truth_dependency':'existing017/DD/Arb and accepted080/094/102/113/121 and g151/155/157/180 source chain;not a new independent H enclosure;167 only transport','new_H_log_gradient_band_q_corner_evaluations':0,'does_not_authorize_dispatch':True,'requires_separate_primary_actual_dispatch_authorization':True,'primary_retains_final_proof_and_global_bound_responsibility':True}
    put(H/'qualified172_174.json',gate)
    final={'status':'FROZEN_both_prerequisite_gates_accepted_v3_only_wait_primary_dispatch','gate_files':{p:sha(H/p) for p in ('qualified172_174.json','qualified_donor_endpoint_provenance_v2.json')},'source_audit_actual_exit1_metadata_vocabulary_fault_preserved':True,'source_workers_all15_exit0':True,'source_cases_replayed':0,'metadata_only_continuation_actual_exit0':True,'new_H_q_band_corner_log_gradient_evaluations':0,'read_only_byte_lock_checks_for_versions':versions,'old174_manifest_and9_files_unchanged':True,'freeze_transport_only':True,'seconds':time.monotonic()-began,'pid':os.getpid()}
    put(H/'final_check.json',final)
    files={p.relative_to(H).as_posix():sha(p) for p in sorted(H.rglob('*')) if p.is_file() and '__pycache__' not in p.parts}
    put(H/'manifest.json',{'task':'audit_bandcurvature174_binding020','status':final['status'],'files':files,'whitelist_files':len(files),'whitelist_bytes':sum((H/p).stat().st_size for p in files),'gates':final['gate_files'],'source_plan_sha256':sha(H/'plan.json'),'continuation_plan_sha256':sha(H/'continuation_plan.json'),'production_v3_manifest_sha256':gate['reviewed_manifest_sha256'],'old174_theory_manifest_sha256':sha(root/'manifest.json'),'new_real_math_calls':0,'future_math_not_authorized':True,'primary_retains_final_proof_responsibility':True})
    print(json.dumps({'manifest_sha256':sha(H/'manifest.json'),'gates':final['gate_files'],'files':len(files),'bytes':sum((H/p).stat().st_size for p in files),'pid':os.getpid(),'seconds':time.monotonic()-began},sort_keys=True))
if __name__=='__main__':main()
