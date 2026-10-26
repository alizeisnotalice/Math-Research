"""UNRUN saved-artifact checker; independent parent release required.

Checks bindings, recorded execution, exact dyadic signs and saved arithmetic
consistency. Does not rerun the generator, physical formulas, logs, or peaks.
The default177-second fuse remains armed through final receipt and exit.
"""
import signal
signal.signal(signal.SIGALRM,signal.SIG_DFL)
signal.setitimer(signal.ITIMER_REAL,177.0)

import sys
sys.dont_write_bytecode=True
import datetime,hashlib,json,os,resource,time,traceback
from pathlib import Path
from fractions import Fraction as Q
resource.setrlimit(resource.RLIMIT_CPU,(60,60))
P=Path(__file__).resolve().parent
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
release_path=P/'checker_release.json'
if not release_path.exists():raise SystemExit('104 saved checker PREPARED_NOT_RELEASED')
release=json.loads(release_path.read_text())
if not(release.get('authorized')is True and release.get('round')=='104-saved-check' and release.get('script_sha256')==sha(Path(__file__)) and release.get('registration_sha256')==sha(P/'checker_registration.json')):
 raise SystemExit('104 saved-check release missing or mismatched')
if any((P/n).exists()for n in ['checked_saved.json','checker_receipt.json','checker_started.json']):raise SystemExit('one saved-check attempt only')
started=time.monotonic();start_cpu=time.process_time();checks=[]
start_utc=datetime.datetime.now(datetime.timezone.utc).isoformat()
(P/'checker_started.json').write_text(json.dumps({'pid':os.getpid(),'started_utc':start_utc,'default177_fuse_active':True,'children':0,'scientific_generator_calls':0},indent=2)+'\n')
def I(d):
 e=int(d['common_binary_exponent']);a=Q(int(d['lower_mantissa']))*Q(2)**e;b=Q(int(d['upper_mantissa']))*Q(2)**e
 assert a<=b
 return a,b
def overlaps(a,b):return max(a[0],b[0])<=min(a[1],b[1])
def scaled(a,x):assert x>=0;return a[0]*x,a[1]*x
def power(a,n):assert a[0]>=0;return a[0]**n,a[1]**n
def all_dyadics(obj):
 if isinstance(obj,dict):
  if {'lower_mantissa','upper_mantissa','common_binary_exponent'}<=set(obj):I(obj)
  for v in obj.values():all_dyadics(v)
 elif isinstance(obj,list):
  for v in obj:all_dyadics(v)
def check():
 cr=json.loads((P/'checker_registration.json').read_text())
 assert sha(P/'saved_check_binding.json')==cr['binding_sha256']
 assert sha(Path(__file__))==cr['checker_script_sha256']
 bind=json.loads((P/'saved_check_binding.json').read_text())
 for name,h in bind['local_sha256'].items():assert sha(P/name)==h
 pre=json.loads((P/'pre_run_manifest.json').read_text())
 for row in pre['files']:assert sha(P/row['name'])==row['sha256']
 checks.append('all bound local saved/pre-run bytes unchanged')
 r=json.loads((P/'results.json').read_text());rec=json.loads((P/'run_receipt.json').read_text());rel=json.loads((P/'parent_release.json').read_text())
 tool=json.loads((P/'execution_evidence/actual_tool_receipt.json').read_text())
 gate=json.loads((P/'execution_evidence/release_receipt.json').read_text())
 up=json.loads((P/'upstream_constants.json').read_text())
 assert r['source_code_sha256']==rel['script_sha256']==sha(P/'run_prepared.py')
 assert r['registration_sha256']==rel['registration_sha256']==sha(P/'registration.json')
 assert r['upstream_constants_sha256']==sha(P/'upstream_constants.json')
 assert rec['release_sha256']==gate['release_sha256']==sha(P/'parent_release.json')
 assert rel['pre_run_manifest_sha256']==sha(P/'pre_run_manifest.json')
 assert rel['authorized'] and rel['scientific_executions_before_release']==0
 assert rec['status']=='COMPLETED' and rec['completed_configurations']==4 and rec['failure']is None
 assert rec['SIGALRM_handler_is_default'] and rec['fuse_not_cancelled']
 assert Q(str(rec['fuse_remaining_before_final_receipt_write']))>0
 assert Q(str(rec['CPU_process_total_before_final_receipt_write']))<60
 assert rec['children_created']==0 and rec['retries']==0
 assert tool['tool_result']['exit_code']==0 and tool['tool_result']['chunk_id']=='d5b280'
 out=json.loads(tool['tool_result']['output'])
 assert out['status']=='COMPLETED' and out['completed_configurations']==4 and out['pid']==rec['pid']==42765
 assert tool['scientific_invocations']==1 and tool['native_calls']==0
 assert 0<Q(str(tool['full_tool_roundtrip_seconds']))<180
 assert gate['only_hash_and_AST_static_checks'] and gate['science_calls']==0
 checks.append('unique producer release/PID/exit/default fuse and complete outer roundtrip evidence consistent')
 all_dyadics(r)
 cfg={x['configuration']:x for x in r['configurations']};assert sorted(cfg)==[1,2,3,4]
 assert cfg[1]['ordered_pair_counts']==[180,450,300] and sum(cfg[1]['ordered_pair_counts'])==cfg[1]['pair_count']==930
 assert cfg[1]['J2_squared']==[20,32,36]
 assert all(I(x['all930_pair_distance_lower'])[0]>Q(1,64)for x in cfg[1]['all_pair_gaps'])
 assert cfg[1]['Crem_missing_peak_loss_and_third_deletion_used']is False
 checks.append('saved complete930-pair structure and every recorded positive gap')
 core=cfg[2];full=cfg[3]
 assert overlaps(power(I(core['I_one_core']),5),power(I(core['full_G_core_norm']),6))
 assert overlaps(I(core['I_ten_cores']),scaled(I(core['I_one_core']),10))
 ratio=I(core['ratio_to_saved100_same_core_norm']);assert 0<ratio[0]<=ratio[1]<1
 n100=I(up['norm100_core'])[1]
 assert overlaps(scaled(ratio,n100),I(core['full_G_core_norm']))
 assert I(core['S_paid_CEH_plus_core_other_unpaid_zero'])[0]>0
 paid=Q(up['paid_CEH_I_upper']);target5=Q(64)**5*(Q(1,65000)/6)**6
 assert (paid+I(core['I_ten_cores'])[1])**5<target5
 checks.append('core saved power/ten-region/ratio consistency and conditional remaining allowance sign')
 rows=full['cases'];assert [x['k']for x in rows]==[1,2,3,4,5]
 assert [x['multiplicity']for x in rows]==[1,5,10,10,5]
 for row in rows:
  assert overlaps(power(I(row['weighted_I']),5),scaled(power(I(row['full_norm_one']),6),Q(row['multiplicity'])**5))
 totals=(sum((I(x['weighted_I'])[0]for x in rows),Q(0)),sum((I(x['weighted_I'])[1]for x in rows),Q(0)))
 assert overlaps(totals,I(full['new_middle_I_upper']))
 newI=I(full['new_middle_I_upper']);oldI=I(full['old051short_plus074long_same_middle_I_upper_expression'])
 assert newI[0]>oldI[1]>0
 quot=(newI[0]/oldI[1],newI[1]/oldI[0]);assert overlaps(quot,I(full['new_over_old_middle_expression']))
 assert (paid+newI[0])**5>target5
 assert I(full['S_new_plan_far_optimistically_zero'])[1]<0
 assert I(full['S_old_plan_far_optimistically_zero'])[1]<0
 assert full['g0_status']=='unproved globally' and full['actual_I_or_C6_lower_bound']is False
 checks.append('five-type saved power/sum/ratio consistency and negative complete-middle allocation diagnostic')
 assert cfg[4]['three_nonparallel_normal_contributions_sum_zero'] and cfg[4]['no_two_normal_contributions_cancel']
 assert r['Cn_improved']is False and r['new_native_kernel_peak_grid_replays']==0 and r['children']==0
 checks.append('scope flags preserve no-C6 claim; triangle witness flags bound by producer hash')
 return {'status':'PASS_LIMITED_SAVED_BINDING_AND_CONSISTENCY_CHECK','checks':checks,'scope':'saved hashes, reported execution, exact dyadic signs and power/addition consistency only; no proof or log/radial/norm-formula regeneration','generator_invocations':0,'new_scientific_configurations':0,'body_kernel_peak_grid_calls':0,'proof_audit_separate':True,'global_g_I_S_unpaid':True}
status='FIRST_FAILURE_STOP';failure=None
try:
 checked=check();(P/'checked_saved.json').write_text(json.dumps(checked,indent=2)+'\n');status='COMPLETED'
except Exception:failure=traceback.format_exc()
finally:
 receipt={'status':status,'pid':os.getpid(),'started_utc':start_utc,'receipt_constructed_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'elapsed_before_final_receipt_write':time.monotonic()-started,'CPU_process_total_before_final_receipt_write':time.process_time(),'default_SIGALRM_still_active':signal.getsignal(signal.SIGALRM)==signal.SIG_DFL,'remaining_fuse_seconds_before_final_receipt_write':signal.getitimer(signal.ITIMER_REAL)[0],'fuse_never_cancelled':True,'fuse_scope':'continues through receipt/print/normal exit','children':0,'retries':0,'generator_invocations':0,'checks_completed':checks,'failure':failure}
 (P/'checker_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
 print(json.dumps({'status':status,'saved_checks_completed':len(checks),'generator_invocations':0,'default_fuse_active':receipt['default_SIGALRM_still_active']}))
if status!='COMPLETED':raise SystemExit(1)
