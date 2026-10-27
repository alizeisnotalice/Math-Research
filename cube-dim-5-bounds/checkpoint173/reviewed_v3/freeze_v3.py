"""Preserve v2, freeze minimal v3 and the sole new integrated test receipt."""
from pathlib import Path
import ast,difflib,hashlib,json,os,time
H=Path(__file__).resolve().parent;V2=H.parent/'reviewed_v2';W=Path('/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/upper_bandcurvature173')
def sha(p):
 h=hashlib.sha256()
 with open(p,'rb') as f:
  while b:=f.read(8*1024*1024):h.update(b)
 return h.hexdigest()
def write(p,x):
 assert not p.exists();p.write_text(json.dumps(x,sort_keys=True,separators=(',',':'))+'\n')
def main():
 began=time.monotonic();assert not (H/'prepared_manifest.json').exists() and not (W/'run001').exists() and not (W/'run_v2_001').exists() and not (W/'run_v3_001').exists();old=V2/'prepared_manifest.json';assert sha(old)=='db85430d700ddefb04c002c899ada1d1aed05547f3777497649ffef793340448';m=json.load(open(old))
 for n,h in m['repository_files'].items():assert sha(V2/n)==h
 for p,h in {**m['source_locks'],**m['prepared_external_files']}.items():assert sha(p)==h
 for n in ('core173.py','worker173.py','util173.py','fault_cleanup173.py'):assert sha(H/n)==sha(V2/n)
 diff=''.join(difflib.unified_diff((V2/'run173.py').read_text().splitlines(keepends=True),(H/'run173.py').read_text().splitlines(keepends=True),fromfile='v2/run173.py',tofile='v3/run173.py'));p=H/'run173.py.vs_v2.diff';assert not p.exists();p.write_text(diff)
 for p in H.glob('*.py'):ast.parse(p.read_text(),filename=str(p))
 T=W/'reviewed_v3_synthetic';r=json.load(open(T/'test_receipt.json'));assert r['status']=='ONE_NEW_ZERO_H_INTEGRATED_REGRESSION_PASSED' and r['master_actual_exitcode']==2 and r['summaries_at_cleanup_entry']==1 and r['did_not_wait_for_all15_summaries'] and r['all15_workers_exited_OS_confirmed'] and not r['old_v2_three_toys_repeated'] and r['new_real_H_jet_donor_band_q_corner_calls']==0;O=T/'fault_report_and_stuck_peer/work/run_v3_001';cleanup=json.load(open(O/'actual_attempt0001_mass_fault_cleanup.json'));assert cleanup['all_joined'] and len(cleanup['workers'])==15
 for w in cleanup['workers']:
  assert not w['alive'] and w['actual_join_called']
  try:os.kill(w['pid'],0);raise AssertionError('synthetic worker remains alive')
  except ProcessLookupError:pass
 synthetic={str(T/'test_receipt.json'):sha(T/'test_receipt.json')}
 for n in ['actual_attempt0001_mass_exit_receipt.json','actual_attempt0001_mass_fault_cleanup.json','actual_attempt0001_summary.json','actual_attempt0001_mass_pool_registration.json']:synthetic[str(O/n)]=sha(O/n)
 locks={**m['source_locks'],str(old):sha(old)}
 for n,h in m['repository_files'].items():locks[str(V2/n)]=h
 receipt={'pid':os.getpid(),'v2_prepared_bytes_unchanged':True,'v3_core_worker_util_cleanup_identical_to_v2':True,'actual_math_namespace_absent':True,'sole_new_synthetic_test':'fault summary plus stop/SIGTERM-ignoring peer;actual run_stage immediate finally','received_summaries_before_cleanup':1,'all15_synthetic_workers_actual_join_and_OS_absence':True,'master_test_exitcode':2,'new_real_H_jet_donor_band_q_corner_calls':0,'normal_numeric_wall_cap':None,'wall_seconds':time.monotonic()-began};rp=W/'round020/v3_prepared_validation.json';write(rp,receipt)
 manifest={'schema':'173_reviewed_v3_prepared','repository_files':{p.name:sha(p) for p in sorted(H.iterdir()) if p.is_file()},'source_locks':locks,'prepared_external_files':m['prepared_external_files'],'inherited_v2_prepared_sha256':sha(old),'inherited_original_prepared_sha256':m['inherited_original_prepared_sha256'],'known_qualified168':m['known_qualified168'],'sources':455,'weighted_band_terms':22295,'fixed_pairs':35872,'receivers':1120,'max_concurrent_workers':15,'output_namespace':str(W/'run_v3_001'),'sole_new_synthetic_not_math_or_qualification':synthetic,'v3_validation':{'path':str(rp),'sha256':sha(rp)},'dispatch_authorized':False,'new_real_H_jet_donor_band_q_corner_calls':0};write(H/'prepared_manifest.json',manifest);print(json.dumps({'v3_manifest_sha256':sha(H/'prepared_manifest.json'),'v3_run_sha256':sha(H/'run173.py'),'v3_diff_sha256':sha(H/'run173.py.vs_v2.diff'),'new_synthetic_receipt_sha256':sha(T/'test_receipt.json'),'v3_validation_sha256':sha(rp),'old_v2_manifest_unchanged':sha(old)},sort_keys=True),flush=True)
if __name__=='__main__':main()
