"""Portable evidence replay: copied frozen scripts, compact parts, zero Green nodes."""
from pathlib import Path
import tempfile,shutil,subprocess,json,sys,hashlib
OUT=Path(__file__).resolve().parent
STAGE=OUT/'cb11_integral_cheb_final_staging'
def run():
 for script,expected in [('cb11_integral_cheb_check.py','4f89a377fbdf263061c2d6525d2101a773f7a6ba4127bb367bd3b7328ed57eb5'),('cb11_integral_cheb_selfcheck.py','368ef2a0705f6adcdc606f0621a46f9921eab605487b25d384e7b1d9e7dd6dcd'),('cb11_cones.json','23258e3d0236700ff442d602ccbb738a17c92728f3b40d95ddf9dd5d9ef924f7')]:
  assert hashlib.sha256((STAGE/script).read_bytes()).hexdigest()==expected
 with tempfile.TemporaryDirectory(prefix='cb11_integral_cheb_final_verify_',dir=OUT) as tmp:
  tmp=Path(tmp)
  for name in ['cb11_integral_cheb_check.py','cb11_integral_cheb_selfcheck.py','cb11_cones.json']+[f'cb11_integral_cheb_part_{p:02}.json' for p in range(12)]:shutil.copy2(STAGE/name,tmp/name)
  for script,args in [('cb11_integral_cheb_selfcheck.py',['--parts','12']),('cb11_integral_cheb_check.py',['--stage','aggregate','--parts','12'])]:subprocess.run([sys.executable,str(tmp/script),*args],check=True,capture_output=True,text=True)
  old=json.loads((STAGE/'cb11_integral_cheb_all_result.json').read_text());new=json.loads((tmp/'cb11_integral_cheb_all_result.json').read_text())
  fields=['status','I_certified_locally','config','patch_count','IH_upper_exact','IH_upper','Etime','I_upper','I_upper_display','triangle_weighted_area_exact','checker_sha256']
  assert all(old[k]==new[k] for k in fields)
  selfcheck=json.loads((tmp/'cb11_integral_cheb_all_selfcheck.json').read_text());assert selfcheck['patches_checked']==8256 and selfcheck['full_coverage'] and selfcheck['new_point_evaluations']==0
 print(json.dumps({'status':'PASS compact8256 exact chord/coverage/aggregate replay','new_point_evaluations':0,'IH_upper_exact':new['IH_upper_exact'],'I_upper_packed':new['I_upper'],'I_upper_display':new['I_upper_display'],'independent_audit':'this is author evidence replay, not independent full-node acceptance'},indent=2))
if __name__=='__main__':run()
