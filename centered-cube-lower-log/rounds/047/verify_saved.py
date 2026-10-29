from pathlib import Path
import json,hashlib,math,os,shutil,subprocess,sys,tempfile,argparse,time
import numpy as np
H=Path(__file__).resolve().parent
SKIP={'elapsed_seconds','elapsed_main_v2_seconds','input_sha256','local_full_result_file_sha256'}
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def cmp(a,b):
 if isinstance(a,dict):
  assert set(a)-SKIP==set(b)-SKIP
  for k in set(a)-SKIP:cmp(a[k],b[k])
 elif isinstance(a,list):
  assert len(a)==len(b)
  for x,y in zip(a,b):cmp(x,y)
 elif isinstance(a,float):assert math.isclose(a,b,rel_tol=1e-8,abs_tol=2e-10),(a,b)
 else:assert a==b,(a,b)
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,required=True);args=ap.parse_args();assert not args.output.exists();st=time.time();prov=json.loads((H/'provenance.json').read_text());runs=[]
 for n,h in prov['payload_SHA256'].items():assert sha(H/n)==h,n
 for n,h in prov['reference_SHA256'].items():assert sha(H.parent/n)==h,n
 with tempfile.TemporaryDirectory(prefix='cube47-')as td:
  root=Path(td)
  for d in ['round047_complex','round047_profile','round047_growth']:shutil.copytree(H/d,root/d)
  def run(d,s,out,compare=True):
   p=subprocess.run([sys.executable,str(root/d/s)],cwd=root/d,env=dict(os.environ,PYTHONDONTWRITEBYTECODE='1'),capture_output=True,text=True,timeout=180)
   if p.returncode:raise RuntimeError(s+'\n'+p.stderr[-6000:])
   a=json.loads((root/d/out).read_text())
   if compare:cmp(json.loads((H/d/out).read_text()),a)
   runs.append(d+'/'+s);return a
  c=run('round047_complex','complex_probe.py','RESULTS.json');p=run('round047_profile','probe.py','PROFILE.json')
  g=run('round047_growth','run_growth.py','RESULTS.json',False)
  ix=json.loads((H/'round047_growth/PROFILE_INDEX.json').read_text())
  assert hashlib.sha256(json.dumps(g['rows'],sort_keys=True,separators=(',',':')).encode()).hexdigest()==ix['canonical_rows_sha256']
  check=run('round047_growth','verify_and_report.py','CHECKS.json')
  run('round047_growth','package_minimal.py','PROFILE_INDEX.json',False)
  with np.load(H/'round047_growth/PROFILE_MIN.npz')as a,np.load(root/'round047_growth/PROFILE_MIN.npz')as b:
   assert a.files==b.files
   for k in a.files:assert np.array_equal(a[k],b[k]),k
  out=dict(status='PASS',runs=runs,complex_anchor_cases=len(c['anchor_normalized_cases']),complex_log_cases=len(c['complex_log_cases']),single_profile_v_cases=len(p['rows']),growth_configurations=len(g['rows']),all_profile_points=ix['profile_points'],all_profile_arrays_exact=True,canonical_growth_rows_exact=True,deterministic_target_cap=check['deterministic_target_cap'],independent_selected_Hoeffding_radius=check['independent_test_Hoeffding_radius'],density_cases=check['independent_density_checks'],max_density_error=check['max_absolute_log_density_error'],max_R1_empirical_SE_z=check['max_R1_profile_absolute_z'],new_independent_queries=0,actual_HL_growth_input=False,new_Cn_log_gain_proved=False,seconds=time.time()-st)
 args.output.write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
if __name__=='__main__':main()
