"""Portable Round045 exact-model and frozen-test replay, using temporary copies only."""
from pathlib import Path
import argparse,hashlib,json,math,os,shutil,subprocess,sys,tempfile,time
H=Path(__file__).resolve().parent
SKIP={'seconds','python','numpy','scipy','child_wall_seconds','wrapper_wall_seconds'}
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def compare(a,b):
 if isinstance(a,dict):
  assert set(a)-SKIP==set(b)-SKIP,(set(a),set(b))
  for k in set(a)-SKIP:compare(a[k],b[k])
 elif isinstance(a,list):
  assert len(a)==len(b)
  for x,y in zip(a,b):compare(x,y)
 elif isinstance(a,float):assert math.isclose(a,b,abs_tol=2e-9,rel_tol=2e-8),(a,b)
 else:assert a==b,(a,b)
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,required=True);args=ap.parse_args();assert not args.output.exists()
 start=time.time();p=json.loads((H/'provenance.json').read_text());runs=[]
 for name,hh in p['payload_SHA256'].items():assert sha(H/name)==hh,name
 pl=json.loads((H/'round045_prior/PLAN.json').read_text());assert pl['script_SHA256']==sha(H/'round045_prior/probe.py') and pl['cost_SHA256']==sha(H/'round045_prior/COST.json')
 summary=json.loads((H/'round045_prior/SUMMARY.json').read_text())
 for r in summary['results']:
  frozen=json.loads((H/f"round045_prior/FROZEN_n{r['n']}.json").read_text());assert frozen['plan_SHA256']==sha(H/'round045_prior/PLAN.json') and frozen['pi']==r['pi']
 with tempfile.TemporaryDirectory(prefix='cube-round045-')as td:
  root=Path(td)
  for folder in ['round045_prior','round045_baseline','round045_is_audit','accepted']:shutil.copytree(H/folder,root/folder)
  prior=root/'round045_prior'
  for row in summary['results']:
   name=f"RESULT_n{row['n']}.json";data=json.dumps(row,indent=2)+'\n';assert hashlib.sha256(data.encode()).hexdigest()==p['reconstructed_RESULT_SHA256'][name];(prior/name).write_text(data)
  def run(folder,script,output,arguments=(),expected=True):
   out=root/folder/output;out.unlink(missing_ok=True)
   rr=subprocess.run([sys.executable,str(root/folder/script),*map(str,arguments)],cwd=root/folder,env=dict(os.environ,PYTHONDONTWRITEBYTECODE='1'),capture_output=True,text=True,timeout=180)
   if rr.returncode:raise RuntimeError(script+'\n'+rr.stdout[-2000:]+'\n'+rr.stderr[-6000:])
   got=json.loads(out.read_text())
   if expected:compare(json.loads((H/folder/output).read_text()),got)
   runs.append(dict(script=folder+'/'+script,passed=True));return got
  finite=run('round045_is_audit','finite_probe.py','FINITE.json')
  trap=run('round045_is_audit','mass_trap.py','MASS_TRAP.json')
  ci=run('round045_is_audit','ci_profile_audit.py','CI_PROFILE.json')
  baseline=run('round045_baseline','portable_verify.py','PORTABLE_VERIFY.json',('--paper',root/'accepted/paper.tex','--labels',root/'accepted/source_d10.npz','--output',root/'round045_baseline/PORTABLE_VERIFY.json'),False)
  assert baseline['status']=='PASS'and baseline['byte_identical_result'] and baseline['exact_json_equal'] and baseline['original_inputs_unchanged']
  sampler=run('round045_prior','validate_sampler.py','SAMPLER_VALIDATION.json')
  plateau=run('round045_prior','plateau_replay.py','PLATEAU_REPLAY.json',expected=False)
  assert all(r['guarded_classifications_unchanged']and r['max_weak_change']<1e-12 for r in plateau['results'])
  verifier=run('round045_prior','verify_saved.py','VERIFICATION.json')
  assert not any(prior.glob('NEAR_INDICES*'))
  replay=run('round045_prior','replay_summary.py','REPLAY_SUMMARY.json',expected=False)
  assert replay['passed']and not any(prior.glob('NEAR_INDICES*'))
  rr=[]
  for r in replay['results']:
   assert r['guarded_counts_exact']and r['near_ids_exact']
   rr.append({k:r[k]for k in ['n','near_total','comparison_errors','whole_profile_mixture_max','whole_profile_bestfixed_max','whole_profile_bestfixed_h']})
  out=dict(status='PASS',payload_hashes=len(p['payload_SHA256']),runs=runs,
   finite_exact_cases=len(finite['rows']),finite_wrong_RN_mismatches=sum(r['true_weak']!=r['wrong_per_coordinate_RN_mean']for r in finite['rows']),
   mass_trap_true=trap['true_weak_exact'],mass_trap_wrong=trap['wrong_RN_weak_exact'],
   Decimal80_CI_intervals=ci['models_thresholds_checked'],Decimal80_CI_endpoints=2*ci['models_thresholds_checked'],probability_direction_budget=ci['directions'],Decimal80_max_endpoint_difference=ci['max_EB_endpoint_difference'],
   baseline_parameter_rows=baseline['parameter_rows_replayed'],baseline_result_byte_identical=baseline['byte_identical_result'],source_CDF=baseline['sampler_replayed'],
   rational_kernel_cases=verifier['rational_kernel_cases'],kernel_max_error=verifier['max_kernel_float_error'],frozen_replays=rr,
   exact_canonical_near_hashes=3,near_indices_reconstructed=sum(r['near_total']for r in replay['results']),stored_near_arrays=0,
   training_reexecuted=False,new_independent_test_queries=0,FP_probability_or_continuum_certificate=False,new_Cn_log_gain_proved=False,
   seconds=time.time()-start,python=sys.version,numpy=__import__('numpy').__version__,scipy=__import__('scipy').__version__)
 args.output.parent.mkdir(parents=True,exist_ok=True);args.output.write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
if __name__=='__main__':main()
