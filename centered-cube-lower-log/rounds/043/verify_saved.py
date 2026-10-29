"""Portable round043 audit; optional all-frozen-test replay, no retraining."""
from pathlib import Path
import argparse,hashlib,json,math,os,shutil,subprocess,sys,tempfile,time
HERE=Path(__file__).resolve().parent
SKIP={'seconds','python','python_version','numpy','numpy_version'}
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def compare(a,b):
 if isinstance(a,dict):
  assert set(a)-SKIP==set(b)-SKIP,(set(a)-set(b),set(b)-set(a))
  for k in set(a)-SKIP:
   if not k.endswith('_seconds'):compare(a[k],b[k])
 elif isinstance(a,list):
  assert len(a)==len(b)
  for x,y in zip(a,b):compare(x,y)
 elif isinstance(a,float):assert math.isclose(a,b,rel_tol=2e-8,abs_tol=2e-9),(a,b)
 else:assert a==b,(a,b)
def prepare(root):
 for folder in ['round043_tilt','round043_tilt_audit','round043_binary']:shutil.copytree(HERE/folder,root/folder)
 helper=HERE.parent/'042/round042_channel/probe.py';plan=json.loads((HERE/'round043_tilt/PLAN.json').read_text())
 assert sha(helper)==plan['channel42_SHA256']
 target=root/'round042_channel/probe.py';target.parent.mkdir();shutil.copy2(helper,target)
 hashes=json.loads((HERE/'shared_binary_result_hashes.json').read_text());mapped=[]
 for prefix in ['', 'BALANCED_']:
  summary=json.loads((root/'round043_binary'/f'{prefix}SUMMARY.json').read_text())
  for row in summary['results']:
   name=f'{prefix}RESULT_n{row["n"]}.json';p=root/'round043_binary'/name
   p.write_text(json.dumps(row,indent=2)+'\n');assert sha(p)==hashes[name],name;mapped.append(name)
 return mapped

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--output',required=True);ap.add_argument('--replay-test',action='store_true');args=ap.parse_args()
 outpath=Path(args.output);assert not outpath.exists();started=time.time()
 # The binary pilot's CP arithmetic needs SciPy; do not silently omit that check.
 import numpy,scipy
 prov=json.loads((HERE/'provenance.json').read_text())
 for name,h in prov['payload_SHA256'].items():assert sha(HERE/name)==h,name
 runs=[]
 with tempfile.TemporaryDirectory(prefix='cube-round043-')as td:
  root=Path(td);mapped=prepare(root)
  def run(folder,script,output):
   target=root/folder/output;expected=json.loads((HERE/folder/output).read_text());target.unlink(missing_ok=True)
   p=subprocess.run([sys.executable,str(root/folder/script)],cwd=root/folder,
    env=dict(os.environ,PYTHONDONTWRITEBYTECODE='1'),capture_output=True,text=True)
   if p.returncode:raise RuntimeError(script+'\n'+p.stdout+'\n'+p.stderr)
   got=json.loads(target.read_text());compare(expected,got);runs.append(dict(script=folder+'/'+script,matched=True));return got
  (root/'round043_tilt/PLAN.json').unlink()
  run('round043_tilt','probe.py','SUMMARY.json')
  assert sha(root/'round043_tilt/PLAN.json')==sha(HERE/'round043_tilt/PLAN.json')
  run('round043_tilt','normalization_probe.py','NORMALIZATION.json')
  run('round043_tilt','exact_tilt_example.py','EXACT_TILT_EXAMPLE.json')
  run('round043_tilt','price_probe.py','PRICE.json')
  audit=run('round043_tilt_audit','checks.py','CHECKS.json')
  binary=run('round043_binary','verify_saved.py','VERIFICATION.json')
  run('round043_binary','combine.py','COMBINED_SUMMARY.json')
  if args.replay_test:run('round043_binary','replay_guards.py','GUARD_REPLAY.json')
 out=dict(status='PASS',payload_hash_count=len(prov['payload_SHA256']),runs=runs,
  shared_duplicate_RESULT_files_mapped_with_exact_SHA=mapped,
  numerical_query_laws_reexecuted=audit['numerical_replay_count'],
  exact_physical_source_recomputed=True,new_prior_trace_bound_proved=True,
  original_prior_price_not_omitted=True,model_training_reexecuted=False,
  full_frozen_guard_test_reexecuted=args.replay_test,new_independent_samples=0,
  binary_continuum_or_full_FP_certificate=False,PRNG_iid_proved=False,Cn_log_gain_proved=False,
  python=sys.version,numpy=numpy.__version__,scipy=scipy.__version__,seconds=time.time()-started)
 outpath.write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
if __name__=='__main__':main()
