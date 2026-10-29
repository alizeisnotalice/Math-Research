"""Round044 exact checks and stored-data replay; never re-solves time-limited MILP."""
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
def source_checks():
 src=HERE/'round044_symmetry';records=[]
 plan=json.loads((src/'PLAN.json').read_text())
 assert plan['script_SHA256']==sha(src/'probe.py')
 assert json.loads((src/'SCAN_PLAN.json').read_text())['script_SHA256']==sha(src/'probe.py')
 for n in [32,64,128]:
  f=json.loads((src/f'FROZEN_n{n}.json').read_text());s=f['selector']
  assert f['plan_SHA256']==sha(src/'PLAN.json')and len(s)==n+1 and all(b in range(4)for b in s)
  nums=[sum(math.comb(n,k)for k,b in enumerate(s)if b==j)for j in range(4)]
  assert nums==f['prior_exact_numerators']and sum(nums)==f['prior_exact_denominator']==2**n
  assert all(math.isclose(a,b/2**n,abs_tol=5e-15)for a,b in zip(f['prior_mass'],nums))
  records.append(dict(n=n,exact_mass='1',prior_numerators=nums,denominator=2**n))
 m=HERE/'round044_milp';f=json.loads((m/'FROZEN.json').read_text());pl=json.loads((m/'PLAN.json').read_text())
 assert pl['script_SHA256']==sha(m/'probe.py')and pl['cost_SHA256']==sha(m/'COST.json')
 assert f['training_SHA256']==sha(m/'TRAINING.json')and f['plan_SHA256']==sha(m/'PLAN.json')
 assert json.loads((m/'SUMMARY.json').read_text())['frozen_SHA256']==sha(m/'FROZEN.json')
 records0=json.loads((m/'TRAINING.json').read_text())
 assert all(r['status']==1 and not r['success']for r in records0)
 best=max(records0,key=lambda r:r['recomputed_train_weak'])
 assert f['selector']==best['selector']and f['lambda_']==best['lambda_']
 nums=[sum(math.comb(32,k)for k,b in enumerate(f['selector'])if b==j)for j in range(4)]
 assert nums==f['prior_exact_numerators']and sum(nums)==f['prior_exact_denominator']==2**32
 from scipy.stats import beta
 summary=json.loads((m/'SUMMARY.json').read_text())
 for row in summary['rows']:
  k,N=row['k'],row['N'];a=.05/12
  expected=[0. if k==0 else summary['lambda_']*float(beta.ppf(a,k,N-k+1)),1.*summary['lambda_']if k==N else summary['lambda_']*float(beta.ppf(1-a,k+1,N-k))]
  assert all(math.isclose(x,y,abs_tol=2e-12)for x,y in zip(expected,row['CP_weak']))
 return records

def prepare(root):
 for folder in ['round044_symmetry','round044_symmetry_audit','round044_milp']:shutil.copytree(HERE/folder,root/folder)
 helper43=HERE.parent/'043/round043_binary/probe.py';helper42=HERE.parent/'042/round042_selector/probe.py'
 expected43=json.loads((HERE/'round044_symmetry/PLAN.json').read_text())['fit_helper_SHA256']
 assert sha(helper43)==expected43==json.loads((HERE/'round044_milp/PLAN.json').read_text())['source_SHA256']
 assert sha(helper42)==json.loads((HERE/'round044_symmetry/JOINT_K4_PLAN.json').read_text())['helper42_SHA256']
 for name,src in [('round043_binary',helper43),('round042_selector',helper42)]:
  dst=root/name/'probe.py';dst.parent.mkdir();shutil.copy2(src,dst)

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--output',required=True);ap.add_argument('--replay-mirror-test',action='store_true');args=ap.parse_args()
 target=Path(args.output);assert not target.exists();start=time.time();import numpy,scipy
 prov=json.loads((HERE/'provenance.json').read_text())
 for name,h in prov['payload_SHA256'].items():assert sha(HERE/name)==h,name
 for name,h in prov['reference_SHA256'].items():assert sha(HERE.parent/name)==h,name
 masses=source_checks();runs=[]
 with tempfile.TemporaryDirectory(prefix='cube-round044-')as td:
  root=Path(td);prepare(root)
  def run(folder,script,output,args0=(),remove_plan=None):
   p=root/folder/output;p.unlink(missing_ok=True)
   if remove_plan:(root/folder/remove_plan).unlink()
   proc=subprocess.run([sys.executable,str(root/folder/script),*args0],cwd=root/folder,
    env=dict(os.environ,PYTHONDONTWRITEBYTECODE='1'),capture_output=True,text=True)
   if proc.returncode:raise RuntimeError(script+'\n'+proc.stdout+'\n'+proc.stderr)
   got=json.loads(p.read_text());compare(json.loads((HERE/folder/output).read_text()),got)
   runs.append(dict(script=folder+'/'+script,matched=True));return got
  run('round044_symmetry','probe.py','SCAN.json',('scan',),'SCAN_PLAN.json')
  run('round044_symmetry','joint_probe.py','JOINT_SCAN.json',remove_plan='JOINT_PLAN.json')
  run('round044_symmetry','joint_k4_probe.py','JOINT_K4.json',remove_plan='JOINT_K4_PLAN.json')
  run('round044_symmetry','scale_budget_probe.py','SCALE_BUDGET.json')
  check=run('round044_symmetry_audit','checks.py','CHECKS.json')
  joint=run('round044_symmetry_audit','JOINT_CHECK.py','JOINT_CHECK.json')
  milp=run('round044_milp','verify_saved.py','VERIFICATION.json')
  run('round044_milp','greedy_control.py','GREEDY_SUMMARY.json')
  assert sha(root/'round044_milp/GREEDY_PLAN.json')==sha(HERE/'round044_milp/GREEDY_PLAN.json')
  compare(json.loads((root/'round044_milp/GREEDY_FROZEN.json').read_text()),json.loads((HERE/'round044_milp/GREEDY_FROZEN.json').read_text()))
  if args.replay_mirror_test:run('round044_symmetry','witness_probe.py','WITNESS.json')
 out=dict(status='PASS',payload_hash_count=len(prov['payload_SHA256']),reference_hash_count=len(prov['reference_SHA256']),
  exact_mirror_frozen_mass_checks=masses,runs=runs,
  exact_mirror_kernel_checks=check['exact_kernel_reflection_checks'],
  exact_joint_channel_integrals=joint['exact_joint_abstract_channel_integrals'],
  MILP_solver_reexecuted=False,all_MILP_frozen_test_reexecuted=True,greedy_control_reexecuted=True,
  mirror_model_training_reexecuted=False,all_mirror_frozen_test_reexecuted=args.replay_mirror_test,
  new_independent_samples=0,full_high_dimensional_FP_or_continuum_certificate=False,Cn_log_gain_proved=False,
  python=sys.version,numpy=numpy.__version__,scipy=scipy.__version__,seconds=time.time()-start)
 target.write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
if __name__=='__main__':main()
