"""Portable round042 checks. Default: no training/full Monte Carlo replay."""
from pathlib import Path
from fractions import Fraction as F
import argparse,hashlib,json,math,os,shutil,subprocess,sys,tempfile,time
HERE=Path(__file__).resolve().parent
SKIP={'runtime','python','python_version','numpy','numpy_version','seconds','elapsed_seconds'}
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def compare(a,b):
 if isinstance(a,dict):
  assert set(a)-SKIP==set(b)-SKIP,(set(a)-set(b),set(b)-set(a))
  for k in set(a)-SKIP:
   if not k.endswith('_seconds'):compare(a[k],b[k])
 elif isinstance(a,list):
  assert len(a)==len(b)
  for x,y in zip(a,b):compare(x,y)
 elif isinstance(a,float):assert math.isclose(a,b,rel_tol=2e-7,abs_tol=2e-9),(a,b)
 else:assert a==b,(a,b)
def mass_checks():
 result=[];src=HERE/'round042_selector'
 for folder in [src,src/'balanced/n7',src/'balanced/n15']:
  plan=json.loads((folder/'PLAN.json').read_text());f=json.loads((folder/'FROZEN_SELECTOR.json').read_text());c=plan['config'];n=c['n']
  assert f['plan_SHA256']==sha(folder/'PLAN.json')
  assert (plan['script_SHA256']if n==8 else plan['core_SHA256'])==sha(src/'probe.py')
  if n!=8:assert plan['driver_SHA256']==sha(src/'balanced.py')
  index=dict(zip(map(tuple,f['types']),f['selector']));pi=[F(0)for _ in range(4)];orbits=set()
  assert len(index)==math.comb(n+3,3)
  for a,b in index.items():
   assert len(a)==4 and all(isinstance(v,int)and v>=0 for v in a)and sum(a)==n and b in range(4)
   pi[b]+=F(math.factorial(n),math.prod(math.factorial(v)for v in a)*4**n)
   if n%2:
    orb=[a[r:]+a[:r]for r in range(4)];assert len(set(orb))==4
    assert all(index[t]==(b+r)%4 for r,t in enumerate(orb));orbits.add(min(orb))
  assert sum(pi)==1 and list(map(float,pi))==f['prior_mass']
  if n%2:assert pi==[F(1,4)]*4
  result.append(dict(n=n,types=len(index),exact_mass='1',exact_pi=list(map(str,pi)),orbits=len(orbits)if n%2 else None))
 return result

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--output',required=True);ap.add_argument('--replay-cells',action='store_true');args=ap.parse_args()
 outpath=Path(args.output);assert not outpath.exists();start=time.time()
 prov=json.loads((HERE/'provenance.json').read_text())
 for name,h in prov['payload_SHA256'].items():assert sha(HERE/name)==h,name
 masses=mass_checks();runs=[]
 with tempfile.TemporaryDirectory(prefix='cube-round042-')as temp:
  root=Path(temp)
  for folder in ['round042_channel','round042_selector','round042_audit']:shutil.copytree(HERE/folder,root/folder)
  def run(folder,script,output,compare_expected=True):
   target=root/folder/output;expected=json.loads((HERE/folder/output).read_text()) if compare_expected else None
   target.unlink(missing_ok=True)
   proc=subprocess.run([sys.executable,str(root/folder/script)],cwd=root/folder,
    env=dict(os.environ,PYTHONDONTWRITEBYTECODE='1'),capture_output=True,text=True)
   if proc.returncode:raise RuntimeError(script+'\n'+proc.stdout+'\n'+proc.stderr)
   got=json.loads(target.read_text())
   if compare_expected:compare(expected,got)
   runs.append(dict(script=folder+'/'+script,recomputed=True,compared=compare_expected))
   return got
  run('round042_channel','probe.py','SUMMARY.json')
  run('round042_audit','audit.py','AUDIT.json')
  point=run('round042_audit','INTERVAL_CHECK.py','INTERVAL_CHECK.json')
  exact=run('round042_selector','exact_ci.py','EXACT_CI.json')
  # Runtime metadata makes a regenerated JSON a different byte sequence. The
  # independent frozen audit binds the original bytes already hash-checked above.
  (root/'round042_selector/EXACT_CI.json').write_bytes((HERE/'round042_selector/EXACT_CI.json').read_bytes())
  cell=run('round042_audit','CELL_CI_CHECK.py','CELL_CI_CHECK.json')
  if args.replay_cells:
   (root/'round042_selector/CELL_PLAN.json').unlink()
   replay=run('round042_selector','cell_recheck.py','CELL_RESULTS.json')
   # Integer endpoints are recomputed against these replayed counts, too.
   (root/'round042_selector/EXACT_CI.json').unlink()
   proc=subprocess.run([sys.executable,str(root/'round042_selector/exact_ci.py')],cwd=root,
     env=dict(os.environ,PYTHONDONTWRITEBYTECODE='1'),capture_output=True,text=True)
   assert proc.returncode==0,proc.stdout+proc.stderr
 out=dict(status='PASS',payload_hash_count=len(prov['payload_SHA256']),mass_checks=masses,runs=runs,
  model_training_reexecuted=False,full_query_cell_replay_reexecuted=args.replay_cells,
  exact_integer_CI_comparisons=cell['independent_integer_comparisons'],
  local_fraction_kernel_checks=point['cases'][0]['kernel_fraction_enclosures']+sum(x['kernel_fraction_enclosures']for x in point['cases'][1:]),
  cell_corner_kernel_checks=cell['exact_kernel_corner_checks'],
  new_independent_samples=0,PRNG_iid_proved=False,Cn_log_gain_proved=False,
  scope='Hash, exact normalization, deterministic channel re-integration, representative rational interval checks and exact CI certificates; full MC replay only if explicitly requested.',
  python=sys.version,seconds=time.time()-start)
 outpath.write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
if __name__=='__main__':main()
