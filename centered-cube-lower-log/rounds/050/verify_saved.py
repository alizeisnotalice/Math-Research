from pathlib import Path
import argparse,hashlib,json,math,os,shutil,subprocess,sys,tempfile,time
H=Path(__file__).resolve().parent
SKIP={'elapsed_seconds','seconds','runtime_python'}
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
 ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,required=True);args=ap.parse_args();assert not args.output.exists();st=time.time();pr=json.loads((H/'provenance.json').read_text());runs=[]
 for n,hh in pr['payload_SHA256'].items():assert sha(H/n)==hh,n
 for n,hh in pr['reference_SHA256'].items():assert sha(H.parent/n)==hh,n
 with tempfile.TemporaryDirectory(prefix='cube50-')as td:
  root=Path(td)
  for d in ['round050_joint_grid','round050_mixture_budget','round050_gate_probe','round050_cutoff']:shutil.copytree(H/d,root/d)
  old=root/'round048_banks';old.mkdir();shutil.copyfile(H.parent/'048/round048_banks/MANIFEST.json',old/'MANIFEST.json')
  def run(d,s):
   p=subprocess.run([sys.executable,str(root/d/s)],cwd=root/d,env=dict(os.environ,PYTHONDONTWRITEBYTECODE='1'),capture_output=True,text=True,timeout=180)
   if p.returncode:raise RuntimeError(s+'\n'+p.stderr[-6000:])
   data=json.loads((root/d/'RESULTS.json').read_text());cmp(json.loads((H/d/'RESULTS.json').read_text()),data);runs.append(d+'/'+s);return data
  joint=run('round050_joint_grid','probe.py');budget=run('round050_mixture_budget','check_budget.py');gate=run('round050_gate_probe','probe.py');cut=run('round050_cutoff','probe.py')
  out=dict(status='PASS',runs=runs,joint_budget_rows=len(joint['budget_rows']),same_endpoint_checks=joint['round048_four_real_banks_same_endpoint_check'],finite_atom_models=len(budget['finite_atom_models']),E_lower=budget['E_lower'],gate_proxy_cases=len(gate['rows']),cutoff_probe_cases=len(cut['rows']),new_independent_queries=0,new_actual_HL_log_gain_proved=False,seconds=time.time()-st)
 args.output.write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
if __name__=='__main__':main()
