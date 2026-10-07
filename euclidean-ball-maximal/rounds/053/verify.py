#!/usr/bin/env python3
"""Re-run round53 tests; print a compact reproducibility receipt."""
from pathlib import Path
import hashlib,json,os,subprocess,sys,tempfile
HERE=Path(__file__).resolve().parent
PROJECT=HERE.parents[1]
def main():
 env=dict(os.environ,PYTHONDONTWRITEBYTECODE='1')
 with tempfile.TemporaryDirectory(prefix='euclidean53_') as folder:
  out=Path(folder)/'records.json'
  subprocess.run([sys.executable,str(HERE/'record_probe.py'),'--output',str(out)],check=True,capture_output=True,text=True,env=env)
  data=json.loads(out.read_text())
  obs=json.loads(subprocess.run([sys.executable,str(HERE/'obstruction_probe.py')],check=True,capture_output=True,text=True,env=env).stdout)
  assert data['status']==obs['status']=='passed'
  assert data['batch_counts']==[3,3,3] and len(data['cases'])==9
  assert data['total_positive_cells']==84 and data['total_source_pairs']==226
  assert data['total_acceptance_probability_checks']==15812
  assert all(c['independent_round49_O_check']=='passed' for c in data['cases'])
  compact=[]
  for c in data['cases']:
   compact.append({k:v for k,v in c.items() if k not in ('positive_cell_records','elapsed_seconds','ln2_rational_lower','atoms')})
  deps=sorted((PROJECT/'runtime').rglob('*.py'))+[PROJECT/('rounds/'+p) for p in ['041/square_probe.py','042/row_probe.py','043/global_probe.py','044/tail_probe.py','046/prefix_probe.py','046/threshold_probe.py','047/flow_probe.py','049/global_probe.py','049/verification.json']]+sorted(HERE.glob('*.py'))
  receipt=dict(status='passed',round=53,date='2026-10-07',python=sys.version.split()[0],record_tests=dict(batch_counts=[3,3,3],positive_cells=84,source_pairs=226,acceptance_checks=15812,cases=compact,entropy_precision_digits=80,entropy_is_interval_certificate=False),obstructions=obs,main_theorem_proved=False,actual_global_outer_budget_proved=False,sha256={str(p.relative_to(PROJECT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in deps})
  print(json.dumps(receipt,ensure_ascii=False,indent=2))
if __name__=='__main__':main()
