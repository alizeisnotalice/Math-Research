#!/usr/bin/env python3
"""Run exact lag integrations and constant certificates; emit compact receipt."""
from pathlib import Path
import hashlib,json,os,subprocess,sys,tempfile
HERE=Path(__file__).resolve().parent;PROJECT=HERE.parents[1]
def main():
 env=dict(os.environ,PYTHONDONTWRITEBYTECODE='1')
 with tempfile.TemporaryDirectory(prefix='euclidean54_') as folder:
  out=Path(folder)/'lag.json'
  subprocess.run([sys.executable,str(HERE/'lag_probe.py'),'--output',str(out)],check=True,capture_output=True,text=True,env=env)
  lag=json.loads(out.read_text())
  geo=json.loads(subprocess.run([sys.executable,str(HERE/'geometry_certificate.py')],check=True,capture_output=True,text=True,env=env).stdout)
  assert lag['status']==geo['status']=='passed' and lag['batch_counts']==[6,3,6]
  assert len(lag['cases'])==15 and lag['original_round53_cases_checked']==9
  rows=sum(len(c['record_rows']) for c in lag['cases']);breaks=sum(c['breakpoint_checks'] for c in lag['cases'])
  assert rows==1685 and breaks==6852
  compact=[{k:v for k,v in c.items() if k not in ('record_rows','elapsed_seconds','record_row_witness')} for c in lag['cases']]
  deps=sorted((PROJECT/'runtime').rglob('*.py'))+[PROJECT/('rounds/'+p) for p in ['041/square_probe.py','042/row_probe.py','043/global_probe.py','044/tail_probe.py','046/prefix_probe.py','046/threshold_probe.py','047/flow_probe.py','049/global_probe.py','053/verification.json']]+sorted(HERE.glob('*.py'))
  receipt=dict(status='passed',round=54,date='2026-10-08',python=sys.version.split()[0],lag_tests=dict(status='passed',batch_counts=lag['batch_counts'],record_rows=rows,breakpoint_checks=breaks,reference_cases=9,arithmetic='Fraction',cases=compact),geometry=geo,scope='Finite n=1 integration plus exact constants supporting a separate analytic high-dimensional proof. No sampled high-dimensional volume.',main_weak_type_theorem_proved=False,actual_outer_budget_proved=False,sha256={str(p.relative_to(PROJECT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in deps})
  print(json.dumps(receipt,ensure_ascii=False,indent=2))
if __name__=='__main__':main()
