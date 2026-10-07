#!/usr/bin/env python3
"""Reproduce shared-source diagnostics and rational local obstruction checks."""
from pathlib import Path
import hashlib,json,os,subprocess,sys,tempfile
HERE=Path(__file__).resolve().parent;PROJECT=HERE.parents[1]
def main():
 env=dict(os.environ,PYTHONDONTWRITEBYTECODE='1')
 with tempfile.TemporaryDirectory(prefix='euclidean55_') as folder:
  out=Path(folder)/'shared.json'
  subprocess.run([sys.executable,str(HERE/'shared_probe.py'),'--output',str(out)],check=True,capture_output=True,text=True,env=env)
  data=json.loads(out.read_text())
  cert=json.loads(subprocess.run([sys.executable,str(HERE/'gram_certificate.py')],check=True,capture_output=True,text=True,env=env).stdout)
  assert data['status']==cert['status']=='passed'
  assert data['batch_counts']==[7,3,6] and len(data['cases'])==16
  assert data['total_positive_cells']==2749 and data['total_independent_area_checks']==1731
  keys=['batch','label','alpha','D','O','O_over_X','T_sigma_over_E','O_over_alpha_Tsigma','sigma_min','sigma_max','tau_min','tau_max','sigma_cost','sigma_requested_groups','sigma1_theta1_cost','E_sigma_endpoint_entropy','S_sigma_square','positive_cells','three_independent_area_checks']
  compact=[{k:c[k] for k in keys} for c in data['cases']]
  depnames=['041/square_probe.py','042/row_probe.py','043/global_probe.py','044/tail_probe.py','046/prefix_probe.py','046/threshold_probe.py','047/flow_probe.py','048/shear_probe.py','049/global_probe.py','053/obstruction_probe.py','054/lag_probe.py','054/verification.json']
  deps=sorted((PROJECT/'runtime').rglob('*.py'))+[PROJECT/'rounds'/p for p in depnames]+sorted(HERE.glob('*.py'))
  receipt=dict(status='passed',round=55,date='2026-10-08',python=sys.version.split()[0],shared_tests=dict(batch_counts=data['batch_counts'],positive_cells=data['total_positive_cells'],independent_area_checks=data['total_independent_area_checks'],round54_reference_cases=15,new_independent_full_E_case=1,cases=compact),smooth_local_checks=data['local_smooth_sigma1_theta1_certificates'],gram_certificate=cert,scope='Exact finite n=1 atomic integration; local smooth interval certificates; separate analytic all-k proof. Entropy diagnostics use Decimal, not interval enclosures.',main_weak_type_theorem_proved=False,actual_outer_budget_proved=False,sha256={str(p.relative_to(PROJECT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in deps})
  print(json.dumps(receipt,ensure_ascii=False,indent=2))
if __name__=='__main__':main()
