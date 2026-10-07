#!/usr/bin/env python3
"""Reproduce round50 using standard-library exact arithmetic and scalar Decimal checks."""
from pathlib import Path
from fractions import Fraction as F
import hashlib,json,os,subprocess,sys,tempfile
sys.dont_write_bytecode=True
HERE=Path(__file__).resolve().parent;PROJECT=HERE.parents[1]
def main():
 with tempfile.TemporaryDirectory(prefix='euclidean50_') as folder:
  env=dict(os.environ,PYTHONDONTWRITEBYTECODE='1');out=Path(folder)/'activation.json'
  subprocess.run([sys.executable,str(HERE/'activation_probe.py'),'--output',str(out)],check=True,capture_output=True,text=True,env=env)
  data=json.loads(out.read_text())
  counter=json.loads(subprocess.run([sys.executable,str(HERE/'backward_probe.py')],check=True,capture_output=True,text=True,env=env).stdout)
  assert data['status']==counter['status']=='passed'
  assert data['case_count']==12 and data['batch_counts']==[4,3,5]
  assert counter['batch_counts']==[3,3,3] and counter['strict_clock_extreme_checks']==288
  assert counter['activation_scalar_checks']==27 and len(counter['deep_factor_checks'])==9
  assert all(c['independent_flow_R_B_check']=='passed' for c in data['cases'])
  assert sum(c.get('independent_round49_CG_crosscheck')=='passed' for c in data['cases'])==11
  assert max(F(c['C_D_over_X']) for c in data['cases'])==F(512639,16777216)
  assert max(F(c['T_over_X']) for c in data['cases'])==F(297151,8388608)
  keys=['batch','label','N','D','alpha','X','C_D_over_X','T_over_X','R_over_X','B_over_X','T_over_R','extra_T_over_X','activation_first_moment_over_X','C_D_later_old_label','deep_budget_check','outer_identity_check','retained_bridge_checks','independent_flow_R_B_check','exact_error']
  dep=sorted((PROJECT/'runtime').rglob('*.py'))+[PROJECT/('rounds/'+p) for p in ['041/square_probe.py','042/row_probe.py','043/global_probe.py','044/tail_probe.py','046/prefix_probe.py','046/threshold_probe.py','047/flow_probe.py','048/shear_probe.py','049/global_probe.py','049/verification.json']]+sorted(HERE.glob('*.py'))
  receipt=dict(status='passed',round=50,date='2026-10-07',python=sys.version.split()[0],arithmetic='Fraction exact spatial and threshold integrals; Decimal scalar checks are not interval certificates',activation=dict(case_count=12,batch_counts=data['batch_counts'],round49_crosschecks=11,independent_flow_crosschecks=12,counts={k:sum(c[k] for c in data['cases']) for k in ['T_terms','target_pieces','three_record_intersections']},cases=[{k:c[k] for k in keys} for c in data['cases']]),backward=counter,main_theorem_proved=False,general_T_bound_proved=False,sha256={str(p.relative_to(PROJECT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in dep})
  print(json.dumps(receipt,ensure_ascii=False,indent=2))
if __name__=='__main__':main()
