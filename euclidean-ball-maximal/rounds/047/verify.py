#!/usr/bin/env python3
"""Run round47 from final code; emit a compact, exact verification receipt."""
from pathlib import Path
from fractions import Fraction as F
import hashlib,json,os,subprocess,sys,tempfile
sys.dont_write_bytecode=True
HERE=Path(__file__).resolve().parent;PROJECT=HERE.parents[1]
def main():
 with tempfile.TemporaryDirectory(prefix='euclidean47_') as folder:
  out=Path(folder)/'flow.json';env=dict(os.environ,PYTHONDONTWRITEBYTECODE='1')
  subprocess.run([sys.executable,str(HERE/'flow_probe.py'),'--output',str(out)],check=True,capture_output=True,text=True,env=env)
  stress=json.loads(out.read_text())
  local=json.loads(subprocess.run([sys.executable,str(HERE/'local_law_probe.py')],check=True,capture_output=True,text=True,env=env).stdout)
  assert stress['status']==local['status']=='passed'
  assert [b['count'] for b in stress['batches']]==[12,16,15]
  assert [b['positive_count'] for b in stress['batches']]==[3,7,10]
  assert local['batch_counts']==[9,9,15] and len(local['law_tests'])==24
  full=[c for c in stress['cases'] if 'record_kernel' in c]
  assert [sum(c['batch']==b for c in full) for b in range(3)]==[2,3,8]
  assert sum(c['record_kernel']['exclusive_second_source_annular_pairs_checked'] for c in full)==427
  assert sum(c['record_energy']['forward_record_collapse_checks'] for c in full)==102751
  assert sum(c['record_energy']['reverse_record_collapse_checks'] for c in full)==102751
  assert sum(c['record_energy']['weighted_symmetric_kernel_pairs_checked'] for c in full)==40033
  negative=next(c for c in stress['cases'] if c['label']=='geometric_contact_clouds_s491101')
  assert negative['beta_negative_flow_probability']=='46/201'
  assert negative['first_negative_beta_interval']['F']=='-12493/2338733088'
  minimal=next(c for c in stress['cases'] if c['label']=='minimal_four_atoms')
  assert F(minimal['mean']['F'])==F(117,8388608)
  dep=sorted((PROJECT/'runtime').rglob('*.py'))+[PROJECT/('rounds/'+p) for p in ['041/square_probe.py','042/row_probe.py','043/global_probe.py','044/tail_probe.py','046/prefix_probe.py','046/threshold_probe.py']]+sorted(HERE.glob('*.py'))
  kernel_keys=['positive_over_X','negative_over_X','record_overlap_removed_positive_over_X','density_difference_charged_positive_over_X','exclusive_second_source_charged_positive_over_X','exclusive_second_source_annular_pairs_checked']
  energy_keys=['Ef_over_X','Er_over_X','energy_sum_over_X','Bbar_over_X','Tobsbar_over_X','Gsym_over_X','imbalance_over_X','forward_record_collapse_checks','reverse_record_collapse_checks','weighted_symmetric_kernel_pairs_checked']
  receipt=dict(status='passed',round=47,date='2026-10-07',python=sys.version.split()[0],arithmetic='Fraction exact',
   stress=dict(batches=stress['batches'],exact_checks=stress['exact_checks'],
    cases=[dict(batch=c['batch'],label=c['label'],N=c['N'],D=c['D'],alpha=c['alpha'],mean_F_over_X=c['ratio']['F_over_X'],mean_Q_over_X=c['ratio']['Q_over_X'],beta_interval_count=c['beta_interval_count'],beta_negative_flow_probability=c['beta_negative_flow_probability']) for c in stress['cases']],
    negative_threshold_witness={k:negative[k] for k in ['label','atoms','alpha','X','beta_negative_flow_probability','first_negative_beta_interval']},
    full_spatial_checks=[dict(batch=c['batch'],label=c['label'],kernel={k:c['record_kernel'][k] for k in kernel_keys},energy={k:c['record_energy'][k] for k in energy_keys}) for c in full]),
   local_law=dict(case_count=local['case_count'],batch_counts=local['batch_counts'],rectangle_halfwidth=local['rectangle_halfwidth'],smoothing_radius=local['smoothing_radius'],continuous_lower=local['continuous_lower'],continuous_upper=local['continuous_upper'],
    parameter_ranges=[dict(batch=b,t_min=str(min(F(c['t']) for c in local['cases'] if c['batch']==b)),t_max=str(max(F(c['t']) for c in local['cases'] if c['batch']==b)),epsilon_min=str(min(F(c['epsilon']) for c in local['cases'] if c['batch']==b)),epsilon_max=str(max(F(c['epsilon']) for c in local['cases'] if c['batch']==b))) for b in range(3)],
    law_test_count=len(local['law_tests']),law_tests=local['law_tests']),
   proved_general_record_budget=True,main_theorem_proved=False,outer_flow_budget_proved=False,
   sha256={str(p.relative_to(PROJECT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in dep})
  print(json.dumps(receipt,ensure_ascii=False,indent=2))
if __name__=='__main__':main()
