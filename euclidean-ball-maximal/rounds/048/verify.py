#!/usr/bin/env python3
"""Recompute round48, keeping only a compact receipt; standard library only."""
from pathlib import Path
import hashlib,json,os,subprocess,sys,tempfile
from fractions import Fraction as F
sys.dont_write_bytecode=True
HERE=Path(__file__).resolve().parent;PROJECT=HERE.parents[1]
def main():
 with tempfile.TemporaryDirectory(prefix='euclidean48_') as folder:
  env=dict(os.environ,PYTHONDONTWRITEBYTECODE='1');out=Path(folder)/'shear.json'
  subprocess.run([sys.executable,str(HERE/'shear_probe.py'),'--output',str(out)],check=True,capture_output=True,text=True,env=env)
  shear=json.loads(out.read_text())
  geo=json.loads(subprocess.run([sys.executable,str(HERE/'geometry_probe.py')],check=True,capture_output=True,text=True,env=env).stdout)
  assert shear['status']==geo['status']=='passed' and shear['batch_counts']==[3,2,4]
  assert geo['exact_n1_area_checks']==6739 and geo['scalar_cap_checks']==2535
  assert len(geo['high_dimension'])==len(geo['relaxed_box'])==9
  assert sum(c['positive_outer_source_polygons'] for c in shear['cases'])==493
  assert sum(c['target_classification_pieces'] for c in shear['cases'])==1335
  assert sum(c['beta_target_record_pairs'] for c in shear['cases'])==462
  assert sum(c.get('archived_outer_kernel_crosscheck')=='passed' for c in shear['cases'])==8
  new=next(c for c in shear['cases'] if c['label']=='minimal_cloud_with_two_coarse_atoms')
  assert new['class_fraction']['band_J1']=='55/104' and new['remainder_over_O']=='95/104'
  for c in shear['cases']:
   assert c['exact_error']=='0' and F(c['pushforward_density_integral'])==F(c['O'])
   if c['label'] in ['chain_L4_alpha2','chain_L8_alpha2','chain_L16_alpha2']:
    assert c['class_fraction']['MP_le_alpha']=='1'
  dep=sorted((PROJECT/'runtime').rglob('*.py'))+[PROJECT/('rounds/'+p) for p in ['041/square_probe.py','042/row_probe.py','043/global_probe.py','044/tail_probe.py','046/prefix_probe.py','046/threshold_probe.py','047/flow_probe.py','047/verification.json']]+sorted(HERE.glob('*.py'))
  keys=['batch','label','N','alpha','O_over_X','class_fraction','compatible_over_O','remainder_over_O','pushforward_density_max_over_alpha','positive_outer_source_polygons','target_classification_pieces','beta_target_record_pairs','exact_error']
  receipt=dict(status='passed',round=48,date='2026-10-07',python=sys.version.split()[0],
   shear=dict(batch_counts=shear['batch_counts'],cases=[{k:c[k] for k in keys} for c in shear['cases']],independent_archived_crosschecks=8,new_J1_certificate={k:new[k] for k in ['atoms','alpha','O','escape_witness_by_kind']}),
   geometry=dict(exact_n1_area_checks=geo['exact_n1_area_checks'],n1_batches=[{k:b[k] for k in ['batch','max_gap','source_distance_denominator','exact_area_checks']} for b in geo['n1']],high_dimension=geo['high_dimension'],scalar_cap_checks=geo['scalar_cap_checks'],relaxed_box_parameters=[{k:c[k] for k in ['batch','j','k','per_scale_lower']} for c in geo['relaxed_box']],high_dimension_scope=geo['high_dimension_scope'],analytic_all_dimension_sum_bound=83,analytic_rational_upper='497/6',actual_outer_cost='0',relaxed_outer_lower=geo['relaxed_outer_lower']),
   main_theorem_proved=False,actual_outer_budget_proved=False,
   sha256={str(p.relative_to(PROJECT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in dep})
  print(json.dumps(receipt,ensure_ascii=False,indent=2))
if __name__=='__main__':main()
