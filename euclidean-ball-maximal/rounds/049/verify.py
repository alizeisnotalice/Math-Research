#!/usr/bin/env python3
"""Recompute global records and continuous-background obstruction; standard library."""
from pathlib import Path
from fractions import Fraction as F
import hashlib,json,os,subprocess,sys,tempfile
sys.dont_write_bytecode=True
HERE=Path(__file__).resolve().parent;PROJECT=HERE.parents[1]
def main():
 with tempfile.TemporaryDirectory(prefix='euclidean49_') as folder:
  env=dict(os.environ,PYTHONDONTWRITEBYTECODE='1');out=Path(folder)/'global.json'
  subprocess.run([sys.executable,str(HERE/'global_probe.py'),'--output',str(out)],check=True,capture_output=True,text=True,env=env)
  glob=json.loads(out.read_text())
  counter=json.loads(subprocess.run([sys.executable,str(HERE/'counterexample_probe.py')],check=True,capture_output=True,text=True,env=env).stdout)
  assert glob['status']==counter['status']=='passed' and glob['batch_counts']==[7,6,9]
  positive=[c for c in glob['cases'] if F(c['O'])>0];assert len(positive)==13
  assert all(F(c['new_paid'])==F(c['O'])==F(c['same_label'])<=F(c['Bglobal']) for c in glob['cases'])
  assert sum(c.get('round48_crosscheck')=='passed' for c in glob['cases'])==9
  assert sum(c.get('independent_eligible_flow_B_check')=='passed' for c in glob['cases'])==9
  assert counter['batch_counts']==[3,3,3] and counter['exact_global_label_checks']==7155
  assert counter['smooth_extreme_shift_checks']==2120 and counter['activation_scalar_checks']==81
  maximum=max(glob['cases'],key=lambda c:F(c['Bglobal_over_X']))
  assert maximum['label']=='geometric_contact_clouds_s491107'
  keys=['batch','label','N','D','alpha','O_over_X','Bglobal_over_X','Beligible_over_X','new_paid_over_O','caught_original_escape_fraction','remainder_over_O','global_record_relation','exact_error']
  dep=sorted((PROJECT/'runtime').rglob('*.py'))+[PROJECT/('rounds/'+p) for p in ['041/square_probe.py','042/row_probe.py','043/global_probe.py','044/tail_probe.py','046/prefix_probe.py','046/threshold_probe.py','047/flow_probe.py','048/shear_probe.py','048/verification.json']]+sorted(HERE.glob('*.py'))
  receipt=dict(status='passed',round=49,date='2026-10-07',python=sys.version.split()[0],arithmetic='Fraction exact; smooth family checked by worst-shift rational bounds',
   global_records=dict(case_count=22,batch_counts=glob['batch_counts'],positive_case_count=13,round48_crosschecks=9,independent_eligible_B_crosschecks=9,counts={k:sum(c[k] for c in glob['cases']) for k in ['comparison_pairs','source_polygons','target_pieces','target_record_pairs']},cases=[{k:c[k] for k in keys} for c in glob['cases']]),
   counterexample=counter,main_theorem_proved=False,retained_activation_bridge_bound_proved=False,
   sha256={str(p.relative_to(PROJECT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in dep})
  print(json.dumps(receipt,ensure_ascii=False,indent=2))
if __name__=='__main__':main()
