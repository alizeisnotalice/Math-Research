#!/usr/bin/env python3
"""Run final round51 probes; retain only the necessary verification summary."""
from pathlib import Path
from decimal import Decimal as D
import hashlib,json,os,subprocess,sys,tempfile
HERE=Path(__file__).resolve().parent

def main():
 env=dict(os.environ,PYTHONDONTWRITEBYTECODE='1')
 with tempfile.TemporaryDirectory(prefix='euclidean51_') as directory:
  output=Path(directory)/'geometry.json'
  subprocess.run([sys.executable,str(HERE/'geometry_probe.py'),'--output',str(output)],check=True,capture_output=True,text=True,env=env)
  geometry=json.loads(output.read_text())
  def run(name):return json.loads(subprocess.run([sys.executable,str(HERE/name)],check=True,capture_output=True,text=True,env=env).stdout)
  budget=run('budget_probe.py');radial=run('radial_probe.py')
  assert geometry['status']==budget['status']==radial['status']=='passed'
  assert geometry['batch_counts']==[3,3,3] and len(geometry['whole_ball'])==9
  assert geometry['exact_ray_checks']==2346 and len(geometry['mask_obstruction'])==9
  assert budget['batch_counts']==[5,4,4] and len(budget['cases'])==13
  assert budget['cases'][-1]['n']==4096 and budget['cases'][-1]['L']==18
  assert radial['clock_checks']==96 and radial['D']==7 and radial['half_threshold_labels']==[6,7]
  keys=['batch','n','outer_probability','conditional_depth_mean','conditional_repetition_partial','repetition_partial_terms','cross_quadrature_absolute_error','interval_certificate']
  receipt=dict(status='passed',round=51,date='2026-10-07',python=sys.version.split()[0],budget=budget,
   geometry=dict(batch_counts=geometry['batch_counts'],exact_ray_checks=2346,arbitrary_mask_certificates=geometry['mask_obstruction'],whole_ball=[{key:c[key] for key in keys} for c in geometry['whole_ball']],arithmetic='120-digit Decimal, 64/96-node quadrature comparison; no interval certification',max_quadrature_difference=str(max(D(c['cross_quadrature_absolute_error']) for c in geometry['whole_ball']))),
   actual_radial_obstruction=radial,main_theorem_proved=False,dimension_free_outer_budget_proved=False,
   sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(HERE.glob('*.py'))})
  print(json.dumps(receipt,ensure_ascii=False,indent=2))
if __name__=='__main__':main()
