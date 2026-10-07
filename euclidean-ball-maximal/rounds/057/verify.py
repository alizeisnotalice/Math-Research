#!/usr/bin/env python3
"""Reproduce exact capture widths, subtree capacities and smooth-grid constants."""
from pathlib import Path
from fractions import Fraction as F
import hashlib,json,os,subprocess,sys,tempfile
HERE=Path(__file__).resolve().parent;PROJECT=HERE.parents[1]
def main():
 env=dict(os.environ,PYTHONDONTWRITEBYTECODE='1')
 with tempfile.TemporaryDirectory(prefix='euclidean57_') as folder:
  out=Path(folder)/'poset.json'
  subprocess.run([sys.executable,str(HERE/'poset_probe.py'),'--output',str(out)],check=True,capture_output=True,text=True,env=env)
  data=json.loads(out.read_text())
  grid=json.loads(subprocess.run([sys.executable,str(HERE/'grid_certificate.py')],check=True,capture_output=True,text=True,env=env).stdout)
  assert data['status']==grid['status']=='passed'
  assert data['batch_counts']==[7,3,10] and len(data['cases'])==20
  assert data['source_width_certificates_checked']==176 and data['original_receipts_checked']==16
  assert sum(c['laminar'] for c in data['cases'])==9
  assert len(grid['cases'])==9 and sum(c['actual_positive_cells'] for c in grid['cases'])==1031
  abstract=[]
  for batch,depths in enumerate(((1,2,4),(8,16,32),(64,128,256))):
   for L in depths:
    # Exact algebraic tree model, not Euclidean capture data.
    path_ratio=F(3,2)*(1-F(1,2**(L+1)))
    assert path_ratio==F(3,4)*sum((F(1,2**k) for k in range(L+1)),F(0))<F(3,2)
    root_ratio=F(3,4)*(L+1)
    assert root_ratio>=path_ratio
    abstract.append(dict(batch=batch,depth=L,chain_tail_ratio=str(path_ratio),root_capacity=str(root_ratio)))
  keys=['batch','label','N','D','capture_classes','X','M','C_alloc','natural_max','atom_width_lower_bound','greedy_kappa','width_over_C','strongest_source','explicit_antichain','tight_source_set','tight_source_mass','tight_cut_demand','ratio_iteration_count','laminar','subtree_C','crossing_witness']
  compact=[{k:c[k] for k in keys} for c in data['cases']]
  names=['041/square_probe.py','042/row_probe.py','043/global_probe.py','044/tail_probe.py','046/prefix_probe.py','046/threshold_probe.py','047/flow_probe.py','053/obstruction_probe.py','054/lag_probe.py','056/allocation_probe.py','056/verification.json']
  deps=sorted((PROJECT/'runtime').rglob('*.py'))+[PROJECT/'rounds'/p for p in names]+sorted(HERE.glob('*.py'))
  receipt=dict(status='passed',round=57,date='2026-10-08',python=sys.version.split()[0],poset=dict(batch_counts=data['batch_counts'],source_width_certificates=176,standalone_checks=4,previous_receipts=16,laminar_full_capture_cases=9,cases=compact),grid=grid,abstract_tree=dict(scope='Symbolic finite tree model only; no Euclidean realization claimed',cases=abstract),scope='Full-E finite n=1 exact integration and allocation; separate analytic all-h smooth counterexample to unweighted chain cover. No general dimension-free weak bound proved.',main_weak_type_theorem_proved=False,sha256={str(p.relative_to(PROJECT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in deps})
  print(json.dumps(receipt,ensure_ascii=False,indent=2))
if __name__=='__main__':main()
