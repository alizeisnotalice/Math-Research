#!/usr/bin/env python3
"""Exact allocation certificates, chain budgets and original outer columns."""
from pathlib import Path
import json,subprocess,sys,os,tempfile,hashlib
from fractions import Fraction as F
HERE=Path(__file__).resolve().parent;PROJECT=HERE.parents[1]
def main():
 env=dict(os.environ,PYTHONDONTWRITEBYTECODE='1')
 with tempfile.TemporaryDirectory(prefix='euclidean56_') as folder:
  results={}
  for name in ('allocation','column'):
   out=Path(folder)/(name+'.json')
   subprocess.run([sys.executable,str(HERE/(name+'_probe.py')),'--repo-root',str(PROJECT),'--output',str(out)],check=True,capture_output=True,text=True,env=env)
   results[name]=json.loads(out.read_text());assert results[name]['status']=='passed'
  alloc,cols=results['allocation'],results['column']
  assert alloc['batch_counts']==[7,3,6] and cols['batch_counts']==[8,4,8]
  assert len(alloc['cases'])==16 and len(cols['cases'])==20
  assert cols['total_positive_columns']==428 and cols['total_active_edges']==2765 and cols['total_breakpoint_checks']==3662
  bylabel={c['label']:c for c in cols['cases']}
  for a in alloc['cases']:
   c=bylabel[a['label']];assert F(a['M'])==F(c['fine_record_mass_M']) and F(a['X'])==F(c['X'])
  acases=[{k:v for k,v in c.items() if k not in ('primal_rows','chain_cover','history')} for c in alloc['cases']]
  ckeys=['batch','label','N','D','alpha','X','O_over_X','M_over_X','O_over_M_weighted_mean_conditional_column','Jaccard_O_over_X','tau_O_over_X','maximum_raw_column','maximum_conditional_column','positive_columns','active_edges','breakpoint_checks']
  ccases=[{k:c[k] for k in ckeys} for c in cols['cases']]
  deps=sorted((PROJECT/'runtime').rglob('*.py'))+[PROJECT/'rounds'/p for p in ['041/square_probe.py','042/row_probe.py','043/global_probe.py','044/tail_probe.py','046/prefix_probe.py','046/threshold_probe.py','047/flow_probe.py','049/global_probe.py','053/obstruction_probe.py','054/lag_probe.py','055/gram_certificate.py','055/verification.json']]+sorted(HERE.glob('*.py'))
  receipt=dict(status='passed',round=56,date='2026-10-08',python=sys.version.split()[0],arithmetic='Fraction exact; no floating point feasibility decisions',allocation=dict(batch_counts=alloc['batch_counts'],cases=acases,primal_and_cut_dual_certificates_checked=16,exhaustive_cut_crosschecks=sum(c['exhaustive_crosscheck'] for c in alloc['cases']),chain_tail_constant='3/2'),columns=dict(batch_counts=cols['batch_counts'],cases=ccases,positive_columns=428,active_edges=2765,breakpoint_checks=3662,previous_O_receipts=16,independent_new_O_checks=4,Jaccard_budget_constant='3/8',tau_budget_constant='3/4'),scope='Finite n=1 full-E diagnostics; independent analytic chain lemma applies in arbitrary n with its stated structure. No general chain cover, column bound or source-capacity bound proved.',main_weak_type_theorem_proved=False,sha256={str(p.relative_to(PROJECT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in deps})
  print(json.dumps(receipt,ensure_ascii=False,indent=2))
if __name__=='__main__':main()
