import json,sys,re,hashlib,itertools,math
from pathlib import Path
from decimal import Decimal
from collections import Counter
import numpy as np
p=Path(sys.argv[1]); out=Path(sys.argv[2])
c=json.load(open(p/'python_columns/columns71.json')); w=json.load(open(p/'wolfram/verify71.json')); g=json.load(open(p/'python_grid/grid71.json')); a=json.load(open(p/'point_formal/points71.json'))
def dec(s):
 s=re.sub(r'`[^*]*','',s).replace('*^','E');return Decimal(s)
def pairs(rows,keys,val):
 d={};errs=[]
 for r in rows:d.setdefault(tuple(r[k] for k in keys),[]).append(r)
 for v in d.values():
  assert len(v)==2;errs.append(abs(dec(v[0][val])-dec(v[1][val])))
 return len(errs),str(max(errs))
cp,ce=pairs(w['columns'],['N','regime'],'D');kp,ke=pairs(w['cells'],['N','d','M','mode'],'probability')
base={(r['N'],r['regime']):r for r in c['rows'] if r['refinement']==1}
cross=max(abs(float(dec(r['D']))-base[r['N'],r['regime']]['D']) for r in w['columns'] if r['precision']==100)
ref=[abs(r['D']-base[r['N'],r['regime']]['D']) for r in c['rows'] if r['refinement']!=1]
cellerr=[]
for r in w['cells']:
 if r['precision']!=100:continue
 row=next(x for x in g['rows'] if x['spec']['N']==r['N'] and x['spec']['d']==r['d'] and x['spec']['mode']==r['mode'])
 cellerr.append(abs(float(dec(r['probability']))-row['probabilities'][r['M']])/float(dec(r['probability'])))
res=[]
for r in g['rows']:
 s=r['spec'];n=s['N'];delta=s['step'];src=np.array(list(itertools.product(range(s['source_side']),repeat=n)));rec=np.array(list(itertools.product(range(-s['receiver_outer'],s['source_side']),repeat=n)));mass=np.array(r['source_mass']);sw=np.exp(-delta*src.sum(1))*(-math.expm1(-delta))**n
 worst=0
 for j in r['selected_indices']:
  diff=src-rec[j];ok=np.all(diff>=1,axis=1);q=np.array(r['probabilities'])[diff[ok].sum(1)];q=q*(2.**np.sum(diff[ok]>=s['d']+1,axis=1) if s['mode']=='finite' else 1);coef=s['budget']*q/sw[ok];coef[coef<=1e-14]=0
  worst=min(worst,float(coef@mass[ok]-1))
 assert mass.min()>=-1e-7 and mass.sum()<=1+1e-7 and worst>=-1e-7
 res.append(worst)
assert all(r['numeric'] and r['mode_certified'] for r in c['rows'])
assert all(r['numeric'] for r in w['columns']+w['cells'])
assert all(r['right_pass'] and r['left_pass'] and r['Python_rational_certificate'] for r in w['modes'])
assert hashlib.sha256((p/'python_columns/columns71.json').read_bytes()).hexdigest()==w['candidate_file_SHA256']
assert cross<2e-11 and max(cellerr)<2e-9 and max(ref)<2e-11 and Decimal(ce)<Decimal('1e-14') and Decimal(ke)<Decimal('1e-14')
r={'pass':True,'column_precision_pairs':cp,'column_precision_max_abs':ce,'cell_precision_pairs':kp,'cell_precision_max_abs':ke,'cross_column_count':32,'cross_column_max_abs':cross,'cross_cell_count':len(cellerr),'cross_cell_max_relative':max(cellerr),'Python_refinement_pairs':len(ref),'Python_refinement_max_abs':max(ref),'rational_mode_certificates':len(c['rows']),'Wolfram_mode_checks':len(w['modes']),'grid_primal_checks':len(res),'grid_worst_threshold_residual':min(res),'grid_statuses':dict(Counter(x['status'] for x in g['rows'])),'grid_max_ratio':max(x['collocation_ratio'] for x in g['rows']),'grid_max_solver_gap':max(x['solver_gap'] for x in g['rows']),'atomic_max_ratio':max(x['ratio'] for x in a['rows']),'atomic_evaluations':sum(x['evaluations'] for x in a['rows']),'atomic_invalid_evaluations':sum(x['invalid_objective_evaluations'] for x in a['rows']),'atomic_converged':sum(x['solver_success'] for x in a['rows']),'scope':'Floating consistency gates; only stated posterior mode/logarithm brackets and geometric lemmas are rigorous certificates.'}
out.write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(r,indent=2))
