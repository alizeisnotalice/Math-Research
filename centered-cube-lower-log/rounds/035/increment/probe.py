"""Frozen local-increment diagnostic; exact finite-atom pair sum, ordinary float."""
from pathlib import Path
import json,hashlib,math
import numpy as np
out=Path(__file__).resolve().parent
plan={'seed':202610290351,'K':[63,255,1023,4095], 'replicates':16,
      'tau':np.geomspace(1e-5,1,81).tolist(),'scope':'finite source kernel increments, not maximal lower bound'}
reg=out/'registration.json'
if reg.exists(): raise RuntimeError('Refusing repeated registered run')
plan['source_sha256']=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
reg.write_text(json.dumps(plan,indent=2))
rng=np.random.default_rng(plan['seed']);rows=[]
for K in plan['K']:
  for rep in range(plan['replicates']):
    u=rng.random(K);gap=1+np.roll(u,-1)-u
    tau=np.array(plan['tau']); pair=np.maximum(0,tau[:,None]-gap).sum(axis=1)/K
    decrement=tau/2-tau*tau/2+pair
    expected=tau/2-tau*tau/2+tau**3/6
    relative=(decrement-expected)/(tau/2)
    rows.append({'K':K,'rep':rep,'min_scaled_gap':float(gap.min()),
      'maximum_relative_increment_error':float(abs(relative).max()),
      'error_tau_le_001':float(abs(relative[tau<=.01]).max()),
      'at_tau_1':float(relative[-1]),
      'normalized_decrement':(decrement/(tau/2)).tolist()})
res={'numpy':np.__version__,'registration_sha256':hashlib.sha256(reg.read_bytes()).hexdigest(),'rows':rows}
(out/'results.json').write_text(json.dumps(res,indent=2))
for K in plan['K']:
 a=[r for r in rows if r['K']==K]
 print(K,'max_all',max(r['maximum_relative_increment_error'] for r in a),
       'max_local',max(r['error_tau_le_001'] for r in a))
