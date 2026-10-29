"""Independent finite-density checks and analysis of the frozen Gaussian experiment."""
from pathlib import Path
import hashlib,json
import numpy as np
from scipy.special import logsumexp
HERE=Path(__file__).resolve().parent
prior_index=json.loads((HERE/'PROFILE_INDEX.json').read_text()) if (HERE/'PROFILE_INDEX.json').exists() else {}
obj=json.loads((HERE/'RESULTS.json').read_text());rows=obj['rows'];rng=np.random.default_rng(47990001)
maxerr=0.;checks=0
for v in [8,64]:
 for R in [1,4]:
  for c in [0,.25,.5,.9,1]:
   z=rng.standard_normal((20,R+1));base=np.sum(z*z,axis=1)
   logs=[]
   for anchor in range(R):
    mu=np.zeros(R+1);mu[0]=np.sqrt(v*c);mu[anchor+1]=np.sqrt(v*(1-c))
    logs.append(-.5*np.sum((z-mu)**2,axis=1)+base/2)
   density=logsumexp(np.array(logs),axis=0)-np.log(R)
   y=np.sqrt(c)*z[:,0,None]+np.sqrt(1-c)*z[:,1:]
   logS=logsumexp(-v/2+np.sqrt(v)*y,axis=1)-np.log(R)
   maxerr=max(maxerr,float(abs(density-logS).max()));checks+=len(z)
assert maxerr<1e-12
alpha=.01;cap=10/7;N=24000
# Seed1 selects within each profile. Seed2 is independent; union only over 96 selected thresholds.
eps=cap*np.sqrt(np.log(2*len(rows)/alpha)/(2*N))
selected=[]
for r in rows:
 W=r['independent_W_at_discovery_h'];fixed=r['fixed_global_best_W']
 selected.append(dict(v=r['v'],R=r['R'],c=r['c'],h=r['discovery_h'],W=W,
                      se=r['independent_SE_at_discovery_h'],gain=W/fixed,
                      simultaneous99_lower_gain=max(0,W-eps)/fixed,
                      simultaneous99_upper_gain=(W+eps)/fixed,
                      high_ratio=r['high_baseline']['ratio_to_fixed']))
controls=[];fixedz=[];zero=0
for r in rows:
 for sample in r['seed_runs']:
  for name in ['RN_mass','mean_one_control']:
   x=sample[name];controls.append(abs(x['mean']-1)/x['se'])
  zero+=sum(n==0 for n in sample['proposal_tail_hits'])
 if r['R']==1:
  exact=np.array(r['fixed_W']);se=np.array(r['SE']);estimate=np.array(r['W'])
  fixedz.extend(abs(estimate-exact)/se)
assert zero==0
summ=dict(scope='Gaussian proxy only; exact density algebra checked numerically, Hoeffding cap derived algebraically; no HL/C_n inference.',
         independent_density_checks=checks,max_absolute_log_density_error=maxerr,
         max_control_absolute_z=float(max(controls)),max_R1_profile_absolute_z=float(max(fixedz)),
         proposal_zero_hit_profile_entries=zero,independent_test_simultaneous_confidence=.99,
         independent_test_Hoeffding_radius=float(eps),deterministic_target_cap=cap,
         independent_selected_thresholds=selected,
         best_by_v=[max([r for r in selected if r['v']==v],key=lambda r:r['gain']) for v in [8,16,32,64]],
         elapsed_main_v2_seconds=obj['elapsed_seconds'],
         elapsed_v1_seconds=(json.loads((HERE/'RESULTS.v1.json').read_text())['elapsed_seconds']
                             if (HERE/'RESULTS.v1.json').exists() else prior_index.get('history_elapsed_v1_seconds')),
         input_sha256={f:hashlib.sha256((HERE/f).read_bytes()).hexdigest() for f in ['FROZEN_PLAN.json','run_growth.py','RESULTS.json']})
(HERE/'CHECKS.json').write_text(json.dumps(summ,indent=2)+'\n')
print(json.dumps({k:summ[k] for k in ['max_absolute_log_density_error','max_control_absolute_z','max_R1_profile_absolute_z','independent_test_Hoeffding_radius','best_by_v']},indent=2))
