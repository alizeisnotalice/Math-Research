"""Independent covariance-density check; retain compressed full threshold summaries."""
from pathlib import Path
import json,hashlib
import numpy as np
from scipy.special import logsumexp
HERE=Path(__file__).resolve().parent
rng=np.random.default_rng(48990001);maxerr=0.;checks=0
for v in [16,32]:
 for c in [.2,.5]:
  R=3;J=5;t=np.arange(J)*.25;K=np.exp(-abs(t[:,None]-t[None,:])/v)
  C=np.zeros((1+R*J,1+R*J));C[0,0]=1
  for b in range(R):C[1+b*J:1+(b+1)*J,1+b*J:1+(b+1)*J]=K
  inv=np.linalg.inv(C);z=rng.standard_normal((40,1+R*J))@np.linalg.cholesky(C).T
  base=np.einsum('ni,ij,nj->n',z,inv,z);logs=[]
  for b in range(R):
   for anchor in range(J):
    mu=np.zeros(1+R*J);mu[0]=np.sqrt(c*v);mu[1+b*J:1+(b+1)*J]=np.sqrt((1-c)*v)*K[:,anchor]
    zz=z-mu;logs.append(-.5*np.einsum('ni,ij,nj->n',zz,inv,zz)+base/2)
  ld=logsumexp(np.array(logs),axis=0)-np.log(R*J)
  la=-v/2+np.sqrt(c*v)*z[:,0,None,None]+np.sqrt((1-c)*v)*z[:,1:].reshape(40,R,J)
  lb=logsumexp(la,axis=(1,2))-np.log(R*J)
  maxerr=max(maxerr,float(abs(ld-lb).max()));checks+=40
assert maxerr<1e-10
obj=json.loads((HERE/'RESULTS.json').read_text());rows=obj['rows'];summary=[]
for r in rows:
 baseline=next(x for x in rows if x['v']==r['v'] and x['c']==r['c'] and x['T']==r['T'] and x['R']==1)
 l=r['selected']['local'];n=r['selected']['node'];b=baseline['selected']['local']
 ratio=lambda x,y:dict(point=x['mean']/y['mean'],lower=x['lower']/y['upper'],
                       upper=x['upper']/y['lower'] if y['lower']>0 else None)
 summary.append({k:r[k] for k in ['v','c','R','T','J','seeds','fixed_single_node_best_h','fixed_single_node_best_W']}|
                dict(selected=r['selected'],local_over_R1_same_grid=ratio(l,b),local_over_own_node=ratio(l,n),
                     local_over_analytic_single=l['mean']/r['fixed_single_node_best_W']))
fields={};offsets=[0]
for i,r in enumerate(rows):
 for key in ['h']:fields.setdefault(key,[]).extend(r[key])
 for si,run in enumerate(r['runs']):
  for target in ['local','node']:
   for stat,values in run[target].items():fields.setdefault(f'{si}_{target}_{stat}',[]).extend(values)
 offsets.append(offsets[-1]+len(r['h']))
np.savez_compressed(HERE/'PROFILE_MIN.npz',offsets=np.array(offsets),
                    **{k:np.array(v,dtype=np.int32 if k.endswith('hits') else np.float64) for k,v in fields.items()})
controls=max(abs(run[k]['mean']-1)/run[k]['SE'] for r in rows for run in r['runs'] for k in ['RN_mass','mean_B'])
out=dict(scope=obj['scope'],elapsed_seconds=obj['elapsed_seconds'],independent_density_checks=checks,
         maximum_log_density_error=maxerr,max_control_abs_z=controls,
         confidence=obj['confidence'],confidence_source='https://www.cs.mcgill.ca/~colt2009/papers/012.pdf theorem4',
         profile_points=offsets[-1],rows=summary,
         source_sha256={name:hashlib.sha256((HERE/name).read_bytes()).hexdigest()
                        for name in ['FROZEN_PLAN.json','run_process.py','check_and_package.py']},
         canonical_numeric_rows_sha256=hashlib.sha256(json.dumps(rows,sort_keys=True,separators=(',',':')).encode()).hexdigest(),
         profile_npz_sha256=hashlib.sha256((HERE/'PROFILE_MIN.npz').read_bytes()).hexdigest(),
         profile_fields='Slices use offsets; 0/1 prefixes are independent seedA/seedB, local/node targets; variance is unbiased sample variance of the importance-weighted target, not variance of its mean. participation is weighted conditional softmax sum of squares at witnessed max/node0.',
         archive_files=['FROZEN_PLAN.json','run_process.py','check_and_package.py','SUMMARY.json','PROFILE_MIN.npz','REPORT.md'])
(HERE/'SUMMARY.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps({'density_checks':checks,'max_log_density_error':maxerr,'configs':len(rows),
                   'compressed_bytes':(HERE/'PROFILE_MIN.npz').stat().st_size,'elapsed_main':obj['elapsed_seconds']}))
