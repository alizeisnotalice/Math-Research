"""Read-only identity/shape/chain audit of accepted protocol replay, no RNG calls."""
from pathlib import Path
import json,hashlib,math
import numpy as np
HERE=Path(__file__).resolve().parent;REP=HERE.parent/'round041_replay';OLD=HERE.parent/'round040_stream'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
receipt=json.loads((REP/'RECONSTRUCTION_RECEIPT.json').read_text());plan=json.loads((REP/'PLAN.json').read_text())
assert receipt['plan_SHA256']==sha(REP/'PLAN.json')and receipt['all_acceptance_checks_passed']
assert receipt['new_independent_samples']==0
for field,path in [('old_plan_SHA256',OLD/'PLAN.json'),('base_driver_SHA256',OLD/'probe.py'),
 ('observer_SHA256',REP/'reconstruct.py'),('upper_math_SHA256',REP/'upper_math.py'),
 ('postprocessing_plan_SHA256',HERE.parent/'round041_post/PLAN.json'),
 ('MCD_proof_SHA256',HERE.parent/'round041_mcd/PROOF.md')]:assert plan[field]==sha(path)
rows=[];arrays=0
for rec in receipt['cases']:
    R=rec['R'];old=json.loads((OLD/f'R{R}_RESULTS.json').read_text());re=json.loads((REP/f'base_replay/R{R}_RESULTS.json').read_text())
    for key in ['count_stream_SHA256','sampled_anchor_nominal_stream_SHA256','stream_statistics_SHA256']:
        assert re[key]==old[key]
    with np.load(OLD/f'R{R}_paired.npz',allow_pickle=False)as a,np.load(REP/f'base_replay/R{R}_paired.npz',allow_pickle=False)as b:
        assert set(a.files)==set(b.files)
        for key in a.files:assert np.array_equal(a[key],b[key]);arrays+=1
        oldupper=a['upper'];oldlower=a['lower'];h=a['h'];logD=a['logD']
    assert rec['upper_inputs_SHA256']==sha(REP/rec['upper_inputs_file'])
    with np.load(REP/rec['upper_inputs_file'],allow_pickle=False)as p:
        assert np.array_equal(p['h'],h)and np.array_equal(p['logD'],logD)
        assert p['logmean_upper'].shape==p['logvariance_old'].shape==(3,256,233)
        assert p['max_table_logproduct'].shape==p['table0_logproduct'].shape==(256,233)
        assert np.all(p['max_table_logproduct']>=p['table0_logproduct'])and np.all(p['table0_logproduct']>=0)
        upper=p['new_upper'];var=p['variance_only_upper']
        B=oldlower*np.exp(logD[None,:,None]-h[None,None,:])
        lower=3*oldlower/(1+2*B)
        assert upper.shape==var.shape==oldupper.shape
        assert np.all(np.isfinite(upper))and np.all(upper>=0)
        assert np.max(lower-upper)<1e-9 and np.max(upper-var)<1e-9 and np.max(var-oldupper)<1e-9
        maxratio=0.
        for mi,model in enumerate(plan['targets']['models']):
            for hi,hr in enumerate(h):
                cap=old['models'][model]['profile'][hi]['cap']
                maxratio=max(maxratio,float(upper[mi,:,hi].max()/cap));assert upper[mi,:,hi].max()<=cap*(1+1e-8)
        rows.append({'R':R,'three_stream_hashes_equal':True,'all_old_payload_arrays_equal':True,
                     'upper_inputs_hash_matches':True,'max_checked_source_r':rec['max_checked_source_r'],
                     'max_sharp_lower_minus_new_upper':float(np.max(lower-upper)),
                     'max_new_upper_minus_variance_only':float(np.max(upper-var)),
                     'max_variance_only_minus_old_upper':float(np.max(var-oldupper)),
                     'new_upper_strict_improvement_over_variance_only_count':int((upper<var-1e-10).sum()),
                     'max_new_upper_over_old_cap':maxratio})
out={'cases':rows,'old_payload_arrays_compared':arrays,'all_replay_identity_checks_pass':True,
     'replay_seconds':receipt['total_seconds'],'new_independent_samples':0,
     'scope':'No RNG/source recomputation; verifies stored replay identity, bounds chain and caps. Offline per-scale new-upper recomputation assigned to separate reviewer.'}
target=HERE/'replay_checks.json'
if target.exists():raise RuntimeError('Refusing replacement')
target.write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
