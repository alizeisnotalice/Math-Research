"""Read replay logs/invariants and a synthetic underflow example; no MC."""
from pathlib import Path
import json,hashlib,math
import numpy as np

HERE=Path(__file__).resolve().parent;SRC=HERE.parent/'cell'
old=json.loads((SRC/'RESULTS.json').read_text());new=json.loads((SRC/'RESULTS_stable_ESS.json').read_text())
rev=json.loads((SRC/'ESS_revision.json').read_text())
mapping={'probe.py':'original_source_SHA256','registration.json':'registration_SHA256',
         'RESULTS.json':'original_RESULTS_SHA256','replay_stable_ess.py':'fix_script_SHA256',
         'ESS_FIX.patch':'minimal_patch_SHA256'}
hashes={n:hashlib.sha256((SRC/n).read_bytes()).hexdigest()==rev[k]for n,k in mapping.items()}
assert all(hashes.values())
keys=['metadata','N','seed','mass_weight_mean','mass_SE_descriptive','weight_max','uniform_rows',
      'source_b_counts_by_model','max_actual_coefficient_contribution','contribution_cap_all_thresholds',
      'max_PB_probability_mass_residual','positive_subnormal_probability_entries','subnormal_loss_magnitude_budget_only']
same=True;hold_same=True;raw=[]
for a,b,check in zip(old['cases'],new['cases'],new['same_data_protocol_replay_checks']):
    assert all(a[k]==b[k]for k in keys)
    hold_same &= all(a['models'][m]['holdout']==b['models'][m]['holdout']for m in a['models'])
    path=SRC/check['local_raw_logdata']
    assert hashlib.sha256(path.read_bytes()).hexdigest()==check['raw_logdata_SHA256']
    data=np.load(path);lm=data['logM'];ld=data['logD'];N=a['N']
    assert lm.shape==(N,a['metadata']['R']+2)and ld.shape==(N,)
    stream=hashlib.sha256()
    for start in range(0,N,64):
        stream.update(np.ascontiguousarray(lm[start:start+64]).tobytes())
        stream.update(np.ascontiguousarray(ld[start:start+64]).tobytes())
    assert stream.hexdigest()==check['logM_logD_batch_stream_SHA256']
    bad={}
    for i,name in enumerate(a['models']):
        order=np.argsort(lm[:,i])[::-1];w=np.exp(-ld[order])
        with np.errstate(invalid='ignore',divide='ignore',under='ignore'):
            ess=np.cumsum(w)**2/np.cumsum(w*w)
        bad[name]=int(np.count_nonzero(~np.isfinite(ess)))
    raw.append({'n':a['metadata']['n'],'R':a['metadata']['R'],
                'raw_hash_and_stream_passed':True,'old_formula_nonfinite_prefixes':bad})
assert new['macroscopic_diagnostics']==old['macroscopic_diagnostics']
lw=np.array([-400.,-401.,-402.]);w=np.exp(lw)
scaled=np.exp(lw-lw.max())
expected=math.fsum(scaled)**2/math.fsum(x*x for x in scaled)
stable=float(np.exp(2*np.logaddexp.accumulate(lw)[-1]-np.logaddexp.accumulate(2*lw)[-1]))
assert abs(stable-expected)<1e-12
out={'scope':'No sampling: saved replay-log/invariant checks plus deterministic underflow arithmetic.',
     'hashes':hashes,'nine_case_invariants_equal':True,'all_120_holdout_CI_entries_unchanged':bool(hold_same),
     'original_main_seconds':old['runtime_seconds'],'same_protocol_replay_seconds':new['runtime_seconds'],
     'two_execution_seconds_total':old['runtime_seconds']+new['runtime_seconds'],
     'original_raw_logs_missing':not rev['original_raw_logM_logD_saved'],
     'raw_checks':raw,'synthetic_scaled_ESS':expected,'synthetic_log_ESS':stable,
     'caution':'Matching protocol/invariants plus new replay hashes do not constitute an original-vs-replay full raw hash comparison; original raw logs were not saved.'}
(HERE/'ess_audit_checks.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps({k:v for k,v in out.items()if k!='raw_checks'},indent=2))
print(json.dumps(raw[-1],indent=2))
