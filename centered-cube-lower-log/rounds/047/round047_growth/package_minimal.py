"""Keep full raw summary local; emit exact compressed per-threshold numeric evidence."""
from pathlib import Path
import hashlib,json
import numpy as np
HERE=Path(__file__).resolve().parent
prior_index=json.loads((HERE/'PROFILE_INDEX.json').read_text()) if (HERE/'PROFILE_INDEX.json').exists() else {}
obj=json.loads((HERE/'RESULTS.json').read_text());rows=obj['rows']
offsets=[0];fields={k:[] for k in ['h','aW','aSE','aESS','aHits','bW','bSE','bESS','bHits']};configs=[]
for r in rows:
 a,b=r['seed_runs'];data=dict(h=r['h'],aW=a['W'],aSE=a['SE'],aESS=a['tail_ess'],aHits=a['proposal_tail_hits'],
                            bW=b['W'],bSE=b['SE'],bESS=b['tail_ess'],bHits=b['proposal_tail_hits'])
 for k in fields:fields[k].extend(data[k])
 offsets.append(offsets[-1]+len(r['h']));configs.append(dict(v=r['v'],R=r['R'],c=r['c'],seeds=r['seeds']))
np.savez_compressed(HERE/'PROFILE_MIN.npz',offsets=np.array(offsets,dtype=np.int64),
                    **{k:np.array(v,dtype=np.int32 if k.endswith('Hits') else np.float64) for k,v in fields.items()})
M=offsets[-1];epsall=(10/7)*np.sqrt(np.log(2*M/.01)/(2*48000));epshigh=(10/7)*np.sqrt(np.log(2*96/.01)/(2*48000))
index=dict(scope='Exact binary64 per-threshold summaries; full RESULTS and all v1 artifacts remain local, not recommended for Git.',
           samples_per_seed=24000,configurations=configs,profile_points=M,
           field_definition='Each configuration slices offsets[i]:offsets[i+1]. a/b are independent seeds; W mean, SE empirical standard error, ESS weighted-tail effective sample size, Hits proposal event hits. Combined mean=(aW+bW)/2 and SE=sqrt(aSE²+bSE²)/2. Fixed Gaussian curve and same-h ratio are reconstructed analytically.',
           simultaneous99_all_profile_Hoeffding_radius=float(epsall),simultaneous99_high96_Hoeffding_radius=float(epshigh),
           canonical_rows_sha256=hashlib.sha256(json.dumps(rows,sort_keys=True,separators=(',',':')).encode()).hexdigest(),
           local_full_result_file_sha256=hashlib.sha256((HERE/'RESULTS.json').read_bytes()).hexdigest(),
           local_v1_sha256={f:(hashlib.sha256((HERE/f).read_bytes()).hexdigest() if (HERE/f).exists()
                               else prior_index.get('local_v1_sha256',{}).get(f))
                            for f in ['FROZEN_PLAN.v1.json','run_growth.v1.py','RESULTS.v1.json','SUMMARY.v1.json']},
           history_elapsed_v1_seconds=(json.loads((HERE/'RESULTS.v1.json').read_text())['elapsed_seconds']
                                       if (HERE/'RESULTS.v1.json').exists() else prior_index.get('history_elapsed_v1_seconds')),
           v2_revision='v1 ran6.643s with two approximate tilts. v2 ran3.101s with exact size-biased t=sqrt(v), new disjoint base seeds47030001/47040001. v1 is not counted as independent replication of v2 or pooled.',
           minimum_archive=['run_growth.py','verify_and_report.py','package_minimal.py','FROZEN_PLAN.json','CHECKS.json','PROFILE_MIN.npz','PROFILE_INDEX.json','REPORT.md'])
index['PROFILE_MIN_sha256']=hashlib.sha256((HERE/'PROFILE_MIN.npz').read_bytes()).hexdigest()
(HERE/'PROFILE_INDEX.json').write_text(json.dumps(index,indent=2)+'\n')
print(json.dumps({'profile_points':M,'compressed_bytes':(HERE/'PROFILE_MIN.npz').stat().st_size,'all_profile_radius':epsall,'high96_radius':epshigh}))
