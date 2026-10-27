#!/usr/bin/env python3
"""R063 final identity/data lock only; no threshold calculation repeated."""
import hashlib
import json
from fractions import Fraction as Q
from pathlib import Path

HERE=Path(__file__).resolve().parent
BASE=Path('/Users/zhengzhihao/Desktop/T')
SOURCE=BASE/'cube-r030-parameters/cube-dim-7-bounds/upper/round063/wideband_threshold/THRESHOLD_PROOF.md'
EXPECTED='70fc76d864f233900baadd7f177b7532b654355c3ce3b515a822fe6b65a19133'
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
assert sha(SOURCE)==EXPECTED
receipt=json.loads((HERE/'compressed_receipt.json').read_text())
assert receipt['status']=='PASS_INDEPENDENT_COMPRESSED_R063_THRESHOLD'
assert receipt['own_script_sha256']==sha(HERE/'check_compressed_threshold.py')
assert receipt['sources']['R063_source']['sha256']==EXPECTED
assert all(sha(Path(v['path']))==v['sha256'] for v in receipt['sources'].values())

# Exact accepted-data identity, without rerunning the old R056 budget or verifier.
G56=Q(172531397335235056462877787662549380239,1643616601687301779182642153860479582208)
old56=Path(receipt['sources']['R056_historical_values']['path'])
assert Q(json.loads(old56.read_text())['gamma'])==G56
assert receipt['base_values']['gamma56_accepted_readonly']==str(G56)
assert str(G56.numerator)+'/' in SOURCE.read_text()
assert str(G56.denominator) in SOURCE.read_text()

files=[SOURCE,HERE/'review.md',HERE/'check_compressed_threshold.py',HERE/'compressed_receipt.json',Path(__file__).resolve()]
out={'status':'R063_FULL_THRESHOLD_ACCEPT_SOURCE_LOCK',
 'author_source_sha256':EXPECTED,'audit_conclusion':'FULL ACCEPT for every integer n with natural log(n)>=2^33',
 'new_threshold_round':'R063 only; R059 historical eventual record unchanged',
 'accepted_gamma56_identity_matches':True,'new_mathematical_run_in_this_lock':False,
 'sha256':{str(p):sha(p) for p in files},
 'preserved_dependencies':receipt['sources'],
 'scope':'Identity/data check after the new independent compressed threshold computation. No author program, accepted verifier, old log/grid, MATLAB or Git operation was executed.'}
(HERE/'source_lock.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps({'status':out['status'],'author_source_sha256':EXPECTED,'review_sha256':sha(HERE/'review.md')},indent=2))
