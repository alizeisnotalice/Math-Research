"""Read saved integer endpoints only; no native calls or writes."""
import sys
sys.dont_write_bytecode=True
from probe import H,sha,check,NS,DEN,T
from pathlib import Path
import json
r=json.loads((H/'registration.json').read_text())
e=json.loads((H/'execution.json').read_text())
assert sha(H/'probe.py')==r['driver_sha256']
assert sha(H/'registration.json')==e['registration_sha256']
for path,h in r['source_sha256'].items():assert sha(Path(path))==h
assert r['numerators']==NS and r['denominator']==DEN and r['T']==str(T)
assert e['returncode']==0 and e['all_children_waited'] and e['elapsed_seconds']<120
assert sha(H/'native.stdout.json')==e['stdout_sha256']
assert sha(H/'native.stderr.txt')==e['stderr_sha256']
assert not(H/'native.stderr.txt').read_text()
out=check(json.loads((H/'native.stdout.json').read_text()))
assert out==json.loads((H/'checked_results.json').read_text())
assert out['strictly_excludes64'] and not out['necessary64_passes']
print('PASS saved-only: new62-frame, exact endpoints/Euler/1891pairS2;0native')
