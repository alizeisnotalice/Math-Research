from pathlib import Path
from fractions import Fraction as F
import hashlib,json,gzip
p=Path(__file__).resolve().parent;H=lambda b:hashlib.sha256(b).hexdigest();read=lambda x:json.loads((p/x).read_text());m=read('archive_manifest.json')
for x,h in m['sha256'].items():assert H((p/x).read_bytes())==h
c=read('candidate_certificate.json');old=read('candidate_initial_sign_error.json');assert c['code_sha256']==H((p/'candidate.py').read_bytes()) and old['code_sha256']==H((p/'candidate_initial_sign_error.py').read_bytes());assert all(c['verdict'].values()) and not old['verdict']['tau_identity'];assert F(c['exact']['tau_gap'])==F(1,20000) and F(c['exact']['raw_minus_1_5'])==F(889,181800)
reg=read('derivative_registration.json');s=read('derivative_summary.json');assert reg['code_sha256']==H((p/'derivative.py').read_bytes()) and s['registration_sha256']==H((p/'derivative_registration.json').read_bytes())
with gzip.open(p/'original_chunks.json.gz','rt') as f:data=json.load(f)
r={x['file']:x for x in s['receipts']};assert len(r)==len(data)==128
count=0;direct=0;baseline=0
for x in data:
 row=r[x['file']];assert H(x['original_utf8'].encode())==x['sha256']==row['sha256'];c=json.loads(x['original_utf8']);assert c['status']==row['status']=='completed' and c['states']==row['states']==1800;count+=c['states'];direct+=sum('direct_scaled_error' in z for z in c['independent_checks']);baseline+=sum('rho_scaled_mean_delta' in z for z in c['independent_checks'])
assert (count,direct,baseline)==(230400,1024,256)
print('128 raw chunks hashed, 230400 stated points, 1024 independent checks and 256 baseline checks; corrected candidate and initial sign transcription error retained')
print('Recorded values and hashes checked only; no new integration or global theorem established.')
