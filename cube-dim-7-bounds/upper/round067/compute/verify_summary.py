from pathlib import Path
import hashlib,json,gzip
p=Path(__file__).resolve().parent;H=lambda b:hashlib.sha256(b).hexdigest();read=lambda n:json.loads((p/n).read_text())
manifest=read('archive_manifest.json')
for n,h in manifest['sha256'].items():assert H((p/n).read_bytes())==h,(n,'hash mismatch')
c=read('candidate_certificate.json');assert c['code_sha256']==H((p/'candidate.py').read_bytes()) and len(c['verdict'])==9 and all(c['verdict'].values())
reg=read('registration.json');s=read('summary.json');assert reg['code_sha256']==H((p/'phi_diagnostic.py').read_bytes()) and s['registration_sha256']==H((p/'registration.json').read_bytes());receipts={x['file']:x for x in s['receipts']}
with gzip.open(p/'original_chunks.json.gz','rt') as f:rows=json.load(f)
assert len(receipts)==len(rows)==256
states=direct=baseline=0
for x in rows:
 r=receipts[x['file']];assert H(x['original_utf8'].encode())==x['sha256']==r['sha256'];j=json.loads(x['original_utf8']);assert j['source_sha256']==reg['code_sha256'] and j['status']==r['status']=='completed' and j['states']==r['states']==2250
 states+=j['states'];direct+=sum('scaled_value_error' in z for z in j['independent_checks']);b=[z for z in j['independent_checks'] if 'rho_scaled_mean_error' in z];assert all('rho2_scaled_second_moment' in z and 'rho2_scaled_variance' in z for z in b);baseline+=len(b)
assert (states,direct,baseline)==(576000,2048,512)
ratio=read('ratio_probe.json');assert ratio['source_sha256']==H((p/'ratio_probe.py').read_bytes()) and ratio['phi_source_sha256']==reg['code_sha256'] and ratio['states']==len(ratio['rows'])==1536
print('256 original raw chunks verified; 576000 states; 2048 direct checks; 512 baseline cases each include mean, second moment and variance')
print('9 scalar verdicts and 1536 post-run H/J records retained; no integrals or derivatives rerun, no global theorem inferred.')
