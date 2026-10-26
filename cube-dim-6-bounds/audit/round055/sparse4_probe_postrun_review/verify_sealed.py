"""Hash-only closure; does not invoke any scientific checker or native process."""
import hashlib,json
from pathlib import Path
P=Path(__file__).resolve().parent
read=lambda p:json.loads(Path(p).read_text())
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
m=read(P/'manifest.json');r=read(P/'registration.json');cfg=read(P/'input_manifest.json')
C=Path(cfg['producer_directory']);B=P.parent/'sparse4_probe_precheck_review'
n=0
def pin(p,h):
 global n
 assert sha(p)==h,('changed sealed file',str(p))
 n+=1
for f,h in m['whitelist_sha256'].items():pin(P/f,h)
pin(C/'full_manifest.json',r['producer_full_manifest_sha256'])
f=read(C/'full_manifest.json')
for p,h in f['whitelist_sha256'].items():pin(C/p,h)
pin(C/'precheck_manifest.json',f['precheck_manifest_sha256'])
pre=read(C/'precheck_manifest.json')
for p,h in pre['whitelist_sha256'].items():pin(C/p,h)
for p,h in pre['external_readonly_sha256'].items():pin(p,h)
pin(B/'manifest.json',r['precheck_manifest_sha256'])
for p,h in read(B/'manifest.json')['whitelist_sha256'].items():pin(B/p,h)
for p,h in read(B/'input_manifest.json')['sha256'].items():pin(p,h)
for p,h in read(B/'header_binding.json')['headers_sha256'].items():pin(p,h)
for p,h in cfg['additional_readonly_sha256'].items():pin(p,h)
e=read(P/'execution_receipt.json')
pin(P/'input_manifest.json',e['input_manifest_sha256']);pin(P/'results.json',e['results_sha256'])
assert e['exit_code']==0 and e['status']=='PASS_SAVED_ONLY_FIXED_FRAME_TIME_OBSTRUCTION'
print(json.dumps({'status':'PASS_HASH_ONLY_POSTRUN_SEAL','hash_entries_verified':n,'new_scientific_configurations':0,'native_calls':0}))
