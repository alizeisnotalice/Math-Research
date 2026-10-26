"""Hash-only reinspection of frozen B062 and author sources, no mathematics rerun."""
import hashlib,json
from pathlib import Path
P=Path(__file__).resolve().parent
rd=lambda p:json.loads(Path(p).read_text())
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
r=rd(P/'registration.json');m=rd(P/'manifest.json');C=Path(r['producer_directory']);count=0
def pin(p,h):
 global count
 assert sha(p)==h,('hash mismatch',str(p));count+=1
for f,h in m['whitelist_sha256'].items():pin(P/f,h)
pin(C/'manifest.json',r['producer_manifest_sha256']);pin(C/'report.md',r['producer_report_sha256'])
a=rd(C/'manifest.json')
for f,h in a['whitelist_sha256'].items():pin(C/f,h)
for f,h in a['source_sha256_unchanged'].items():pin(f,h)
e=rd(P/'execution_receipt.json');g=rd(P/'process_receipt.json')
pin(P/'results.json',e['results_sha256']);pin(P/'input_manifest.json',e['input_manifest_sha256'])
assert e['exit_code']==g['exit_code']==0 and e['completed_configurations']==8
assert g['status']=='COMPLETED' and g['child_reaped'] and not g['active_handles']
print(json.dumps({'status':'PASS_HASH_ONLY_SEAL','hash_checks':count,'new_scientific_configurations':0,'new_kernel_calls':0}))
