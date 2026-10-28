from pathlib import Path
import json,hashlib,gzip
p=Path(__file__).resolve().parent;H=lambda b:hashlib.sha256(b).hexdigest();read=lambda n:json.loads((p/n).read_text());m=read('archive_manifest.json')
for n,h in m['sha256'].items():assert H((p/n).read_bytes())==h
reg=read('frontier4_registration.json');s=read('frontier4_summary.json');assert reg['source_sha256']==H((p/'frontier4.py').read_bytes()) and s['registration_sha256']==H((p/'frontier4_registration.json').read_bytes());receipts={x['file']:x for x in s['receipts']}
with gzip.open(p/'original_chunks.json.gz','rt') as f:rows=json.load(f)
assert len(receipts)==len(rows)==1280;counts={'domain':0,'payment':0,'eligible':0};tops=[]
for x in rows:
 r=receipts[x['file']];assert H(x['original_utf8'].encode())==x['sha256']==r['sha256'];j=json.loads(x['original_utf8']);assert j['source_sha256']==reg['source_sha256'] and j['status']==r['status']=='completed' and j['counts']==r['counts'];tops.extend(j['top32'])
 for k in counts:counts[k]+=j['counts'][k]
assert sum(counts.values())==10485760;tops.sort(key=lambda z:z['cost_float'],reverse=True);assert tops[:64]==s['top64_float_proposals'];x=read('frontier4_recheck.json');assert x['source_sha256']==H((p/'frontier4_recheck.py').read_bytes()) and x['screen_summary_sha256']==H((p/'frontier4_summary.json').read_bytes());print('1280 raw chunks, 10485760 reported combinations, float top64 and post-run recheck identities verified:',counts);print('No scalar searches, high-precision evaluations, or operator theorems rerun.')
