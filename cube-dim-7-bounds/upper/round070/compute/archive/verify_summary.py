from pathlib import Path
import json,hashlib,gzip
p=Path(__file__).resolve().parent;H=lambda b:hashlib.sha256(b).hexdigest();read=lambda n:json.loads((p/n).read_text());m=read('archive_manifest.json')
for n,h in m['sha256'].items():assert H((p/n).read_bytes())==h
with gzip.open(p/'original_chunks.json.gz','rt') as f:allrows=json.load(f)
for stage,N,states,regname,src in [('pilot',14,157500,'registration.json','h2.py'),('production',242,2722500,'registration.json','h2.py'),('width',256,4096,'width_registration.json','finite_width.py'),('width_e',256,4096,'width_e_registration.json','finite_width_e.py')]:
 reg=read(regname);s=read(stage+'_summary.json');assert reg['source_sha256']==H((p/src).read_bytes()) and s['registration_sha256']==H((p/regname).read_bytes());rs={x['file']:x for x in s['receipts']};assert len(rs)==len(allrows[stage])==N;count=0
 for row in allrows[stage]:
  r=rs[row['file']];assert H(row['original_utf8'].encode())==row['sha256']==r['sha256'];j=json.loads(row['original_utf8']);assert j['source_sha256']==reg['source_sha256'] and j['status']==r['status']=='completed';count+=r['states']
 assert count==states;print(stage,N,'chunks',count,'states verified')
assert not ({x['file'] for x in allrows['pilot']}&{x['file'] for x in allrows['production']})
print('Original bytes, identities and reported counts checked; no numerical integrations rerun or operator theorem inferred.')
