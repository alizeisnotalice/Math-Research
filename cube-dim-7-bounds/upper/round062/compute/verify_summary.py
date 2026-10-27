import json,hashlib,gzip
from pathlib import Path
from decimal import Decimal
from fractions import Fraction as F
p=Path(__file__).resolve().parent
H=lambda data:hashlib.sha256(data).hexdigest()
read=lambda name:json.loads((p/name).read_text())
manifest=read('archive_manifest.json')
for name,sha in manifest['sha256'].items():assert H((p/name).read_bytes())==sha,(name,'hash mismatch')
with gzip.open(p/'original_chunks.json.gz','rt') as f:data=json.load(f)
for batch,ntasks,ntuple,nbands in (('frontier',128,4096,245760),('frontier_v2',256,131072,50331648)):
 reg=read(batch+'_registration.json');summary=read(batch+'_summary.json');assert reg['code_sha256']==H((p/(batch+'.py')).read_bytes());assert summary['registration_sha256']==H((p/(batch+'_registration.json')).read_bytes());rr={x['file']:x for x in summary['receipts']};assert len(rr)==len(data[batch])==ntasks
 tops=[];tuples=bands=eligible=0
 for z in data[batch]:
  x=rr[z['file']];assert H(z['original_utf8'].encode())==z['sha256']==x['sha256'];c=json.loads(z['original_utf8']);assert c['status']==x['status']=='completed' and c['tuples']==x['tuples'] and c['bands']==x['bands'] and c['eligible']==x['eligible'];tuples+=c['tuples'];bands+=c['bands'];eligible+=c['eligible'];tops.extend(c['top32'])
 assert tuples==ntuple and bands==nbands;tops.sort(key=lambda x:Decimal(x['c_model']),reverse=True);assert tops[:64]==summary['top64'];print(batch,'chunks',ntasks,'tuples',tuples,'band terms',bands,'eligible',eligible,'top64 checked')
reg=read('frontier_v2_registration.json');assert reg['old_pilot_fraction_overlap']==0
post=read('frontier_v2_postrun_manifest.json');assert post['source_sha256']==post['source_hash_actual']==H((p/'frontier_v2.py').read_bytes());assert post['registration_sha256']==H((p/'frontier_v2_registration.json').read_bytes()) and post['summary_sha256']==H((p/'frontier_v2_summary.json').read_bytes());assert post['chunk_sha256']=={x['file']:x['sha256'] for x in data['frontier_v2']}
c=read('candidate_certificate.json');assert c['code_sha256']==H((p/'candidate_cert.py').read_bytes());assert len(c['rows'])==128
r=c['rows'];assert all(F(x['eta'])>F(x['chi']) and F(x['sigma_margin'])>0 for x in r);g=sum((F(x['gamma_lower']) for x in r),F(0));assert g==F(c['gamma_total_lower']) and g-F(3,47)-F(3,20)==F(c['gamma_minus_3_47_minus_3_20_lower'])>0
gate=read('candidate_gate.json');assert gate['source_sha256']==H((p/'candidate_gate.py').read_bytes()) and gate['input_sha256']==H((p/'candidate_certificate.json').read_bytes());assert F(gate['gamma_lower'])==g and g>F(gate['gamma_threshold'])==F(429,2000);assert F(gate['raw_coefficient'])==F(14163,94000) and F(gate['raw_minus_3_20'])==F(63,94000)
t=read('frontier_v2_top_certificate.json');assert t['source_sha256']==H((p/'frontier_v2_top_cert.py').read_bytes()) and t['input_summary_sha256']==H((p/'frontier_v2_summary.json').read_bytes());assert t['top_parameters']=={k:read('frontier_v2_summary.json')['top64'][0][k] for k in ('d','qhi','J','i','k')};rows=t['rows'];assert len(rows)==t['top_parameters']['J'];assert all(F(x['eta'])-F(x['mu_upper'])==F(x['eta_margin'])>0 for x in rows);c0=sum((F(x['term_lower']) for x in rows),F(0))-2*max(F(x['mu_upper']) for x in rows);assert c0==F(t['c_model_lower'])>0 and min(F(x['eta_margin']) for x in rows)==F(t['min_eta_minus_mu_lower'])>0
print('candidate 128 layer margins, post-run clean gate and retained top certificate arithmetic verified')
print('Only retained bytes, hashes and recorded arithmetic are checked; omitted 50m term calculations are not rerun. No continuous/global optimum or operator theorem follows.')
