import json,hashlib,gzip
from pathlib import Path
from decimal import Decimal
from fractions import Fraction as F
p=Path(__file__).resolve().parent
H=lambda b:hashlib.sha256(b).hexdigest()
read=lambda n:json.loads((p/n).read_text())
manifest=read('archive_manifest.json')
for name,h in manifest['sha256'].items():assert H((p/name).read_bytes())==h,(name,'hash mismatch')
reg=read('registration.json');summary=read('summary.json');assert summary['registration_sha256']==H((p/'registration.json').read_bytes());assert reg['code_sha256']==H((p/manifest['source']).read_bytes())
with gzip.open(p/'compact_receipts.json.gz','rt') as f:data=json.load(f)
n=manifest['round']
if n==56:
 assert len(data)==len(summary['results'])==60
 for x,y in zip(data,summary['results']):assert x['file']==y['file'] and x['original_sha256']==y['sha256'] and x['source_states']==y['source_states'] and x['status']==y['status']
 assert all(x['status']=='completed' for x in data);assert sum(x['source_states'] for x in data)==15360
 print('R056: 60 compact receipts; 15360 source states; original per-state rows omitted')
elif n==57:
 old=data['initial_compact'];assert len(old)==len(summary['results'])==240
 for x,y in zip(old,summary['results']):assert x['file']==y['file'] and x['original_sha256']==y['sha256'] and x['status']==y['status']
 assert all(x['status']=='exception_or_candidate_counterexample' for x in old)
 new=data['corrected_raw'];assert len(new)==240;vreg=read('corrected_v2/registration.json');assert vreg['code_sha256']==H((p/'corrected_v2/r057_kinks_v2.py').read_bytes());byphase={}
 for phase in ('pilot','run'):
  s=read('corrected_v2/summary_'+phase+'.json');assert s['registration_sha256']==H((p/'corrected_v2/registration.json').read_bytes());byphase[phase]={x['file']:x for x in s['results']}
 assert [len(byphase[x]) for x in ('pilot','run')]==[16,224];assert not set(byphase['pilot'])&set(byphase['run'])
 for x in new:
  r=byphase[x['phase']][x['file']];assert H(x['original_utf8'].encode())==r['sha256']==x['original_sha256'];j=json.loads(x['original_utf8']);assert j['status']=='completed' and j['states_completed']==1656 and len(j['direct_checks'])==8
 print('R057: 240 initial exception receipts; corrected 16+224 disjoint, 397440 states and 1920 direct checks; single-point precision diagnosis is author narrative')
elif n==58:
 receipts={x['file']:x for x in summary['receipts']};assert len(data)==len(receipts)==128;alltops=[]
 for x in data:
  r=receipts[x['file']];assert H(x['original_utf8'].encode())==x['original_sha256']==r['sha256'];j=json.loads(x['original_utf8']);assert j['status']=='completed' and j['checked']==16384 and j['checked']==r['checked'] and j['eligible']==r['eligible'] and j['positive']==r['positive'];alltops.extend(j['top32'])
 alltops.sort(key=lambda x:Decimal(x['c_model']),reverse=True);assert alltops[:64]==summary['top64'];assert sum(x['checked'] for x in receipts.values())==2097152
 cert=read('model_certificate.json');v={k:F(x) for k,x in cert['input'].items() if k in ('d','a','b')};d,a,b=v['d'],v['a'],v['b'];q=1-b*d;Q=1-a*d;R=q+d
 def lint(x):
  y=(x-1)/(x+1);s=2*sum((y**(2*k+1)/F(2*k+1) for k in range(48)),F(0));return s,s+2*y**97/(97*(1-y*y))
 l,u=lint(1/(1-d/(2*R)));ml=R/q*l*l;mu=R/q*u*u;lo,hi=lint((Q+d)/R);cl=lo*q/Q-2*mu;cu=hi*q/Q-2*ml;eta=a*d/(2*Q)
 assert (ml,mu,eta,cl,cu)==tuple(F(cert[k]) for k in ('mu_lower','mu_upper','eta','c_lower','c_upper'));assert mu<eta and cl>0
 print('R058: 128 original chunk bytes/hashes, 2097152 reported tuples, retained top32/top64 ranking and rational model certificate checked; omitted tuples not rerun')
print('retained archive file hashes verified; producer records do not replace independent audit')
