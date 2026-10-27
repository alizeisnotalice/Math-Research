#!/usr/bin/env python3
from pathlib import Path
from fractions import Fraction as F
import hashlib,json
H=Path(__file__).resolve().parent
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
def sig(terms):
 out={}
 for n,c in terms:
  a=b=0
  while n%2==0:n//=2;a+=1
  while n%3==0:n//=3;b+=1
  assert n==1
  k=(a%8,b%8);out[k]=out.get(k,F(0))+c*2**(a//8)*3**(b//8)
 return [[a,b,c.numerator,c.denominator] for (a,b),c in sorted(out.items())]
def terms(a):
 k,p,q=a['kind'],a['p'],a['q']
 return [(p,F(1,2)),(q,F(1,2))] if k=='sum' else [(p,F(1,2) if k=='r' else F(1))]
source=H.parent/'round003/theory_event_certificate.json';data=json.loads(source.read_text()); assert len(data['candidates'])==16
old={}
for m in [16,32,48,64,96,128,192,256]:old[tuple(map(tuple,sig([(m,F(1,2)),(2*m,F(1,2))])))]=f'old m={m}'
for n,p,q in [('left',64,256),('right',32,486)]:old[tuple(map(tuple,sig([(p,F(1,2)),(q,F(1,2))])))]=f'R003 {n}'
selected=[];skipped=[]
for c in data['candidates']:
 key=tuple(map(tuple,c['signature']))
 assert all(sig(terms(a))==c['signature'] for a in c['aliases'])
 item={'rank':c['rank'],'signature':c['signature'],'aliases':c['aliases'],'selected_alias':c['aliases'][0]}
 if key in old:item['reason']=old[key];skipped.append(item)
 else:selected.append(item)
assert len(selected)==14 and len(skipped)==2
reg={'scope':'R004 exactly new events from 16 certified recommendations; do not add candidates','theory_certificate_sha256':sha(source),'theory_manifest_sha256':sha(H.parent/'round003/theory_manifest.json'),'baseline_r003_left_sha256':sha(H.parent/'round003/candidate_left.json'),'reference':'R003 candidate_left result.union_interval upper, exact Fraction','inherited_lower':'R003 candidate_left result.finite_L1.lower','memo_limit_per_candidate':2**20,'parallel_limit':14,'execution':'collect all; no reruns; actual subprocess exit and elapsed; per-child RSS reported from getrusage','selected':selected,'skipped':skipped,'old_m16_prescreen':'excluded by R001 cell volume; no recomputation'}
with (H/'registration.json').open('x') as f:json.dump(reg,f,indent=2);f.write('\n')
print('selected',[(c['rank'],c['selected_alias']) for c in selected]);print('skipped',[(c['rank'],c['reason']) for c in skipped])
