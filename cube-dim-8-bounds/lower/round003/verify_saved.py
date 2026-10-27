#!/usr/bin/env python3
"""Verify frozen arithmetic, hashes and child receipts; no union rerun."""
from pathlib import Path
from fractions import Fraction as F
import hashlib,json
H=Path(__file__).resolve().parent
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
lock=json.loads((H/'input_lock.json').read_text());reg=json.loads((H/'registration.json').read_text())
for k,v in lock['local_sha256'].items():assert sha(H/k)==v,k
old=H.parent/'round002/candidate_96.json';assert sha(old)==lock['r002_candidate_96_sha256'];base=json.loads(old.read_text())['candidate']; ref=F(base['union_interval'][1]); inherit=F(base['finite_L1']['lower']);D=10**35; out=[]
for e in reg['events']:
 name=e['name'];d=json.loads((H/f'candidate_{name}.json').read_text());c=d['result'];receipt=json.loads((H/f'receipt_{name}.json').read_text())
 assert d['status']=='complete' and d['event']==e and receipt['returncode']==0 and not receipt['stderr']
 assert d['input_lock_sha256']==sha(H/'input_lock.json') and F(d['comparison_upper'])==ref and F(d['inherited_finite_lower'])==inherit
 roots=list(map(int,d['event_root_lower_numerators']))
 for m,q in zip(e['roots'],roots):assert q**8<=m*D**8<(q+1)**8
 assert list(map(int,c['spacing_endpoints']))==[2*sum(roots),2*(sum(roots)+2)]
 assert not c['cell_prefilter_excluded'] and c['states']==c['inner_states']+c['outer_states']<=2**20
 lo,hi=map(F,c['union_interval']);assert lo<=hi<=F(c['cell_upper'])
 if c.get('winner'):
  x=c['finite_L1'];eps=F(x['epsilon']);theta=F(x['theta']);K=x['K'];loss=256*c['positive_inward_raw_boxes']*16*eps*2**7
  assert loss==F(x['loss_upper']) and F(c['minimum_width'])>2*eps
  finite=theta*(lo-loss)*F(2*K+1,2*K+5)**8
  assert finite==F(x['lower']) and F(x['margin'])==finite-ref and finite>ref and lo>ref
 else:assert hi<ref
 out.append({'event':name,'inner_vs_reference':(lo>ref)-(lo<ref),'outer_vs_reference':(hi>ref)-(hi<ref),'winner':c.get('winner',False),'states':c['states'],'returncode':receipt['returncode']})
assert [r['event'] for r in out if r['winner']]==['left']
print(json.dumps({'status':'PASS saved arithmetic, event roots, input hashes and child receipts; union not rerun','events':out},indent=2))
