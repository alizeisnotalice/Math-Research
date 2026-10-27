#!/usr/bin/env python3
"""Saved arithmetic, inputs, alias roots and child receipts; union is not rerun."""
from pathlib import Path
from fractions import Fraction as F
import hashlib,json
H=Path(__file__).resolve().parent;sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest();reg=json.loads((H/'registration.json').read_text());lock=json.loads((H/'input_lock.json').read_text())
for n,h in lock['files'].items():assert sha(H/n)==h,n
assert sha(H.parent/'round003/theory_event_certificate.json')==reg['theory_certificate_sha256']
assert sha(H.parent/'round003/theory_manifest.json')==reg['theory_manifest_sha256']
basepath=H.parent/'round003/candidate_left.json';assert sha(basepath)==reg['baseline_r003_left_sha256'];base=json.loads(basepath.read_text())['result'];ref=F(base['union_interval'][1]);inherit=F(base['finite_L1']['lower']);D=10**35; rows=[]
assert len(reg['selected'])==14 and len(reg['skipped'])==2
for e in reg['selected']:
 n=e['rank'];file=H/f'candidate_{n:02}.json';d=json.loads(file.read_text());c=d['result'];rc=json.loads((H/f'receipt_{n:02}.json').read_text());assert d['event']==e and d['status']=='complete' and rc['returncode']==0 and not rc['stderr']
 assert rc['peak_rss_bytes']==d['peak_rss_bytes'] and d['input_lock_sha256']==sha(H/'input_lock.json')
 assert F(d['comparison_upper'])==ref and F(d['inherited_finite_lower'])==inherit
 a=e['selected_alias'];terms=[(a['p'],F(1,2)),(a['q'],F(1,2))] if a['kind']=='sum' else [(a['p'],F(1,2) if a['kind']=='r' else F(1))]
 assert len(d['event_roots'])==len(terms)
 lo=hi=F(0)
 for x,(m,k) in zip(d['event_roots'],terms):
  q=int(x['lower_numerator']);upper=int(x['upper_numerator']);exact=q**8==m*D**8
  assert x['radicand']==m and F(x['coefficient'])==k and q**8<=m*D**8<(q+1)**8 and upper==q+(not exact) and x['perfect_eighth_power']==exact
  lo+=4*k*q;hi+=4*k*upper
 assert list(map(int,c['spacing_endpoints']))==[lo,hi]
 assert c['states']==c['inner_states']+c['outer_states']<=reg['memo_limit_per_candidate']
 inner,outer=map(F,c['union_interval']);assert inner<=outer<=F(c['cell_upper'])
 fl=c['finite_L1'];eps=F(fl['epsilon']);theta=F(fl['theta']);K=fl['K'];loss=256*c['positive_inward_raw_boxes']*16*eps*2**7
 assert F(c['minimum_width'])>2*eps and loss==F(fl['loss_upper'])
 finite=theta*(inner-loss)*F(2*K+1,2*K+5)**8
 assert finite==F(fl['lower']) and F(fl['margin'])==finite-ref and c['winner']==(finite>ref)
 rows.append({'rank':n,'finite_lower':str(finite),'finite_minus_baseline':str(finite-ref),'inner_vs_baseline':(inner>ref)-(inner<ref),'outer_vs_baseline':(outer>ref)-(outer<ref),'states':c['states'],'returncode':rc['returncode'],'peak_rss_bytes':rc['peak_rss_bytes'],'candidate_sha256':sha(file)})
rows.sort(key=lambda x:F(x['finite_lower']),reverse=True);assert rows[0]['rank']==16 and F(rows[0]['finite_lower'])>F(rows[1]['finite_lower']) and F(rows[0]['finite_minus_baseline'])>0
print(json.dumps({'status':'PASS saved arithmetic and receipts; union not rerun','winner_rank':rows[0]['rank'],'ordered':rows},indent=2))
