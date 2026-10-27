#!/usr/bin/env python3
"""Read-only, stdout-only saved arithmetic/hash checker; union audit separate."""
from pathlib import Path
from fractions import Fraction as F
import argparse,hashlib,json
def h(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--root',type=Path,default=Path(__file__).resolve().parent)
    root=ap.parse_args().root.resolve()
    manifest=json.loads((root/'manifest.json').read_text())
    for p,d in manifest['files'].items():assert h(root/p)==d,p
    j=json.loads((root/'results.json').read_text());f=json.loads((root/'frozen_round001.json').read_text())
    assert j['self_sha256']==h(root/'scan_certify.py')
    for p,d in j['input_hashes'].items():assert h(root/p)==d,p
    d=int(j['root_denominator'])
    for m,q in j['root_lower_numerators'].items():
        q=int(q);assert q**7<=int(m)*d**7<(q+1)**7
    event=j['event']; low=high=0
    for m,c in event['key']:
        q=int(j['event_radical_root_lower_numerators'][str(m)])
        assert q**7<=m*d**7<(q+1)**7
        low+=c*q;high+=c*(q+1)
    assert low==int(event['a1_lower_numerator']) and high==int(event['a1_upper_numerator'])
    assert 3*d<=low<=high<=4*d
    lo,hi=map(F,j['full_construction_comparison_interval']);base=F(f['full_interval'][1])
    assert base==F(j['complete_round001_upper']) and lo<=hi and lo>base
    assert lo-base==F(j['strict_margin_vs_complete_round001_upper'])
    x=j['finite_L1'];eps=F(x['epsilon']);theta=F(x['theta']);k=x['K']
    assert F(x['eta'])==eps/2 and F(j['minimum_inward_width'])>2*eps
    loss=128*j['positive_inward_raw_boxes']*14*eps*2**6
    assert loss==F(x['union_shrink_loss_upper'])
    finite=theta*(lo-loss)*F(2*k+1,2*k+5)**7
    assert finite==F(x['exact_ratio_lower']) and finite>base
    assert finite>F(2136208981884,10**12)
    assert finite-base==F(x['strict_margin_vs_complete_round001_upper'])
    t=json.loads((root/'trials.json').read_text());assert t['old_seed_evaluations']==0
    assert t['new_event_evaluations']==len(t['trials'])<=512
    for i,z in enumerate(t['trials']):
        assert z['event_key'] not in ([[1,2],[16,1]],[[1,2],[32,1]])
        above=F(z['full_union_interval'][0])>base
        assert above==z['strict_above_complete_round001_upper']
        assert above==(i==len(t['trials'])-1)
        assert z['states']<=2**20
    assert j['global_memo_states']==t['global_states']<=2**20
    print(json.dumps({'status':'PASS saved arithmetic/hash checks; independent union audit separate','root':str(root),'simple_proved_bound':'C7 > 2136208981884/10^12','first_spacing_event':event['key'],'finite_L1_exact_lower':str(finite),'full_construction_comparison_interval':[str(lo),str(hi)],'new_event_evaluations':len(t['trials']),'old_seed_evaluations':0,'global_memo_states':j['global_memo_states']},indent=2))
if __name__=='__main__':main()
