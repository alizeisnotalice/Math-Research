#!/usr/bin/env python3
"""Read-only saved-arithmetic and file-binding checker; stdout only.

This does not claim an independent union recomputation. For a fresh full union
run copy registration.json, frozen_internal_reference.json and certify.py to
an empty directory and run certify.py there. The independent audit is separate.
"""
from pathlib import Path
from fractions import Fraction as F
import argparse, hashlib, json

def h(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--root',type=Path,default=Path(__file__).resolve().parent)
    args=ap.parse_args(); root=args.root.resolve()
    manifest=json.loads((root/'manifest.json').read_text())
    for name,digest in manifest['files'].items(): assert h(root/name)==digest,name
    ref=json.loads((root/'frozen_internal_reference.json').read_text())
    refl,refu=map(F,ref['interval']); report=[]
    for seed in (4,5):
        j=json.loads((root/f'seed{seed}_results.json').read_text())
        assert j['self_sha256']==h(root/'certify.py')
        assert j['registration_sha256']==h(root/'registration.json')
        assert j['reference_sha256']==h(root/'frozen_internal_reference.json')
        assert j['raw_length_vectors']==8192 and j['deduplicated_states']<=2**20
        d=int(j['root_denominator'])
        for m,num in j['root_lower_numerators'].items():
            q=int(num); assert q**7<=int(m)*d**7<(q+1)**7
        q=int(j['first_root_lower_numerator']); assert q**7<=2**seed*d**7<(q+1)**7
        lo,hi=map(F,(j['union_lower'],j['union_upper']))
        assert refu==F(j['complete_R075_upper']) and lo>refu and lo<=hi
        x=j['finite_L1']; eps=F(x['epsilon']); eta=F(x['eta'])
        assert eta==eps/2 and F(j['min_inward_width'])>2*eps
        loss=128*j['positive_raw_boxes']*14*eps*2**6
        assert loss==F(x['union_shrink_loss_upper'])
        k=x['K']; finite=F(x['theta'])*(lo-loss)*F(2*k+1,2*k+5)**7
        assert finite==F(x['finite_ratio_lower']) and finite>refu
        if seed==5: assert finite>F(2135677,10**6)
        report.append({'seed':seed,'full_union_lower':str(lo),'full_union_upper':str(hi),'finite_L1_lower':str(finite),'strictly_above_R075_upper':True})
    print(json.dumps({'status':'PASS saved arithmetic and hash bindings; independent union audit separate','root':str(root),'simple_finite_L1_bound':'2135677/1000000','seeds':report},indent=2))

if __name__=='__main__': main()
