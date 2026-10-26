#!/usr/bin/env python3
"""Independent integer-only cross-products for the same four saved endpoints."""
from pathlib import Path
import hashlib
import json
import os
import time

HERE = Path(__file__).resolve().parent


def pair(s):
    n, sep, d = s.partition('/')
    return int(n), int(d) if sep else 1


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


start = time.monotonic()
reg = json.loads((HERE/'registration.json').read_text())
result = json.loads((HERE/'results.json').read_text())
source = next(p for p in reg['sources'] if p.endswith('fixed_H64_obstruction_batch0009.json'))
assert sha(source) == reg['sources'][source]
raw = json.loads(Path(source).read_text())
gn, gd = pair(raw['target'])
for original, saved in zip(raw['counterexamples'], result['saved_counterexamples']):
    ln, ld = pair(original['H_interval'][0])
    un, ud = pair(original['H_interval'][1])
    rn, rd = pair(saved['strict_target_minus_H_upper'])
    assert 0 < ln*ud <= un*ld
    assert gn*ud - un*gd > 0
    assert (gn*ud - un*gd)*rd == rn*gd*ud
    assert saved['point'] == original['point']
    for coord in original['point'][1:]:
        n,d=pair(coord);assert 0 <= n <= d
assert len(raw['counterexamples']) == len(result['saved_counterexamples']) == 4
for den,row,abstract in zip((64,128),result['epsilon_constants'],result['abstract_two_leaf_cases']):
    assert row['epsilon'] == '1/'+str(den)
    assert int(row['hessian_magnitude'])*8 == den*den
    assert int(row['Taylor_penalty'])*16 == den*den
    assert int(row['log_score_bound']) == den
    assert int(abstract['minimum_log_margin'])*256 == 1024*256 - 13*den*den
pending=json.loads((HERE/'pending194.json').read_text())['nodes']
assert len(pending) == len({tuple(x['node']) for x in pending}) == 194
assert sum(x['origin']=='queue' for x in pending)==190
assert sum(x['origin']=='inflight_counter_H64' for x in pending)==4
assert all(x['new_epsilon_intended']=='1/128' and x['certified_in_099'] is False for x in pending)
receipt={'status':'passed','pid':os.getpid(),'wall_seconds':time.monotonic()-start,
         'method':'integer cross-products; no Fraction/Arb/kernel imports',
         'checked_saved_endpoints':4,'checked_epsilon_constants':2,'checked_abstract_leaves':2,
         'new_H_calls':0,'new_gradient_calls':0,'code_sha256':sha(__file__),
         'results_sha256':sha(HERE/'results.json'),'pending194_sha256':sha(HERE/'pending194.json')}
(HERE/'verification.json').write_text(json.dumps(receipt,sort_keys=True,indent=2)+'\n')
print(json.dumps(receipt,sort_keys=True))
