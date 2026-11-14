#!/usr/bin/env python3
"""Optimizer-free exact replay of posterior-hull capacity certificates.

Only Python standard library and adjacent r732/replay_exact.py are used.
Reads posterior_hull_capacity.json and ../r732/capacity.json; writes
replay_hull.json. A nonpositive gap does not certify equality.
"""
import hashlib
import importlib.util
import json
import sys
import time
from fractions import Fraction as Q
from pathlib import Path

HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location('r732_exact', HERE.parent/'r732'/'replay_exact.py')
exact = importlib.util.module_from_spec(spec)
spec.loader.exec_module(exact)


def decimal(x, upper=False, places=40):
    scale = 10**places
    value = x.numerator*scale // x.denominator
    if upper and Q(value,scale) < x:
        value += 1
    sign = '-' if value < 0 else ''
    value = abs(value)
    return sign+str(value//scale)+'.'+str(value % scale).zfill(places)


def check(base, r):
    exact.check(base)
    assert base['case'] == r['case'] and base['target'] == r['target']
    assert base['capacity_lower'] == r['full_capacity_lower']
    assert base['capacity_upper'] == r['full_capacity_upper']
    mass = list(map(Q,base['atom_masses']))
    vol = list(map(Q,base['row_volumes']))
    masks = base['atom_masks']
    pi = list(map(Q,base['pi']))
    active = base['active_rows']
    pir = [pi[i] for i in active]
    q = [v/sum(mass) for v in mass]
    kernel = [[m/vol[i] if s & (1 << i) else Q(0)
               for m,s in zip(mass,masks)] for i in active]
    pm = [sum(pir[i]*kernel[i][j] for i in range(len(active)))
          for j in range(len(mass))]
    assert min(pm)>0 and sum(pm)==1
    alpha = [[pir[i]*kernel[i][j]/pm[j] for i in range(len(active))]
             for j in range(len(mass))]
    assert all(min(row)>=0 and sum(row)==1 for row in alpha)
    assert [sum(pm[j]*alpha[j][i] for j in range(len(mass)))
            for i in range(len(active))] == pir
    z = list(map(Q,r['space_weights']))
    assert len(z)==len(mass) and min(z)>0 and sum(z)==1
    nu = [sum(z[j]*alpha[j][i] for j in range(len(mass)))
          for i in range(len(active))]
    assert nu == list(map(Q,r['induced_nu']))
    assert min(nu)>0 and sum(nu)==1
    pn = [sum(nu[i]*kernel[i][j] for i in range(len(active)))
          for j in range(len(mass))]
    assert min(pn)>0 and sum(pn)==1
    fl,fh = exact.lin([(x,x/y) for x,y in zip(pn,q)]
                      +[(-x,x/y) for x,y in zip(nu,pir)])
    g = [exact.lin([(kernel[i][j],pn[j]/q[j]) for j in range(len(mass))]
                   +[(Q(-1),nu[i]/pir[i])]) for i in range(len(active))]
    h = [(sum(alpha[j][i]*g[i][0] for i in range(len(active))),
          sum(alpha[j][i]*g[i][1] for i in range(len(active))))
         for j in range(len(mass))]
    upper = max(x[1] for x in h)
    # Euler identity checked as an interval enclosure, in addition to exact nu=z alpha.
    euler_low = sum(z[j]*h[j][0] for j in range(len(mass)))-fh
    euler_high = sum(z[j]*h[j][1] for j in range(len(mass)))-fl
    assert euler_low<=0<=euler_high
    assert Q(r['hull_capacity_lower'])<=fl
    assert Q(r['hull_capacity_upper'])>=upper
    assert Q(r['hull_interval_width_upper'])>=upper-fl
    gap = Q(base['capacity_lower'])-upper
    assert Q(r['full_minus_hull_lower'])<=gap
    if r['certified_strict_loss']:
        assert Q(r['full_minus_hull_lower'])>0 and gap>0
    return dict(case=r['case'],dimension=base['dimension'],target=r['target'],
                exact_base_LF_certificate=True,exact_atom_posteriors=True,
                exact_prior_in_hull=True,exact_positive_space_weights=True,
                exact_induced_nu=True,positive_induced_support=True,
                hull_lower_verified=True,hull_global_upper_verified=True,
                stated_hull_interval_width_verified=True,
                Euler_identity_interval_contains_zero=True,
                hull_lower_64=decimal(fl),hull_upper_64=decimal(upper,True),
                full_minus_hull_lower_64=decimal(gap),
                full_minus_hull_upper_64=decimal(Q(base['capacity_upper'])-fl,True),
                stated_strict_loss_verified=bool(r['certified_strict_loss']),
                strict_loss_certified_by_replay=gap>0)


def main():
    if not __debug__:
        raise RuntimeError('Run without -O or PYTHONOPTIMIZE; assertions are certificate checks.')
    started=time.monotonic()
    raw=(HERE.parent/'r732'/'capacity.json').read_bytes()
    bases=json.loads(raw)['results']
    candidate=json.loads((HERE/'posterior_hull_capacity.json').read_text())
    digest=hashlib.sha256(raw).hexdigest()
    assert digest==candidate['source_sha256']
    byid={r['case']:r for r in bases}
    assert len(byid)==len(bases)
    rows=candidate['results']
    assert len(rows)==len(bases) and {r['case'] for r in rows}==set(byid)
    results=[check(byid[r['case']],r) for r in rows]
    strict=[r['case'] for r in results if r['strict_loss_certified_by_replay']]
    out=dict(status='complete',all_pass=True,cases=len(results),
             strict_loss_cases=strict,log_series_terms=64,
             dependencies='Python standard library plus adjacent r732/replay_exact.py',
             source_sha256=digest,python=sys.version.split()[0],
             scope='Fixed LF pi; global maximum only over the convex hull of all atom posteriors. Does not optimize the full LF face. Nonpositive gaps do not prove equality; q=3/4 example is not a central asymptotic conclusion.',
             elapsed_seconds=time.monotonic()-started,results=results)
    (HERE/'replay_hull.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({k:out[k] for k in ['status','all_pass','cases','strict_loss_cases','log_series_terms','elapsed_seconds']}))


if __name__=='__main__':
    main()
