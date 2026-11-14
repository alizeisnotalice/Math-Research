#!/usr/bin/env python3
"""Replay stored LF and capacity certificates using only the Python standard library.

Run from any working directory: python3 /path/to/replay_exact.py
Reads capacity.json beside this file and writes replay_exact.json beside it.
No optimization is performed. Checks fixed-pi capacity certificates, not a
maximum over the entire LF face. Log enclosures use 64 exact rational terms.
Geometry reconstruction remains the separate replay.wl check; this script
uses the stored finite member-atom geometry.
"""
import json
import sys
import time
from fractions import Fraction as Q
from functools import lru_cache
from pathlib import Path

TERMS = 64
HERE = Path(__file__).resolve().parent


@lru_cache(None)
def unit(x):
    assert 1 <= x <= 2
    z = (x - 1) / (x + 1)
    lower = sum((2 * z ** (2*j + 1) / Q(2*j + 1)
                 for j in range(TERMS)), Q(0))
    upper = lower + 2 * z ** (2*TERMS + 1) / ((2*TERMS + 1) * (1-z*z))
    return lower, upper


@lru_cache(None)
def logi(x):
    assert x > 0
    y, e = x, 0
    while y >= 2:
        y /= 2
        e += 1
    while y < 1:
        y *= 2
        e -= 1
    a, b = unit(y)
    c, d = unit(Q(2))
    return (a + e*c, b + e*d) if e >= 0 else (a + e*d, b + e*c)


def lin(terms):
    low = high = Q(0)
    for coeff, x in terms:
        a, b = logi(x)
        low += coeff * (a if coeff >= 0 else b)
        high += coeff * (b if coeff >= 0 else a)
    return low, high


def check(r):
    mass = list(map(Q, r['atom_masses']))
    vol = list(map(Q, r['row_volumes']))
    masks = r['atom_masks']
    pi = list(map(Q, r['pi']))
    nu = list(map(Q, r['nu']))
    f = list(map(Q, r['source_f']))
    ell = list(map(Q, r['source_ell']))
    om, price, p = sum(mass), Q(r['exact_NP_price']), Q(r['target'])
    k = len(pi)
    assert all(x > 0 for x in mass + vol)
    assert om == Q(r['union_volume'])
    assert len(mass) == len(masks) == len(f) == len(ell)
    assert len(vol) == k
    assert sum(pi) == sum(nu) == 1 and min(pi) >= 0 and min(nu) > 0
    assert r['active_rows'] == [i for i in range(k) if pi[i] > 0]
    assert len(nu) == len(r['active_rows'])
    for subset in range(1, 1 << k):
        lhs = sum(pi[i] for i in range(k) if subset & (1 << i))
        rhs = sum(m for m, s in zip(mass, masks) if s & subset) / om
        assert lhs <= rhs
    av = [sum(m*f[j] for j, (m, s) in enumerate(zip(mass, masks))
              if s & (1 << i)) / vol[i] for i in range(k)]
    assert all(0 <= x <= 1 for x in f + ell)
    assert all(ell[j] <= av[i] for j, s in enumerate(masks)
               for i in range(k) if s & (1 << i))
    assert sum(m*x for m, x in zip(mass, ell)) >= p*om
    assert sum(m*x for m, x in zip(mass, f)) == price
    dens = [sum(pi[i]/vol[i] for i in range(k) if s & (1 << i))
            for s in masks]
    assert min(dens) > 0
    need, cost = p, Q(0)
    for j in sorted(range(len(mass)), key=lambda j: dens[j], reverse=True):
        take = min(need, mass[j]*dens[j])
        cost += take/dens[j]
        need -= take
    assert need == 0 and cost == price
    active = r['active_rows']
    pir = [pi[i] for i in active]
    kernel = [[m/vol[i] if s & (1 << i) else Q(0)
               for m, s in zip(mass, masks)] for i in active]
    assert all(sum(row) == 1 for row in kernel)
    pn = [sum(nu[i]*kernel[i][j] for i in range(len(nu)))
          for j in range(len(mass))]
    q = [m/om for m in mass]
    assert min(pn) > 0
    fl, fh = lin([(v, v/qj) for v, qj in zip(pn, q)]
                 + [(-v, v/pj) for v, pj in zip(nu, pir)])
    gi = [lin([(kernel[i][j], pn[j]/q[j]) for j in range(len(q))]
              + [(Q(-1), nu[i]/pir[i])]) for i in range(len(nu))]
    gh = max(v[1] for v in gi)
    assert Q(r['capacity_lower']) <= fl
    assert Q(r['capacity_upper']) >= gh
    assert Q(r['certified_capacity_gap_upper']) >= gh-fl
    vl, vh = logi(om/price)
    assert Q(r['test_minus_capacity_lower']) <= vl-gh
    assert Q(r['test_minus_capacity_upper']) >= vh-fl
    pl = [sum(pir[i]*kernel[i][j] for i in range(len(nu)))
          for j in range(len(q))]
    kl, kh = lin([(v, v/qj) for v, qj in zip(pl, q)])
    assert Q(r['KL_pi_lower']) <= kl and Q(r['KL_pi_upper']) >= kh
    return dict(case=r['case'], exact_Hall=True, positive_support=True,
                exact_primal_feasible=True, exact_NP_price=str(cost),
                exact_LF_primal_NP_equality=True,
                capacity_interval_verified=True, KL_interval_verified=True,
                test_minus_capacity_interval_verified=True,
                stated_gap_verified=True)


def main():
    # Assertions are the certificate checks; never silently disable them.
    if not __debug__:
        raise RuntimeError('Run without -O or PYTHONOPTIMIZE; assertions are required.')
    start = time.monotonic()
    source = HERE / 'capacity.json'
    data = json.loads(source.read_text())
    results = [check(r) for r in data['results']]
    out = dict(status='complete', all_pass=True, cases=len(results),
               source=source.name, python=sys.version.split()[0],
               dependencies='Python standard library only',
               log_series_terms=TERMS,
               distinct_log_arguments=logi.cache_info().currsize,
               scope='Exact replay of stored finite member-atom LF and fixed-pi capacity certificates; no optimizer, no geometry reconstruction, no optimization over the LF face.',
               method='Exact Fraction arithmetic; independent 64-term atanh-log implementation with positive geometric remainder and signed interval propagation.',
               elapsed_seconds=time.monotonic()-start, results=results)
    target = HERE / 'replay_exact.json'
    target.write_text(json.dumps(out, indent=2) + '\n')
    print(json.dumps({key: out[key] for key in
                      ['status', 'all_pass', 'cases', 'log_series_terms',
                       'distinct_log_arguments', 'elapsed_seconds']}))


if __name__ == '__main__':
    main()
