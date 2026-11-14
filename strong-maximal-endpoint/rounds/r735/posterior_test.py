#!/usr/bin/env python3
"""Deterministic exact posterior-candidate test on r732 real catalogues.

Uses the r732 independent 64-term Fraction log interval implementation.
No optimizer, random sampling, or claim that posterior maxima equal capacity.
"""
import importlib.util
import json
import multiprocessing as mp
import os
import sys
import time
from fractions import Fraction as F
from pathlib import Path

HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location('r732_exact', HERE.parent / 'r732' / 'replay_exact.py')
exact = importlib.util.module_from_spec(spec)
spec.loader.exec_module(exact)


def decimal(x, upper=False, places=40):
    scale = 10**places
    integer = x.numerator*scale // x.denominator
    if upper and F(integer, scale) < x:
        integer += 1
    sign = '-' if integer < 0 else ''
    integer = abs(integer)
    return sign + str(integer//scale) + '.' + str(integer % scale).zfill(places)


def bounds(pair):
    return [decimal(pair[0]), decimal(pair[1], True)]


def subtract(a, b):
    return a[0]-b[1], a[1]-b[0]


def maximum(pairs):
    return max(x[0] for x in pairs), max(x[1] for x in pairs)


def run(r):
    started = time.monotonic()
    exact.check(r)
    mass = list(map(F, r['atom_masses']))
    masks = r['atom_masks']
    vol = list(map(F, r['row_volumes']))
    pi = list(map(F, r['pi']))
    source = list(map(F, r['source_f']))
    target = F(r['target'])
    omega = sum(mass)
    kernel = [[mass[j]/vol[i] if masks[j] & (1 << i) else F(0)
               for j in range(len(mass))] for i in range(len(pi))]
    q = [x/omega for x in mass]
    pp = [sum(pi[i]*kernel[i][j] for i in range(len(pi)))
          for j in range(len(mass))]
    assert min(pp) > 0
    assert sum(source[j]*pp[j] for j in range(len(mass))) == target
    canonical = [F(0)]*len(mass)
    remaining = target
    for j in sorted(range(len(mass)), key=lambda j: pp[j]/mass[j], reverse=True):
        take = min(remaining, pp[j])
        canonical[j] = take/pp[j]
        remaining -= take
    assert remaining == 0
    assert sum(mass[j]*canonical[j] for j in range(len(mass))) == F(r['exact_NP_price'])
    kl = exact.lin([(x, x/y) for x,y in zip(pp,q)])
    atoms = []
    f_intervals = []
    info_intervals = []
    for a in range(len(mass)):
        alpha = [pi[i]*kernel[i][a]/pp[a] for i in range(len(pi))]
        assert sum(alpha) == 1 and min(alpha) >= 0
        mixture = [sum(alpha[i]*kernel[i][b] for i in range(len(pi)))
                   for b in range(len(mass))]
        ft = [(x,x/y) for x,y in zip(mixture,q) if x]
        ft += [(-x,x/pi[i]) for i,x in enumerate(alpha) if x]
        fi = exact.lin(ft)
        it = [(alpha[i]*kernel[i][b],kernel[i][b]/mixture[b])
              for i in range(len(pi)) for b in range(len(mass))
              if alpha[i] and kernel[i][b]]
        information = exact.lin(it)
        logw = exact.logi(pp[a]/q[a])
        residual = (fi[0]-logw[1]+information[0],
                    fi[1]-logw[0]+information[1])
        assert residual[0] <= 0 <= residual[1]
        f_intervals.append(fi)
        info_intervals.append(information)
        atoms.append(dict(atom_index=a,mask=masks[a],mass=str(mass[a]),
                          w=str(pp[a]/q[a]),posterior=[str(x) for x in alpha],
                          in_certified_NP_source=bool(source[a]),
                          canonical_NP_fraction=str(canonical[a]),
                          F_interval=bounds(fi),information_interval=bounds(information),
                          identity_residual_interval=bounds(residual)))
    all_max = maximum(f_intervals)
    source_max = maximum([x for j,x in enumerate(f_intervals) if source[j]])
    canonical_max = maximum([x for j,x in enumerate(f_intervals) if canonical[j]])
    capacity = (F(r['capacity_lower']),F(r['capacity_upper']))
    assert all_max[1] <= capacity[1]
    diff = subtract(source_max,kl)
    avg_info = (sum(source[j]*pp[j]/target*info_intervals[j][0] for j in range(len(mass))),
                sum(source[j]*pp[j]/target*info_intervals[j][1] for j in range(len(mass))))
    return dict(case=r['case'],dimension=r['dimension'],target=str(target),
                certified_NP_source=[str(x) for x in source],
                canonical_NP_source=[str(x) for x in canonical],
                KL_pi_interval=bounds(kl),capacity_interval=bounds(capacity),
                max_certified_NP_posterior_interval=bounds(source_max),
                max_canonical_NP_posterior_interval=bounds(canonical_max),
                max_all_atom_posterior_interval=bounds(all_max),
                certified_NP_posterior_minus_KL_interval=bounds(diff),
                canonical_NP_posterior_minus_KL_interval=bounds(subtract(canonical_max,kl)),
                all_atom_posterior_minus_KL_interval=bounds(subtract(all_max,kl)),
                capacity_minus_KL_interval=bounds(subtract(capacity,kl)),
                capacity_minus_certified_NP_posterior_interval=bounds(subtract(capacity,source_max)),
                capacity_minus_all_atom_posterior_interval=bounds(subtract(capacity,all_max)),
                NP_mean_information_interval=bounds(avg_info),
                NP_posterior_certifiably_improves_KL=(diff[0]>0),
                NP_posterior_certifiably_below_KL=(diff[1]<0),
                full_capacity_certifiably_above_KL=(capacity[0]>kl[1]),
                atom_identity_checks_pass=True,
                atoms=atoms,worker_pid=os.getpid(),elapsed_seconds=time.monotonic()-started)


def main():
    if not __debug__:
        raise RuntimeError('Run without -O/PYTHONOPTIMIZE: certificate assertions are required.')
    started = time.monotonic()
    source = HERE.parent / 'r732' / 'capacity.json'
    records = json.loads(source.read_text())['results']
    with mp.get_context('spawn').Pool(processes=12) as pool:
        spawned_pids = [p.pid for p in pool._pool]
        rows = pool.map(run, records, chunksize=1)
    summary = dict(cases=len(rows),
                   certified_NP_posterior_improves_KL=sum(r['NP_posterior_certifiably_improves_KL'] for r in rows),
                   certified_NP_posterior_below_KL=sum(r['NP_posterior_certifiably_below_KL'] for r in rows),
                   full_capacity_above_KL=sum(r['full_capacity_certifiably_above_KL'] for r in rows),
                   all_atom_identity_checks=all(r['atom_identity_checks_pass'] for r in rows))
    out = dict(status='complete',scope='Finite real r732 catalogues; posterior candidates only, not a maximization over all nu or the LF face; no near-central evidence.',
               workers=12,spawned_worker_pids=spawned_pids,
               used_worker_pids=sorted({r['worker_pid'] for r in rows}),
               deterministic=True,seeds=None,python=sys.version.split()[0],
               dependencies='Python standard library plus adjacent r732/replay_exact.py',
               log_series_terms=64,source='../r732/capacity.json',
               elapsed_seconds=time.monotonic()-started,summary=summary,results=rows)
    (HERE/'posterior_test.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(dict(status=out['status'],workers=12,used_workers=len(out['used_worker_pids']),
                         elapsed_seconds=out['elapsed_seconds'],**summary)))


if __name__ == '__main__':
    main()
