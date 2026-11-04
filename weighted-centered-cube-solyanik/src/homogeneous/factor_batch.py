"""Parallel finite quadrature of mismatched lp-power/cube coefficients.

No quadrature difference here is a rigorous error enclosure. Checkpoints are
task-specific and resumable, and rejected/nonfinite outputs are never clamped.
"""
import os
for name in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS',
             'VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS'):
    os.environ[name] = '1'
import json
import math
import time
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed
from multiprocessing import get_context
from functools import lru_cache
import numpy as np
from scipy.special import roots_legendre

ROOT = Path(__file__).resolve().parent
ORDERS = (32, 64, 128, 256, 512, 1024)
PS = (1, 1.25, 1.5, 2, 3, 4, 8, 16, 'inf')
VERSION = 'cube-lp-v1-N1024'


def atomic_json(path, value):
    temp = path.with_suffix(f'.{os.getpid()}.tmp')
    temp.write_text(json.dumps(value, allow_nan=False, indent=2))
    temp.replace(path)


@lru_cache(None)
def nodes(order):
    x, w = roots_legendre(order)
    return (x+1)/2, w/2


@lru_cache(maxsize=96)
def face_log(n, p, order):
    x, _ = nodes(order)
    xp = x**p
    return np.log1p(xp if n == 2 else xp[:, None] + xp[None, :])


def solve(job):
    n, p, j = job
    beta = j/128
    tau = n*beta
    folder = ROOT/'checkpoints'/f'n{n}'/f'p{p}'
    folder.mkdir(parents=True, exist_ok=True)
    path = folder/f'j{j:03}.json'
    if path.exists():
        try:
            cached = json.loads(path.read_text())
            valid = (cached.get('version') == VERSION and
                     cached.get('status') == 'completed' and
                     (cached.get('n'),cached.get('p'),cached.get('j')) == job and
                     cached.get('orders') == list(ORDERS) and
                     len(cached.get('coefficients',[])) == len(ORDERS) and
                     all(math.isfinite(v) and v>0 for v in cached['coefficients']))
        except (OSError,ValueError,TypeError,KeyError):
            valid = False
        if valid:
            return dict(cached, run_pid=os.getpid(), reused_this_run=True)
    start = time.monotonic()
    result = dict(version=VERSION, status='running', n=n, p=p, j=j,
                  beta=beta, tau=tau, pid=os.getpid(),
                  run_pid=os.getpid(), reused_this_run=False)
    values = []
    if n == 1 or p == 'inf':
        values = [1.0]*len(ORDERS)
        samples = 0
    else:
        log_volume_ratio = n*math.lgamma(1+1/p)-math.lgamma(1+n/p)
        numerator = math.exp((1-beta)*log_volume_ratio)
        samples = 0
        for order in ORDERS:
            _, w = nodes(order)
            integrand = np.exp(((tau-n)/p)*face_log(n, p, order))
            integral = float(w @ integrand if n == 2 else w @ (integrand @ w))
            values.append(numerator/integral)
            samples += order**(n-1)
    if not all(math.isfinite(v) and v > 0 for v in values):
        result.update(status='failed', reason='nonfinite or nonpositive coefficient')
    else:
        result.update(status='completed', orders=list(ORDERS), coefficients=values,
                      last_difference=abs(values[-1]-values[-2]),
                      integrand_samples=samples, seconds=time.monotonic()-start)
    atomic_json(path, result)
    return result


def main():
    config = dict(version=VERSION, orders=list(ORDERS), ps=list(PS),
                  dimensions=[1,2,3], beta_numerators=[1,127], beta_denominator=128,
                  workers=15, pid=os.getpid(), started=time.time(),
                  numpy=np.__version__)
    atomic_json(ROOT/'config.json', config)
    jobs = [(n,p,j) for n in (1,2,3) for p in PS for j in range(1,128)]
    rows = []
    start = time.monotonic()
    with ProcessPoolExecutor(max_workers=15, mp_context=get_context('spawn')) as pool:
        futures = {pool.submit(solve,job):job for job in jobs}
        for future in as_completed(futures):
            job = futures[future]
            try:
                rows.append(future.result())
            except Exception as exc:
                rows.append(dict(status='failed', job=list(job), reason=repr(exc)))
            if len(rows)%250 == 0:
                state = dict(status='running', completed=len(rows), expected=len(jobs),
                             parent_pid=os.getpid(), elapsed=time.monotonic()-start)
                atomic_json(ROOT/'state.json', state)
                print(json.dumps(state), flush=True)
    failed = [r for r in rows if r['status'] != 'completed']
    completed = [r for r in rows if r['status'] == 'completed']
    summary = dict(status='completed' if not failed else 'failed', expected=len(jobs),
                   completed=len(completed), failed=failed,
                   workers_observed=sorted({r['run_pid'] for r in completed}),
                   computed_jobs=sum(not r['reused_this_run'] for r in completed),
                   reused_jobs=sum(r['reused_this_run'] for r in completed),
                   integrand_samples=sum(r['integrand_samples'] for r in completed),
                   new_integrand_samples=sum(r['integrand_samples'] for r in completed
                                             if not r['reused_this_run']),
                   max_last_difference=max(r['last_difference'] for r in completed),
                   below_one=[r for r in completed if r['coefficients'][-1] < 1-1e-10],
                   elapsed=time.monotonic()-start)
    atomic_json(ROOT/'summary.json', summary)
    atomic_json(ROOT/'state.json', summary)
    atomic_json(ROOT/'grid.json', sorted(completed,key=lambda r:(r['n'],str(r['p']),r['j'])))
    print(json.dumps(summary), flush=True)
    assert len(completed) == len(jobs) and not failed


if __name__ == '__main__':
    main()
