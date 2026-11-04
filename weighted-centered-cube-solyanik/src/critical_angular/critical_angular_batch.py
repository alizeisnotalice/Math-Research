#!/usr/bin/env python3
"""120-digit diagnostics for the critical angular logarithmic envelope.

This generates a full table locally (not committed).  It checks only the
stated scalar inverse/envelope and witness-mass formulae; it is not a proof
of the all-input maximal-operator theorem.
"""
import os
for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS'):
    os.environ[key] = '1'
import hashlib, json, platform, sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from decimal import Decimal as D, getcontext
from multiprocessing import Manager
from pathlib import Path

ROOT = Path(__file__).resolve().parent
CACHE = ROOT / 'runtime_cache'
SUMMARY = ROOT / 'round30_summary.json'
VERSION = 'critical-angular-round30-v2-geometry-20261105'
PREC = 120
NS, GS, KS = (2, 3), ('.25', '.5', '.75', '.9'), (8,16,32,64,128,256,512,1024,4096)

def source_hash():
    return hashlib.sha256(Path(__file__).read_bytes()).hexdigest()

def pw(x, exponent):
    return (exponent * x.ln()).exp()

def inverse_L(H, c):
    """Monotone bisection for L-log(1+cL)=H."""
    lo, hi = D(0), H + (1+c*(H+2)).ln() + 2
    def g(L): return L - (1+c*L).ln() - H
    assert g(lo) < 0 < g(hi)
    for _ in range(600):
        mid = (lo + hi) / 2
        if g(mid) < 0: lo = mid
        else: hi = mid
    return (lo + hi) / 2, g((lo + hi) / 2)

def envelope(x, c, beta, b):
    L, _ = inverse_L(-x.ln(), c)
    return pw((-L).exp(), beta) * (1+b*L)

def data(n, gs, k):
    getcontext().prec = PREC
    gamma = D(gs); beta = 1-gamma; tau = D(n)*beta
    c = D(n-1)/n; b = D(n-1)*beta/(D(n)-gamma)
    kappa = b*pw(c, -beta)
    H = D(k)*D(2).ln(); u = (-H).exp(); L, residual = inverse_L(H, c)
    F = pw((-L).exp(), beta)*(1+b*L)
    normalizer = pw(u,beta)*pw(H,gamma)
    C = pw(H,D(1)/(2*n)); delta = 1/H
    upper_u = u*pw(C,n)
    # The legacy line is just an asymptotic comparison.  The strict witness
    # in the current proof is the (1-2 delta)^n line below.
    legacy_u = u*pw(1-delta,n)/pw(delta,n)
    strict_u = u*pw(1-2*delta,n)/pw(delta,n)
    assert upper_u < 1 and legacy_u < 1 and strict_u < 1
    upper = pw(C,-tau)*envelope(upper_u,c,beta,b)
    legacy = pw(delta,tau)*envelope(legacy_u,c,beta,b)
    strict = pw(delta,tau)*envelope(strict_u,c,beta,b)
    # A future loss of bisection accuracy must not become an exit-zero row.
    assert abs(residual) < D('1e-100')
    return {'n':n, 'gamma':gs, 'k':k, 'inverse_residual':str(residual),
      'kappa':str(kappa), 'local_normalized':str(F/normalizer),
      'upper_normalized':str(upper/normalizer),
      'legacy_lower_normalized':str(legacy/normalizer),
      'strict_lower_normalized':str(strict/normalizer),
      'geometric_witness_valid': C>D(2) and delta<D(1)/3 and u<D(1),
      'valid':abs(residual)<D('1e-100')}

def signature(item):
    payload = {'version':VERSION, 'source_sha256':source_hash(), 'item':item,
       'precision':PREC, 'python':sys.version, 'platform':platform.platform()}
    return hashlib.sha256(json.dumps(payload,sort_keys=True).encode()).hexdigest()

def job(item, barrier=None):
    n, gamma, k = item; target = CACHE/f'n{n}_g{gamma.replace(".","p")}_k{k}.json'; sig = signature(item)
    if barrier is not None: barrier.wait(timeout=60)
    if target.exists():
        old = json.loads(target.read_text())
        if old.get('signature') == sig:
            old['cache'] = 'reused'; old['replay_pid'] = os.getpid(); return old
    out = data(n,gamma,k) | {'signature':sig, 'pid':os.getpid(), 'cache':'computed'}
    temporary = target.with_suffix('.tmp'); temporary.write_text(json.dumps(out)); temporary.replace(target)
    return out

def main():
    CACHE.mkdir(exist_ok=True); items=[(n,g,k) for n in NS for g in GS for k in KS]; rows=[]; failures=[]
    with Manager() as manager:
        barrier=manager.Barrier(15)
        with ProcessPoolExecutor(max_workers=15) as executor:
            first={executor.submit(job,x,barrier):x for x in items[:15]}
            for future in as_completed(first):
                try: rows.append(future.result())
                except Exception as exc: failures.append({'item':first[future],'error':repr(exc)})
            rest={executor.submit(job,x):x for x in items[15:]}
            for future in as_completed(rest):
                try: rows.append(future.result())
                except Exception as exc: failures.append({'item':rest[future],'error':repr(exc)})
    rows.sort(key=lambda row:(row['n'],D(row['gamma']),row['k']))
    result={'version':VERSION,'tasks':len(items),'precision':PREC,'rows':rows,'failures':failures,
            'worker_pids':sorted(set(row['pid'] for row in rows))}
    SUMMARY.write_text(json.dumps(result,indent=2))
    print(json.dumps({'tasks':len(rows),'failures':len(failures),'pids':result['worker_pids']}))
    if failures: raise SystemExit(1)

if __name__ == '__main__': main()
