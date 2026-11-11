"""Multiprocess probe for P(a Z_c + b Z_cp + s G <= mu(a+b))."""
from __future__ import annotations
import math
import multiprocessing as mp
import numpy as np
from scipy.integrate import quad
from scipy.stats import exponnorm

LO, HI = 1.0, 2.0
for _ in range(100):
    mid = (LO + HI) / 2
    if mid < 1 + math.exp(-mid):
        LO = mid
    else:
        HI = mid
MU = (LO + HI) / 2


def hypo_cdf(x: float, a: float, b: float) -> float:
    if x <= 0:
        return 0.0
    if b == 0:
        return 1 - math.exp(-x/a)
    if abs(a - b) <= 1e-8 * max(a, b):
        q = x / ((a + b) / 2)
        return 1 - math.exp(-q) * (1 + q)
    return 1 - (a * math.exp(-x / a) - b * math.exp(-x / b)) / (a - b)


def exgauss_exp_cdf(x: float, scale: float, s: float) -> float:
    # Distribution of scale*Exp(1) + s*N(0,1), evaluated by SciPy's stable ex-Gaussian CDF.
    return float(exponnorm.cdf(x, K=scale / s, loc=0.0, scale=s))


def smooth_cdf(x: float, a: float, b: float, s: float) -> float:
    if s == 0:
        return hypo_cdf(x, a, b)
    if b == 0:
        return exgauss_exp_cdf(x, a, s)
    if a == b:
        scale = (a + b) / 2
        # Gamma(shape=2, scale) + Gaussian; smooth nonoscillatory integral.
        def f(t):
            return (t / scale**2) * math.exp(-t / scale) * 0.5 * math.erfc((t - x) / (s * math.sqrt(2)))
        return quad(f, 0, np.inf, epsabs=2e-12, epsrel=2e-12, limit=100)[0]
    # Hypoexponential density is a signed combination of two exponential densities.
    return (a * exgauss_exp_cdf(x, a, s) - b * exgauss_exp_cdf(x, b, s)) / (a - b)


def probability(r: float, c: float, cp: float, s: float) -> float:
    a, b = 1.0, r
    p, pp = (0.0 if c == math.inf else 1 / (1 + math.exp(c))), (0.0 if cp == math.inf else 1 / (1 + math.exp(cp)))
    total = 0.0
    for B in (0, 1):
        for D in (0, 1):
            weight = (p if B else 1-p) * (pp if D else 1-pp)
            if weight == 0:
                continue
            shift = (c if B else 0.0) + r * (cp if D else 0.0)
            total += weight * smooth_cdf(MU*(a+b)-shift, a, b, s)
    if not math.isfinite(total) or not -1e-12 <= total <= 1 + 1e-12:
        raise ArithmeticError((r, c, cp, s, total))
    return total


def worker(seed: int) -> list[tuple]:
    rng = np.random.default_rng(seed)
    rows = []
    # Distinct randomized batch per worker; r spans five decades, controls include large-c limits.
    for _ in range(1500):
        r = 10 ** rng.uniform(-5, 0)
        c, cp = rng.uniform(0, 15, 2)
        if rng.random() < .015: c = math.inf
        if rng.random() < .015: cp = math.inf
        s = 34 * rng.random() ** 2
        rows.append((probability(r,c,cp,s),r,c,cp,s,'random',seed))
    return rows


if __name__ == '__main__':
    # Compact deterministic boundary grid, evaluated once.
    fixed = []
    for r in [0.0,1e-5,1e-3,.01,.1,.5,1.0]:
        for c in [0,1,MU,4,15,math.inf]:
            for cp in [0,1,MU,4,15,math.inf]:
                for s in [0,.1,1,4,16,34]:
                    fixed.append((probability(r,c,cp,s),r,c,cp,s,'boundary-grid',0))
    with mp.Pool(4) as pool:
        rows = fixed + sum(pool.map(worker,[115201,115202,115203,115204]), [])
    print(f'mu={MU:.16g} workers=4 random=6000 fixed={len(fixed)}')
    for name, pred in [('all',lambda x:True),('s<=.5',lambda x:x[4]<=.5),('s<=2',lambda x:x[4]<=2),('s<=4',lambda x:x[4]<=4),('s==0',lambda x:x[4]==0),('r>=.1',lambda x:x[1]>=.1),('r>=.5',lambda x:x[1]>=.5)]:
        row = min(filter(pred,rows))
        print(name, f'{row[0]:.16g}', *row[1:])
