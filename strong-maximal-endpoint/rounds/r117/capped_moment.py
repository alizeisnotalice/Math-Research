"""Integer certificate for all two-weight actions and all Gaussian variances."""
from concurrent.futures import ProcessPoolExecutor
import numpy as np

M, L, T, U = 1024, 512, 2048, 512
Q, ACC, J = 10**6, 10**15, 12000
N = 5200


def marginal(logistic):
    # Upper bound for (1000/1001)^j. The CDF is decreasing in this bound.
    a = ACC
    cdf = []
    for j in range(J + 1):
        cdf.append(Q * (ACC - a) // (ACC + a if logistic else ACC))
        a = (1000 * a + 1000) // 1001
    return np.diff(cdf, prepend=0).astype(np.int64)


def joint(mass, i):
    # Place both variables at right lattice endpoints. Omitted tails go to infinity.
    bins = (np.arange(J + 1) * (i + 1) + M - 1) // M
    second = np.zeros(bins[-1] + 1, dtype=np.int64)
    np.add.at(second, bins, mass)
    return np.cumsum(np.convolve(mass[:N+1], second[:N+1])[:N+1])


def batch(bounds):
    logistic, exponential = marginal(True), marginal(False)
    z = 2 * np.arange(U, dtype=np.int64) - U
    ell = np.arange(L // 2, dtype=np.int64)[:, None]
    factors = L * U + z * np.where(z < 0, ell + 1, ell)
    t = np.arange(T // 2, T, dtype=np.int64)
    small, large = None, None
    for i in range(*bounds):
        cdf = joint(logistic, i)
        pos = (1279 * (M + i) * factors) // (M * L * U)
        rows = cdf[pos].sum(axis=1) - U * Q**2 // 2
        j = int(rows.argmin())
        row = (int(rows[j]), i, j)
        small = min(small, row) if small else row

        cdf = joint(exponential, i)
        primitive = np.cumsum(np.concatenate(([0], cdf)))
        left = (1280 * (M+i+1) * (T-t) + M*T-1) // (M*T)
        right = ((M+i) * (1000*T + 1279*t)) // (M*T)
        rup = (1280 * (M+i+1) * (t+1) + M*T-1) // (M*T)
        rows = primitive[right] - primitive[left] - rup * Q**2
        j = int(rows.argmin())
        row = (int(rows[j]), i, int(t[j]))
        large = min(large, row) if large else row
    return small, large


if __name__ == '__main__':
    with ProcessPoolExecutor(max_workers=12) as pool:
        results = list(pool.map(batch, [(i, i+16) for i in range(0, M, 16)]))
    small, large = (min(r[k] for r in results) for k in (0, 1))
    print('uniform-noise probability surplus, numerator, b-bin, radius-bin:', small,
          'denominator:', U*Q**2)
    print('large-radius capped-moment surplus in log-lattice units:', large,
          'denominator:', Q**2)
    assert small[0] > 0 and large[0] > 0
