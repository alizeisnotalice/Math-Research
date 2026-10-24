"""Vectorised A_alpha of Psi = (r^-alpha - R^-alpha)_+  (numerical, for LP exploration).
Inside (r < R):  A Psi = C_alpha r^-2a + 2[ I(R-r) + I(R+u-v) ] + 2[ I(R-r) + I(R+v-u) ],
   I(d) = int_d^inf [R^-a - (s + R - d)^-a] s^{-1-a} ds = R^-a d^-a / a - d^-2a 2F1(a, 2a; 2a+1; -(R-d)/d) / (2a).
Outside (r > R): A Psi = 2 sum_j int Psi(x + t e_j) |t|^{-1-a} dt (Psi(x) = 0), by graded Gauss-Legendre.
"""
import numpy as np
from scipy.special import hyp2f1
from sing2 import C_alpha

_GX, _GW = np.polynomial.legendre.leggauss(40)


def I_d(alpha, R, d):
    d = np.asarray(d, float)
    z = -(R - d) / d
    return R ** (-alpha) * d ** (-alpha) / alpha - d ** (-2 * alpha) * hyp2f1(alpha, 2 * alpha, 2 * alpha + 1, z) / (2 * alpha)


def _gl(f, a, b):
    """vectorised Gauss-Legendre of f over [a, b] (arrays of equal shape)."""
    mid = 0.5 * (a + b); half = 0.5 * (b - a)
    t = mid[..., None] + half[..., None] * _GX
    return (f(t) * _GW).sum(-1) * half


def _outside_dir(alpha, R, u, v):
    """2 * int_{-(R-v)}^{R-v} [(|s|+v)^-a - R^-a] |s-u|^{-1-a} ds  for u > R - v (point outside), v < R."""
    L = R - v
    out = np.zeros_like(u)
    ok = L > 0
    u, v, L = u[ok], v[ok], L[ok]
    def f(s, u=None, v=None):
        return 0
    tot = np.zeros_like(u)
    # integrate over s in [0, L] with weight (u-s)^{-1-a} + (u+s)^{-1-a}; graded mesh near 0 (scale v) and near L
    # breakpoints: geometric from 0 up (scale v), and geometric toward L (scale u-L)
    nseg = 24
    gap = u - L
    for side in range(2):
        pass
    # build segments: [0, L] split at points v*2^k (k=-2..) and L - gap*2^k
    edges = [np.zeros_like(u)]
    for k in range(-3, 40):
        edges.append(np.minimum(L * 0.5, v * 2.0 ** k))
    for k in range(40, -4, -1):
        edges.append(np.maximum(L * 0.5, L - gap * 2.0 ** k))
    edges.append(L)
    E = np.sort(np.stack(edges, -1), -1)
    for j in range(E.shape[-1] - 1):
        a = E[..., j]; b = E[..., j + 1]
        nz = b > a
        if not nz.any():
            continue
        uu, vv = u[nz], v[nz]
        fn = lambda s: ((s + vv[:, None]) ** (-alpha) - R ** (-alpha)) * \
            ((uu[:, None] - s) ** (-1 - alpha) + (uu[:, None] + s) ** (-1 - alpha))
        tot[nz] += _gl(fn, a[nz], b[nz])
    out[ok] = 2 * tot
    return out


def A_psi_vec(alpha, R, u, v):
    u = np.abs(np.asarray(u, float)); v = np.abs(np.asarray(v, float))
    r = u + v
    out = np.empty_like(r)
    ins = r < R
    if ins.any():
        ui, vi, ri = u[ins], v[ins], r[ins]
        out[ins] = (C_alpha(alpha) * ri ** (-2 * alpha) + 4 * I_d(alpha, R, R - ri)
                    + 2 * I_d(alpha, R, R + ui - vi) + 2 * I_d(alpha, R, R + vi - ui))
    o = ~ins
    if o.any():
        uo, vo = u[o], v[o]
        out[o] = _outside_dir(alpha, R, uo, vo) + _outside_dir(alpha, R, vo, uo)
    on = np.abs(r - R) < 1e-12
    out[on] = np.inf
    return out
