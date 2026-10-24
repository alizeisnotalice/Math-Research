"""Closed form of L_alpha B for the centred cubic B-spline B(x) = (1/12) sum_k w_k |x-k|^3, w = (1,-4,6,-4,1), k=-2..2.
L_alpha |x|^p = gamma(p,alpha) |x|^{p-alpha} (regularised), gamma(3,alpha) = -24/(alpha(1-alpha)(2-alpha)(3-alpha)) * ... so
   alpha != 1: L_alpha B(s) = kappa_alpha * sum_k w_k |s-k|^{3-alpha},  kappa_alpha = -2/(alpha(1-alpha)(2-alpha)(3-alpha));
   alpha == 1: L_1 B(s) = - sum_k w_k (s-k)^2 log|s-k|.
(Checked against bspl.LB to ~1e-12.)"""
import numpy as np

W = np.array([1.0, -4.0, 6.0, -4.0, 1.0])
KS = np.arange(-2, 3)


def LBv(s, alpha):
    s = np.asarray(s, dtype=float)
    d = np.abs(s[..., None] - KS)
    if abs(alpha - 1.0) < 1e-14:
        with np.errstate(divide='ignore', invalid='ignore'):
            t = np.where(d > 0, d * d * np.log(np.where(d > 0, d, 1.0)), 0.0)
        return -(t * W).sum(-1)
    kap = -2.0 / (alpha * (1 - alpha) * (2 - alpha) * (3 - alpha))
    return kap * ((d ** (3 - alpha)) * W).sum(-1)


def Bv(s):
    s = np.abs(np.asarray(s, dtype=float))
    out = np.zeros_like(s)
    m1 = s < 1
    m2 = (s >= 1) & (s < 2)
    out[m1] = (4 - 6 * s[m1] ** 2 + 3 * s[m1] ** 3) / 6
    out[m2] = (2 - s[m2]) ** 3 / 6
    return out
